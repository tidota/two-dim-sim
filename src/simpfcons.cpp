// simpfcons.cpp

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simpfcons.hpp>

using namespace std;
using namespace Eigen;

// =============================================================================
SimPfCons::SimPfCons(const YAML::Node& doc): SimNonParam(doc),
  mode6_omega(doc["mode6_omega"].as<double>()),
  mode6_sigma(doc["mode6_sigma"].as<double>()),
  mode6_Nref(doc["mode6_Nref"].as<int>()),
  mode6_no_frac_exp(doc["mode6_no_frac_exp"].as<bool>()),
  mode6_rateQ(doc["mode6_rateQ"].as<double>())
{}

// =============================================================================
void SimPfCons::evalByOmega(
  const std::vector<VectorXd>& est, std::vector<double>& cumul_weights,
  std::vector<double>& cumul_weights_comp)
{
  for (int i = 0; i < n_particles; ++i)
  {
    cumul_weights[i] = 1.0; // exp(-0/mode6_sigma/mode6_sigma)
  }
  const double sigma_squared_x2 = mode6_sigma * mode6_sigma * 2;
  for (int i = 0; i < n_particles; ++i)
  {
    // estimate the value on the point of the probability distribution
    for (int j = i + 1; j < n_particles; ++j)
    {
      VectorXd diff = est[i] - est[j];
      double val = exp(-diff.squaredNorm()/sigma_squared_x2);
      cumul_weights[i] += val;
      cumul_weights[j] += val;
    }
    double buff = pow(cumul_weights[i], mode6_omega);
    // p^omega / p = p^(omega - 1)
    cumul_weights[i] = buff / cumul_weights[i];
    // p^(1-omega) / p = p^(-omega) = 1 / p^omega
    cumul_weights_comp[i] = 1.0 / buff;
    if (i > 0)
    {
      cumul_weights[i] += cumul_weights[i-1];
      cumul_weights_comp[i] += cumul_weights_comp[i-1];
    }
  }
}

// =============================================================================
void SimPfCons::evalByZ(
  const std::vector<VectorXd>& est_target,
  const std::vector<double>& cumul_weights_target,
  const std::vector<VectorXd>& est_ref,
  const std::vector<double>& cumul_weights_ref,
  std::vector<double>& cumul_weights, const VectorXd& z,
  const int est_ori_target, const int est_ori_ref)
{
  for (int i = 0; i < n_particles; ++i)
  {
    cumul_weights[i] = 0;
    for (int j = 0; j < mode6_Nref; ++j)
    {
      // draw a sample from the reference population
      int indx = drawRandIndx(cumul_weights_ref);

      VectorXd diff = est_target[i] - est_ref[indx];
      VectorXd z_hat(3);
      z_hat(0) = diff.norm();
      z_hat(1) = std::atan2(diff(1), diff(0));

      if (use_relative_bearing)
      {
        z_hat(1) -= est_oris[est_ori_ref][indx];
      }

      VectorXd z_diff = z - z_hat;
      if (z_diff(1) > M_PI)
      {
        z_diff(1) -= 2*M_PI;
      }
      else if (z_diff(1) <= -M_PI)
      {
        z_diff(1) += 2*M_PI;
      }

      // evaluate
      cumul_weights[i]
        += exp(
            -z_diff(0)*z_diff(0)/sigmaGlobalLocR/sigmaGlobalLocR*mode6_rateQ
            -z_diff(1)*z_diff(1)/sigmaGlobalLocT/sigmaGlobalLocT*mode6_rateQ);
    }
    double omega_weight = cumul_weights_target[i];
    if (i > 0)
    {
      omega_weight -= cumul_weights_target[i - 1];
    }
    cumul_weights[i] *= omega_weight;
    if (i > 0)
    {
      cumul_weights[i] += cumul_weights[i - 1];
    }
  }
}

// =============================================================================
void SimPfCons::resample(
  const int& robot_indx, const std::vector<double>& cumul_weights)
{
  // new  population
  std::vector<VectorXd> new_est(n_particles, VectorXd(n_dim));
  std::vector<double> new_est_ori(n_particles);

  for (int i = 0; i < n_particles; ++i)
  {
    // decide the index to pick up
    int indx = drawRandIndx(cumul_weights);

    // add the picked one
    new_est[i] = this->ests[robot_indx][indx];
    if (use_orientation)
      new_est_ori[i] = this->est_oris[robot_indx][indx];
  }

  // swap
  std::swap(this->ests[robot_indx], new_est);
  if (use_orientation)
    std::swap(this->est_oris[robot_indx], new_est_ori);
}

// =============================================================================
void SimPfCons::mutualLocImpl(const VectorXd& z, const std::pair<int,int>& edge)
{
  const int r1 = edge.first;
  const int r2 = edge.second;

  // calcualte cumulative confidence reduction weights for each population
  std::vector<double> cumul_weights_omega1(n_particles);
  std::vector<double> cumul_weights_omega1_comp(n_particles);
  if (mode6_no_frac_exp)
  {
    for (int i = 0; i < n_particles; ++i)
    {
      cumul_weights_omega1[i] = (i + 1) * (i + 2) / 2;
      cumul_weights_omega1_comp[i] = (i + 1) * (i + 2) / 2;
    }
  }
  else
  {
    evalByOmega(ests[r1], cumul_weights_omega1, cumul_weights_omega1_comp);
  }

  std::vector<double> cumul_weights_omega2(n_particles);
  std::vector<double> cumul_weights_omega2_comp(n_particles);
  if (mode6_no_frac_exp)
  {
    for (int i = 0; i < n_particles; ++i)
    {
      cumul_weights_omega2[i] = (i + 1) * (i + 2) / 2;
      cumul_weights_omega2_comp[i] = (i + 1) * (i + 2) / 2;
    }
  }
  else
  {
    evalByOmega(ests[r2], cumul_weights_omega2, cumul_weights_omega2_comp);
  }

  // calculate cumulative mutual measurement weights for each population
  std::vector<double> cumul_weights1(n_particles);
  if (enable_bidirectional)
  {
    VectorXd z_reversed = z;
    if (use_orientation)
    {
      z_reversed(1) = z_reversed(2);
      evalByZ(
        ests[r1], cumul_weights_omega1, ests[r2], cumul_weights_omega2_comp,
        cumul_weights1, z_reversed, r1, r2);
    }
    else
    {
      if (z_reversed(1) >= 0)
        z_reversed(1) -= M_PI;
      else
        z_reversed(1) += M_PI;
      evalByZ(
        ests[r1], cumul_weights_omega1, ests[r2], cumul_weights_omega2_comp,
        cumul_weights1, z_reversed);
    }
  }
  std::vector<double> cumul_weights2(n_particles);
  if (use_orientation)
  {
    evalByZ(
      ests[r2], cumul_weights_omega2, ests[r1], cumul_weights_omega1_comp,
      cumul_weights2, z, r2, r1);
  }
  else
  {
    evalByZ(
      ests[r2], cumul_weights_omega2, ests[r1], cumul_weights_omega1_comp,
      cumul_weights2, z);
  }

  // resample
  if (enable_bidirectional)
  {
    resample(r1, cumul_weights1);
  }
  resample(r2, cumul_weights2);
}
