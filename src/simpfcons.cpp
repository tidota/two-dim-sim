// simpfcons.cpp

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simpfcons.hpp>

using namespace std;
using namespace Eigen;

// =============================================================================
SimPfCons::SimPfCons(const YAML::Node& doc): SimNonParam(doc),
  mode6_omega(doc["mode6_omega"].as<double>()),
  mode6_sigma(doc["mode6_sigma"].as<double>())
{}

// =============================================================================
void SimPfCons::evalByOmega(
  const std::vector<VectorXd>& est, std::vector<double>& cumul_weights)
{
  for (int i = 0; i < n_particles; ++i)
    cumul_weights[i] = 1.0; // exp(-0/mode6_sigma/mode6_sigma)
  for (int i = 0; i < n_particles; ++i)
  {
    // estimate the value on the point of the probability distribution
    for (int j = i + 1; j < n_particles; ++j)
    {
      VectorXd diff = est[i] - est[j];
      double val = exp(-diff.squaredNorm()/mode6_sigma/mode6_sigma);
      cumul_weights[i] += val;
      cumul_weights[j] += val;
    }
    cumul_weights[i] = pow(cumul_weights[i], mode6_omega - 1);
    if (i > 0)
      cumul_weights[i] += cumul_weights[i-1];
  }
}

// =============================================================================
void SimPfCons::evalByZ(
  const std::vector<VectorXd>& est_target,
  const std::vector<double>& cumul_weights_target,
  const std::vector<VectorXd>& est_ref,
  const std::vector<double>& cumul_weights_ref,
  std::vector<double>& cumul_weights, const VectorXd& z)
{
  for (int i = 0; i < n_particles; ++i)
  {
    // draw a sample from the reference population
    int indx = drawRandIndx(cumul_weights_ref);

    VectorXd diff = est_target[i] - est_ref[indx];
    VectorXd z_hat(2);
    z_hat(0) = diff.norm();
    z_hat(1) = std::atan2(diff(1), diff(0));

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
      = exp(
          -z_diff(0)*z_diff(0)/sigmaGlobalLocR/sigmaGlobalLocR
          -z_diff(1)*z_diff(1)/sigmaGlobalLocT/sigmaGlobalLocT);
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
  std::vector<VectorXd>& est, const std::vector<double>& cumul_weights)
{
  // new  population
  std::vector<VectorXd> new_est(n_particles, VectorXd(n_dim));

  for (int i = 0; i < n_particles; ++i)
  {
    // decide the index to pick up
    int indx = drawRandIndx(cumul_weights);

    // add the picked one
    new_est[i] = est[indx];
  }

  // swap
  std::swap(est, new_est);
}

// =============================================================================
void SimPfCons::mutualLocImpl(const VectorXd& z, const std::pair<int,int>& edge)
{
  const int r1 = edge.first;
  const int r2 = edge.second;

  // calcualte cumulative confidence reduction weights for each population
  std::vector<double> cumul_weights_omega1(n_particles);
  evalByOmega(ests[r1], cumul_weights_omega1);

  std::vector<double> cumul_weights_omega2(n_particles);
  evalByOmega(ests[r2], cumul_weights_omega2);

  if (enable_bidirectional)
  {
    // calculate cumulative mutual measurement weights for each population
    std::vector<double> cumul_weights1(n_particles);
    evalByZ(
      ests[r1], cumul_weights_omega1, ests[r2], cumul_weights_omega2,
      cumul_weights1, -z);
    // resample
    resample(ests[r1], cumul_weights1);
  }

  // calculate mutual measurement weights for each population
  std::vector<double> cumul_weights2(n_particles);
  evalByZ(
    ests[r2], cumul_weights_omega2, ests[r1], cumul_weights_omega1,
    cumul_weights2, z);
  // resample
  resample(ests[r2], cumul_weights2);
}
