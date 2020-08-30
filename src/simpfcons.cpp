// simpfcons.cpp

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simpfcons.hpp>

using namespace std;
using namespace Eigen;

// =============================================================================
SimPfCons::SimPfCons(const YAML::Node& doc): SimNonParam(doc),
  mode6_omega(doc["mode6_omega"].as<double>())
{}

// =============================================================================
void SimPfCons::evalByOmega(
  const std::vector<VectorXd>& est, std::vector<double>& cumul_weights)
{
  for (int i = 0; i < n_particles; ++i)
  {
    // estimate the value on the point of the probability distribution
    for (int j = 0; j < n_particles; ++j)
    {

    }

    // apply omega
    // this->mode6_omega

    // calculate weight
    // cumul_weights[i]
    //   = exp(
    //       -z_diff(0)*z_diff(0)/sigmaGlobalLocR/sigmaGlobalLocR
    //       -z_diff(1)*z_diff(1)/sigmaGlobalLocT/sigmaGlobalLocT);
  }
}

// =============================================================================
void SimPfCons::evalByZ(
  const std::vector<VectorXd>& est_target,
  const std::vector<double>& cumul_weights_target,
  const std::vector<VectorXd>& est_ref,
  const std::vector<double>& cumul_weights_ref,
  std::vector<double>& weights, const VectorXd& z)
{
  for (int i = 0; i < n_particles; ++i)
  {
    for (int j = 0; j < n_particles; ++j)
    {
      VectorXd diff = est_target[i] - est_ref[j];
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
    }

    // evaluate
    // weights[i]
    //   = exp(
    //       -z_diff(0)*z_diff(0)/sigmaGlobalLocR/sigmaGlobalLocR
    //       -z_diff(1)*z_diff(1)/sigmaGlobalLocT/sigmaGlobalLocT);
  }
}

// =============================================================================
void SimPfCons::resample(
  std::vector<VectorXd>& est, const std::vector<double>& cumul_weights)
{
  // new  population
  std::vector<VectorXd> new_est;

  std::uniform_real_distribution<>
    dist(0, cumul_weights[cumul_weights.size() - 1]);
  for (int i = 0; i < n_particles; ++i)
  {
    // get a random number from a uniform distribution.
    double val = dist(gen_pf);

    // decide the index to pick up
    int indx = 0;
    double buff = cumul_weights[0];
    while (indx < n_particles - 1 && val > buff)
    {
      ++indx;
      buff += cumul_weights[indx];
    }

    // add the picked one
    new_est.push_back(est[indx]);
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
  std::vector<double> cumul_weights_omega1(n_particles, 0);
  evalByOmega(ests[r1], cumul_weights_omega1);

  std::vector<double> cumul_weights_omega2(n_particles, 0);
  evalByOmega(ests[r2], cumul_weights_omega2);

  if (enable_bidirectional)
  {
    // calculate cumulative mutual measurement weights for each population
    std::vector<double> cumul_weights1(n_particles, 0);
    evalByZ(
      ests[r1], cumul_weights_omega1, ests[r2], cumul_weights_omega2,
      cumul_weights1, -z);
    // resample
    resample(ests[r1], cumul_weights1);
  }

  // calculate mutual measurement weights for each population
  std::vector<double> cumul_weights2(n_particles, 0);
  evalByZ(
    ests[r2], cumul_weights_omega2, ests[r1], cumul_weights_omega1,
    cumul_weights2, z);
  // resample
  resample(ests[r2], cumul_weights2);
}
