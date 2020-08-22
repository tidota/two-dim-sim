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
  const std::vector<VectorXd>& est, std::vector<double>& weights, double& w_sum)
{
  w_sum = 0;
  for (int i = 0; i < n_particles; ++i)
  {
    // estimate the value on the point of the probability distribution
    for (int j = 0; j < n_particles; ++j)
    {

    }

    // apply omega
    // this->mode6_omega

    // calculate weight
    // weights[i]
    //   = exp(
    //       -z_diff(0)*z_diff(0)/sigmaGlobalLocR/sigmaGlobalLocR
    //       -z_diff(1)*z_diff(1)/sigmaGlobalLocT/sigmaGlobalLocT);

    w_sum += weights[i];
  }
}

// =============================================================================
void SimPfCons::evalByZ(
  const std::vector<VectorXd>& est_target,
  const std::vector<double>& weights_target, const double& w_target_sum,
  const std::vector<VectorXd>& est_ref,
  const std::vector<double>& weights_ref, const double& w_ref_sum,
  std::vector<double>& weights, double& w_sum, const VectorXd& z)
{
  w_sum = 0;
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

    w_sum += weights[i];
  }
}

// =============================================================================
void SimPfCons::resample(
  std::vector<VectorXd>& est, const std::vector<double>& weights,
  const double& w_sum)
{
  // new  population
  std::vector<VectorXd> new_est;

  std::uniform_real_distribution<> dist(0, w_sum);
  for (int i = 0; i < n_particles; ++i)
  {
    // get a random number from a uniform distribution.
    double val = dist(gen);

    // decide the index to pick up
    int indx = 0;
    double buff = weights[0];
    while (indx < n_particles - 1 && val > buff)
    {
      ++indx;
      buff += weights[indx];
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

  // calcualte confidence reduction weights for each population
  std::vector<double> weights_omega1(0, n_particles);
  double weights_omega1_sum;
  evalByOmega(ests[r1], weights_omega1, weights_omega1_sum);

  std::vector<double> weights_omega2(0, n_particles);
  double weights_omega2_sum;
  evalByOmega(ests[r2], weights_omega2, weights_omega2_sum);

  if (enable_bidirectional)
  {
    // calculate mutual measurement weights for each population
    std::vector<double> weights1(0, n_particles);
    double weights1_sum;
    evalByZ(
      ests[r1], weights_omega1, weights_omega1_sum,
      ests[r2], weights_omega2, weights_omega2_sum,
      weights1, weights1_sum, -z);
    // resample
    resample(ests[r1], weights1, weights1_sum);
  }

  // calculate mutual measurement weights for each population
  std::vector<double> weights2(0, n_particles);
  double weights2_sum;
  evalByZ(
    ests[r2], weights_omega2, weights_omega2_sum,
    ests[r1], weights_omega1, weights_omega1_sum,
    weights2, weights2_sum, z);
  // resample
  resample(ests[r2], weights2, weights2_sum);
}
