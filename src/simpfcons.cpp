// simpfcons.cpp

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simpfcons.hpp>

using namespace std;
using namespace Eigen;

// =============================================================================
SimPfCons::SimPfCons(const YAML::Node& doc): SimNonParam(doc),
  mode6_omega_variable(doc["mode6_omega_variable"].as<bool>()),
  mode6_omega(doc["mode6_omega"].as<double>())
{}

// =============================================================================
void SimPfCons::mutualLocImpl(const VectorXd& z, const std::pair<int,int>& edge)
{
  const int r1 = edge.first;
  const int r2 = edge.second;

  // TODO
  // calcualte confidence reduction weights for each population
  // calculate mutual measurement weights for each population
  // resampling

  if (enable_bidirectional)
  {
  }

  // create weights for resampling
  std::vector<double> weights(0, n_particles);
  double w_sum = 0;

  // for all particles of the first robot
  for (int i = 0; i < n_particles; ++i)
  {
    VectorXd z_hat(2);
    z_hat(0) = ests[0][i].norm();
    z_hat(1) = std::atan2(ests[0][i](1), ests[0][i](0));

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
    weights[i]
      = exp(
          -z_diff(0)*z_diff(0)/sigmaGlobalLocR/sigmaGlobalLocR
          -z_diff(1)*z_diff(1)/sigmaGlobalLocT/sigmaGlobalLocT);

    w_sum += weights[i];
  }

  // new  population
  std::vector<VectorXd> new_est;

  // resample
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
    new_est.push_back(this->ests[0][indx]);
  }

  // swap
  std::swap(this->ests[0], new_est);
}
