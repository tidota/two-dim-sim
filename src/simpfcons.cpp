// simpfcons.cpp

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simpfcons.hpp>

using namespace std;
using namespace Eigen;

// =============================================================================
SimPfCons::SimPfCons(const YAML::Node& doc): SimNonParam(doc),
  mode3_omega_variable(doc["mode3_omega_variable"].as<bool>()),
  mode3_omega(doc["mode3_omega"].as<double>()),
  mode3_rateQ(doc["mode3_rateQ"].as<double>())
{}

// =============================================================================
void SimPfCons::mutualLocModeImpl(
  const VectorXd& z, const std::pair<int,int>& edge,
  const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
  const VectorXd& z_diff)
{
  // St matrices for decentralized updates
  MatrixXd St1;
  MatrixXd St2;

  if (mode3_omega_variable)
  {
    double sig1 = std::sqrt(vars[edge.first].determinant());
    double sig2 = std::sqrt(vars[edge.second].determinant());
    mode3_omega = 0.9 - std::fabs(sig1 - sig2)/(sig1 + sig2)*2*0.4;
  }
  St1 = H1 * (vars[edge.first] / mode3_omega) * H1.transpose()
      + H2 * (vars[edge.second] / (1 - mode3_omega)) * H2.transpose()
      + (Q * mode3_rateQ);
  St2 = H1 * (vars[edge.first] / (1 - mode3_omega)) * H1.transpose()
      + H2 * (vars[edge.second] / mode3_omega) * H2.transpose()
      + (Q * mode3_rateQ);
  vars[edge.first] /= mode3_omega;
  vars[edge.second] /= mode3_omega;

  if (enable_bidirectional)
  {
    MatrixXd K1 = vars[edge.first] * H1.transpose() * St1.inverse();
    means[edge.first] += K1 * z_diff;
    vars[edge.first]
      = (MatrixXd::Identity(n_dim, n_dim) - K1 * H1)
      * vars[edge.first];
  }

  MatrixXd K2 = vars[edge.second] * H2.transpose() * St2.inverse();
  means[edge.second] += K2 * z_diff;
  vars[edge.second]
    = (MatrixXd::Identity(n_dim, n_dim) - K2 * H2)
    * vars[edge.second];
}
