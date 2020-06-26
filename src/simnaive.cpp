// simnaive.cpp

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simnaive.hpp>

using namespace std;
using namespace Eigen;

// =============================================================================
SimNaive::SimNaive(const YAML::Node& doc): SimParam(doc),
  mode1_rateQ(doc["mode1_rateQ"].as<double>())
{}

// =============================================================================
void SimNaive::mutualLocModeImpl(
  const VectorXd& z, const std::pair<int,int>& edge,
  const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
  const VectorXd& z_diff)
{
  // St matrices for decentralized updates
  MatrixXd St1;
  MatrixXd St2;

  St1 = H1 * vars[edge.first] * H1.transpose()
      + H2 * vars[edge.second] * H2.transpose() + (Q * mode1_rateQ);
  St2 = St1;
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

// =============================================================================
void SimNaiveOrig::mutualLocModeImpl(
  const VectorXd& z, const std::pair<int,int>& edge,
  const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
  const VectorXd& z_diff)
{
  // St matrices for decentralized updates
  MatrixXd St1;
  MatrixXd St2;

  St1 = H1 * vars[edge.first] * H1.transpose() + Q;
  St2 = H2 * vars[edge.second] * H2.transpose() + Q;
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
