// simcent.cpp

#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simcent.hpp>

// =============================================================================
SimCent::SimCent(const YAML::Node& doc): SimParam(doc),
  cent_vars(MatrixXd::Zero(n_dim * n_robots, n_dim * n_robots))
{}

// =============================================================================
void SimCent::predict()
{
  SimParam::predict();

  // update diagonal matrices
  for (int i = 0; i < n_robots; ++i)
  {
    cent_vars.block(i * n_dim, i * n_dim, n_dim, n_dim) = vars[i];
  }
}

// =============================================================================
void SimCent::globalLocModeImpl(
    const MatrixXd& H, const MatrixXd& Q, const VectorXd& z_diff)
{
  MatrixXd cent_H = MatrixXd::Zero(2, n_robots*2);
  cent_H.block(0,0,2,2) = H;
  MatrixXd St = cent_H * cent_vars * cent_H.transpose() + Q;
  MatrixXd K = cent_vars * cent_H.transpose() * St.inverse();
  cent_vars
    = (MatrixXd::Identity(n_robots*n_dim, n_robots*n_dim) - K * cent_H)
    * cent_vars;
  VectorXd cent_means_diff = K * z_diff;
  for (int i = 0; i < n_robots; ++i)
  {
    means[i] += cent_means_diff.segment(i*n_dim, n_dim);
    vars[i] = cent_vars.block(i*n_dim,i*n_dim,n_dim,n_dim);
  }
}

// =============================================================================
void SimCent::mutualLocModeImpl(
    const VectorXd& z, const std::pair<int,int>& edge,
    const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
    const VectorXd& z_diff)
{
  // Matrices for centralized updates
  MatrixXd H;
  MatrixXd St;

  H = MatrixXd::Zero(2, n_robots*2);
  H.block(0, edge.first * 2, 2, 2)
    = H1;
  H.block(0, edge.second * 2, 2, 2)
    = H2;
  St = H * cent_vars * H.transpose() + Q;

  MatrixXd K = cent_vars * H.transpose() * St.inverse();
  VectorXd cent_means = VectorXd::Zero(n_robots * n_dim);
  cent_means.segment(edge.first*n_dim, n_dim)
    = means[edge.first];
  cent_means.segment(edge.second*n_dim, n_dim)
    = means[edge.second];
  cent_means += K * z_diff;
  means[edge.first]
    = cent_means.segment(edge.first*n_dim, n_dim);
  means[edge.second]
    = cent_means.segment(edge.second*n_dim, n_dim);
  cent_vars
    = (MatrixXd::Identity(n_robots*n_dim, n_robots*n_dim) - K*H)
    * cent_vars;
  vars[edge.first]
    = cent_vars.block(
        edge.first*n_dim, edge.first*n_dim, n_dim, n_dim);
  vars[edge.second]
    = cent_vars.block(
        edge.second*n_dim, edge.second*n_dim, n_dim, n_dim);
}
