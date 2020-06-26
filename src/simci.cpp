// simci.cpp

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simci.hpp>

using namespace std;
using namespace Eigen;

// =============================================================================
SimCi::SimCi(const YAML::Node& doc): SimParam(doc),
  mode4_original(doc["mode4_original"].as<bool>()),
  mode4_omega_mode(doc["mode4_omega_mode"].as<int>()),
  mode4_rateQ(doc["mode4_rateQ"].as<double>()),
  mode4_original_omega(doc["mode4_original_omega"].as<double>()),
  mode4_original_force(doc["mode4_original_force"].as<bool>())
{}

// =============================================================================
void SimCi::mutualLocModeImpl(
    const VectorXd& z, const std::pair<int,int>& edge,
    const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
    const VectorXd& z_diff)
{
  // St matrices for decentralized updates
  MatrixXd St1;
  MatrixXd St2;

  double omega;
  VectorXd mean1 = means[edge.first];
  MatrixXd var1 = vars[edge.first];
  VectorXd mean2 = means[edge.second];
  MatrixXd var2 = vars[edge.second];
  if (mode4_original)
  {
    VectorXd diff_hat = means[edge.second] - means[edge.first];
    double q = diff_hat.squaredNorm();
    MatrixXd H1_mutable = H1;
    MatrixXd H2_mutable = H2;
    VectorXd z_diff_mutable = z_diff;

    // change dz/dx to dx/dz
    // assuming the original jacobians were already calculated
    //H1 = H1.transpose().eval();
    H1_mutable(1, 0) = H1_mutable(1, 0) * q;
    H1_mutable(1, 1) = H1_mutable(1, 1) * q;
    if (var1.determinant() != 0)
    {
      MatrixXd A = vars[edge.first];
      MatrixXd B
        = vars[edge.second]
        + H1_mutable.transpose() * (Q * mode4_rateQ) * H1_mutable;
      VectorXd a = means[edge.first];
      VectorXd b = means[edge.second];
      b(0) = b(0) - z(0) * std::cos(z(1));
      b(1) = b(1) - z(0) * std::sin(z(1));

      if (mode4_original_force)
        omega = mode4_original_omega;
      else
        getOmega(A, B, omega, mode4_omega_mode);
      MatrixXd C = (omega*A.inverse() + (1-omega)*B.inverse()).inverse();
      VectorXd c
        = C * (omega*A.inverse()*a + (1-omega)*B.inverse()*b);
      mean1 = c;
      var1 = C;
    }
    // change dz/dx to dx/dz
    // assuming the original jacobians were already calculated
    //H2 = H2.transpose().eval();
    H2_mutable(1, 0) = H2_mutable(1, 0) * q;
    H2_mutable(1, 1) = H2_mutable(1, 1) * q;
    if (vars[edge.second].determinant() != 0)
    {
      MatrixXd A = vars[edge.second];
      MatrixXd B
        = vars[edge.first]
        + H2_mutable.transpose() * (Q * mode4_rateQ) * H2_mutable;
      VectorXd a = means[edge.second];
      VectorXd b = means[edge.first];
      b(0) = b(0) + z(0) * std::cos(z(1));
      b(1) = b(1) + z(0) * std::sin(z(1));

      if (mode4_original_force)
        omega = mode4_original_omega;
      else
        getOmega(A, B, omega, mode4_omega_mode);
      MatrixXd C = (omega*A.inverse() + (1-omega)*B.inverse()).inverse();
      VectorXd c
        = C * (omega*A.inverse()*a + (1-omega)*B.inverse()*b);
      mean2 = c;
      var2 = C;
    }
  }
  else
  {
    {
      MatrixXd C1 = H1 * vars[edge.first] * H1.transpose();
      MatrixXd C2
        = H2 * vars[edge.second] * H2.transpose()
        + (Q * mode4_rateQ);
      getOmega(C1, C2, omega, mode4_omega_mode);
    }
    if (omega == 0)
    {
      VectorXd corr(n_dim);
      corr(0) = -z(0) * std::cos(z(1));
      corr(1) = -z(0) * std::sin(z(1));
      mean1
        = means[edge.second]
        + corr;
      var1 = vars[edge.second]
           + (H1.transpose() * H1).inverse() * H1.transpose()
           * (Q * mode4_rateQ)
           * H1 * (H1.transpose() * H1).inverse();
    }
    else if (omega < 1)
    {
      St1 = H1 * (vars[edge.first]/omega) * H1.transpose()
          + H2 * (vars[edge.second]/(1-omega)) * H2.transpose()
          + ((Q * mode4_rateQ)/(1-omega));
      MatrixXd K1
        = (vars[edge.first]/omega) * H1.transpose() * St1.inverse();
      mean1 += K1 * z_diff;
      var1
        = (MatrixXd::Identity(n_dim, n_dim) - K1 * H1)
        * (vars[edge.first]/omega);
    }

    {
      MatrixXd C1 = H2 * vars[edge.second] * H2.transpose();
      MatrixXd C2
        = H1 * vars[edge.first] * H1.transpose()
        + (Q * mode4_rateQ);
      getOmega(C1, C2, omega, mode4_omega_mode);
    }
    if (omega == 0)
    {
      VectorXd corr(n_dim);
      corr(0) = z(0) * std::cos(z(1));
      corr(1) = z(0) * std::sin(z(1));
      mean2
        = means[edge.first]
        + corr;
      var2 = vars[edge.first]
           + (H2.transpose() * H2).inverse() * H2.transpose()
           * (Q * mode4_rateQ)
           * H2 * (H2.transpose() * H2).inverse();
    }
    else if (omega < 1)
    {
      St2 = H1 * (vars[edge.first]/(1-omega)) * H1.transpose()
          + H2 * (vars[edge.second]/omega) * H2.transpose()
          + ((Q * mode4_rateQ)/(1-omega));
      MatrixXd K2
        = (vars[edge.second]/omega) * H2.transpose() * St2.inverse();
      mean2 += K2 * z_diff;
      var2
        = (MatrixXd::Identity(n_dim, n_dim) - K2 * H2)
        * (vars[edge.second]/omega);
    }
  }
  means[edge.first] = mean1;
  vars[edge.first] = var1;
  means[edge.second] = mean2;
  vars[edge.second] = var2;
}

// =============================================================================
// helper function for mode4 (Covariance Intersection)
void getOmega(MatrixXd &C1, MatrixXd &C2, double &omega, const int &mode)
{
  if (mode == 0)
  {
    omega = 0.5;
    return;
  }
  double det1 = C1.determinant();
  double det2 = C2.determinant();
  if (det1 == 0 || std::isnan(det1))
  {
    omega = 1;
    return;
  }
  if (det2 == 0 || std::isnan(det2))
  {
    omega = 0;
    return;
  }
  if (mode == 1)
  {
    if (det1 < det2)
      omega = 1;
    else
      omega = 0;
    return;
  }
  if (mode == 2)
  {
    if (det1 <= det2)
      omega = 0.6;
    else
      omega = 0.4;
    return;
  }
  if (mode == 3)
  {
    double sig1 = std::sqrt(det1);
    double sig2 = std::sqrt(det2);
    omega = sig2 / (sig1 + sig2);
    return;
  }
  Eigen::EigenSolver<MatrixXd> slv1(C1);
  auto eigen_values1 = slv1.eigenvalues();
  auto eigen_vector1 = slv1.eigenvectors();
  MatrixXd D1 = MatrixXd::Zero(2,2);
  MatrixXd V1 = MatrixXd::Zero(2,2);
  D1(0, 0) = std::sqrt(eigen_values1[0].real());
  D1(1, 1) = std::sqrt(eigen_values1[1].real());
  V1(0, 0) = eigen_vector1.col(0)[0].real();
  V1(0, 1) = eigen_vector1.col(1)[0].real();
  V1(1, 0) = eigen_vector1.col(0)[1].real();
  V1(1, 1) = eigen_vector1.col(1)[1].real();
  MatrixXd T1 = D1.inverse() * V1.transpose();

  Eigen::EigenSolver<MatrixXd> slv2(T1*C2*T1.transpose());
  auto eigen_values2 = slv2.eigenvalues();
  auto eigen_vector2 = slv2.eigenvectors();
  MatrixXd D2 = MatrixXd::Zero(2,2);
  MatrixXd V2 = MatrixXd::Zero(2,2);
  // D2(0, 0) = std::sqrt(eigen_values2[0].real());
  // D2(1, 1) = std::sqrt(eigen_values2[1].real());
  V2(0, 0) = eigen_vector2.col(0)[0].real();
  V2(0, 1) = eigen_vector2.col(1)[0].real();
  V2(1, 0) = eigen_vector2.col(0)[1].real();
  V2(1, 1) = eigen_vector2.col(1)[1].real();

  double d1 = 1.0/eigen_values2[0].real();
  d1 = d1/(1-d1);
  double d2 = 1.0/eigen_values2[1].real();
  d2 = d2/(1-d2);

  MatrixXd diag = V2*D1*V2.transpose();
  double a1 = diag(0, 0);
  double a2 = diag(1, 1);

  double p = (a1*d2*(1+d1) + a2*d1*(1+d2))
           / (a1*(1+d1) + a2*(1+d2));
  double q = (a1*d2*d2*(1+d1) + a2*d1*d1*(1+d2))
           / (a1*(1+d1) + a2*(1+d2));
  double omega1 = -p + std::sqrt(p*p - q);
  double omega2 = -p - std::sqrt(p*p - q);

  omega = omega1;

  if (omega < 0 || 1 < omega || std::isnan(omega))
    omega = omega2;

  if (omega < 0 || 1 < omega || std::isnan(omega))
  {
    if (det1 < det2)
      omega = 1;
    else
      omega = 0;
  }
}
