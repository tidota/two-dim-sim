#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// for mode4 (Covariance Intersection)
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
