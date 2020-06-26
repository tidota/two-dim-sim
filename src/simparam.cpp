// simparam.cpp

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simparam.hpp>

// =============================================================================
SimParam::SimParam(const YAML::Node& doc): SimBase(doc),
  alpha1M(doc["alpha1M"].as<double>()),
  alpha2M(doc["alpha2M"].as<double>()),
  betaM(doc["betaM"].as<double>()),
  last_loc(n_robots),
  last_mean(n_robots),
  last_var(n_robots)
{
  // for all robots
  for(int i = 0; i < n_robots; ++i)
  {
    // init location x
    VectorXd buff(n_dim);
    for (int j = 0; j < n_dim; ++j)
    {
      buff(j) = doc["robots"][i][j].as<double>();
    }
    means.push_back(buff);
    vars.push_back(MatrixXd::Zero(n_dim, n_dim));
  }
}

// =============================================================================
void SimParam::printSimInfo()
{
  SimBase::printSimInfo();

  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "means[" << i << "]: ";
    for (int idim = 0; idim < 2; ++idim)
      std::cout << means[i](idim) << ((idim == (2 - 1))? "": ", ");
    std::cout << "| ";
    std::cout << "vars[" << i << "]: " << vars[i].determinant() << std::endl;
  }
}

// =============================================================================
bool SimParam::startLog(const std::string& fname)
{
  // start log from parametric
  std::cout << "    t   |";
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "R[" << setw(7 * n_dim - 7) << i << "].x |";
    std::cout << "R[" << setw(7 * n_dim - 7) << i << "].m |";
    std::cout << " det(v) |";
    std::cout << " err |";
  }
  std::cout << std::endl;

  // prep output file
  if (!fout)
    fout.close();
  fout.open(fname);

  return fout.good();
}

// =============================================================================
void SimParam::endLog()
{
  fout.close();

  // output of gnuplot command
  const int off_next_robot = n_dim + n_dim + n_dim*n_dim + 2;
  std::cout << std::endl;
  std::cout << "~~~ gnuplot command (errors vs determinants) ~~~" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "--- ROBOT " << (i + 1) << " ---" << std::endl;
    std::cout << "clear" << std::endl;
    std::cout << "unset object" << std::endl;
    std::cout << "plot \"output.dat\" u 1:"
              << std::to_string(2+(i+1)*off_next_robot-2)
              << " title \"|Sigma|^0.5 of robot" << std::to_string(1+i) << "\" with line, ";
    std::cout << "\"output.dat\" u 1:"
              << std::to_string(2+(i+1)*off_next_robot-1)
              << " title \"err of robot" << std::to_string(1+i) << "\" with line";
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "~~~ gnuplot command (errors as a whole) ~~~" << std::endl;
  std::cout << "clear" << std::endl;
  std::cout << "unset object" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    if (i == 0)
      std::cout << "plot ";
    else
      std::cout << "     ";
    std::cout << "\"output.dat\" u 1:"
              << std::to_string(2+(i+1)*off_next_robot-1)
              << " title \"err" << std::to_string(1+i) << "\" with line";
    if (i < n_robots - 1)
      std::cout << ", \\";
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "~~~ gnuplot command (locations) ~~~" << std::endl;
  std::cout << "clear" << std::endl;
  std::cout << "unset object" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    Eigen::EigenSolver<MatrixXd> s(vars[i]);
    auto eigen_val = s.eigenvalues();
    auto eigen_vec = s.eigenvectors();
    double var_ang
      = std::atan2(eigen_vec.col(0)[1].real(), eigen_vec.col(0)[0].real())
        / M_PI*180.0;
    std::cout << "set object " << std::to_string(i + 1)
              << " ellipse center "
              << std::to_string(last_mean[i](0)) << ","
              << std::to_string(last_mean[i](1))
              << " size "
              << std::to_string(std::sqrt(eigen_val[0].real()) * 2) << ","
              << std::to_string(std::sqrt(eigen_val[1].real()) * 2)
              << " angle "
              << std::to_string(var_ang)
              << " front fillstyle empty border -1" << std::endl;
  }
  std::cout << "h1 = 227/360.0" << std::endl;
  std::cout << "h2 = 40/360.0" << std::endl;
  std::cout << "set palette model HSV functions (1-gray)*(h2-h1)+h1,1,0.68"
            << std::endl;
  std::cout << "set size ratio -1" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    if (i == 0)
      std::cout << "plot   ";
    else
      std::cout << "replot ";
    std::cout << "\"output.dat\" u ";
    for (int j = 0; j < n_dim; ++j)
    {
      std::cout << std::to_string(2+i*off_next_robot+j);
      if (j < n_dim - 1)
        std::cout << ":";
      else
        std::cout << ":1";
    }
    std::cout << " title \"R" << std::to_string(1+i) << "\""
              << " with linespoints lt -1 lw 1.0 ps 3.0"
              << " pt " << std::to_string(i + 1)
              << " lc palette" << std::endl;
  }
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "replot ";
    std::cout << "\"output.dat\" u ";
    for (int j = 0; j < n_dim; ++j)
    {
      std::cout << std::to_string(2+n_dim+i*off_next_robot+j);
      if (j < n_dim - 1)
        std::cout << ":";
      else
        std::cout << ":1";
    }
    std::cout << " title \"Est" << std::to_string(1+i) << "\""
              << " with linespoints lt 1 lw 3.0 ps 3.0"
              << " pt " << std::to_string(i + 1)
              << " lc palette" << std::endl;
  }
  std::cout << std::endl;

  // display errors
  double total_error = 0;
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "robot[" << i << "]'s average error:" << (errors[i]/(max_time*sim_freq)) << std::endl;
    total_error += errors[i];
  }
  std::cout << "overall average error: " << (total_error/(max_time*sim_freq)/n_robots) << std::endl;
}

// =============================================================================
void SimParam::plotImpl()
{
  std::cout << std::fixed << std::setprecision(2);
  std::cout << std::right << std::setw(8) << t << "|";
  fout << std::fixed << std::setprecision(3);
  fout << std::right << std::setw(8) << t << " ";

  for (int i = 0; i < n_robots; ++i)
  {
    for (int j = 0; j < n_dim; ++j)
    {
      std::cout << std::right << std::setw(6) << robots[i](j);
      if (j < n_dim - 1)
        std::cout << ",";
      else
        std::cout << "|";
    }
    for (int j = 0; j < n_dim; ++j)
    {
      std::cout << std::right << std::setw(6) << means[i](j);
      if (j < n_dim - 1)
        std::cout << ",";
      else
        std::cout << "|";
    }
    std::cout << std::right << std::setw(8) << std::sqrt(vars[i].determinant()) << "|";
    std::cout << std::right << std::setw(5) << (robots[i] - means[i]).norm() << "|";

    for (int j = 0; j < n_dim; ++j)
    {
      fout << std::right << std::setw(9) << robots[i](j);
      fout << " ";
    }
    for (int j = 0; j < n_dim; ++j)
    {
      fout << std::right << std::setw(9) << means[i](j);
      fout << " ";
    }
    for (int j = 0; j < n_dim; ++j)
    {
      for (int k = 0; k < n_dim; ++k)
      {
        fout << std::right << std::setw(9) << vars[i](j,k);
        fout << " ";
      }
    }
    fout << std::setw(8) << std::sqrt(vars[i].determinant()) << " ";
    fout << std::right << std::setw(8)
         << (robots[i] - means[i]).norm() << " ";

    // store data to the buffer to plot
    last_loc[i] = robots[i];
    last_mean[i] = means[i];
    last_var[i] = vars[i];
  }

  std::cout << std::endl;
  fout << std::endl;
}

// =============================================================================
void SimParam::predict()
{
  // === prediction ===
  // for all robots
  for (int i = 0; i < n_robots; ++i)
  {
    // motion model
    means[i] = means[i] + MatrixXd::Identity(n_dim, n_dim) * vels[i] * deltaT;

    MatrixXd FloatEffect = MatrixXd::Identity(2, 2) * betaM;
    double v = vels[i].norm();
    MatrixXd EigenVecs(2, 2);
    EigenVecs(0, 0) = vels[i](0) / v;
    EigenVecs(1, 0) = vels[i](1) / v;
    EigenVecs(0, 1) = -vels[i](1) / v;
    EigenVecs(1, 1) = vels[i](0) / v;
    MatrixXd EigenVals = MatrixXd::Identity(2, 2);
    EigenVals(0, 0) = alpha1M * v * v;
    EigenVals(1, 1) = alpha2M * v * v;
    // transpose == inverse because eigen vectors are perpendicular to each other
    MatrixXd M = EigenVecs * EigenVals * EigenVecs.transpose() + FloatEffect;
    MatrixXd V = MatrixXd::Identity(n_dim, n_dim) * deltaT;
    vars[i] = vars[i] + V.transpose() * M * V;
  }
}

// =============================================================================
void SimParam::globalLocImpl(const VectorXd& z)
{
  // update its location estimate
  VectorXd z_hat(2);
  z_hat(0) = means[0].norm();
  z_hat(1) = std::atan2(means[0](1), means[0](0));
  double q = means[0].squaredNorm();
  MatrixXd H(2,n_dim);
  H(0, 0) = means[0](0)/std::sqrt(q);
  H(1, 0) = -means[0](1)/q;
  H(0, 1) = means[0](1)/std::sqrt(q);
  H(1, 1) = means[0](0)/q;
  MatrixXd Q = MatrixXd::Zero(2,2);
  Q(0,0) = sigmaGlobalLocR * sigmaGlobalLocR;
  Q(1,1) = sigmaGlobalLocT * sigmaGlobalLocT;
  VectorXd z_diff = z - z_hat;
  if (z_diff(1) > M_PI)
  {
    z_diff(1) -= 2*M_PI;
  }
  else if (z_diff(1) <= -M_PI)
  {
    z_diff(1) += 2*M_PI;
  }

  globalLocModeImpl(H, Q, z_diff);
}

// =============================================================================
void SimParam::globalLocModeImpl(
  const MatrixXd& H, const MatrixXd& Q, const VectorXd& z_diff)
{
  MatrixXd St = H * vars[0] * H.transpose() + Q;
  MatrixXd K = vars[0] * H.transpose() * St.inverse();
  means[0] += K * z_diff;
  vars[0] = (MatrixXd::Identity(n_dim, n_dim) - K * H) * vars[0];
}

// =============================================================================
void SimParam::mutualLocImpl(
  const VectorXd& z, const std::pair<int,int>& edge)
{
  // update its location estimate
  VectorXd diff_hat = means[edge.second] - means[edge.first];
  VectorXd z_hat(2);
  z_hat(0) = diff_hat.norm();
  z_hat(1) = std::atan2(diff_hat(1), diff_hat(0));
  double q = diff_hat.squaredNorm();

  MatrixXd Q = MatrixXd::Zero(2,2);
  Q(0,0) = sigmaS_R * sigmaS_R;
  Q(1,1) = sigmaS_T * sigmaS_T;

  MatrixXd H1(2,n_dim);
  H1(0, 0) = -diff_hat(0)/std::sqrt(q);
  H1(1, 0) = diff_hat(1)/q;
  H1(0, 1) = -diff_hat(1)/std::sqrt(q);
  H1(1, 1) = -diff_hat(0)/q;

  MatrixXd H2(2,n_dim);
  H2(0, 0) = diff_hat(0)/std::sqrt(q);
  H2(1, 0) = -diff_hat(1)/q;
  H2(0, 1) = diff_hat(1)/std::sqrt(q);
  H2(1, 1) = diff_hat(0)/q;

  VectorXd z_diff = z - z_hat;
  if (z_diff(1) > M_PI)
  {
    z_diff(1) -= 2*M_PI;
  }
  else if (z_diff(1) <= -M_PI)
  {
    z_diff(1) += 2*M_PI;
  }

  mutualLocModeImpl(z, edge, H1, H2, Q, z_diff);
}

// =============================================================================
void SimParam::calcErrors()
{
  // calculate errors
  for (int i = 0; i < n_robots; ++i)
  {
    errors[i] += (means[i] - robots[i]).norm();
  }
}
