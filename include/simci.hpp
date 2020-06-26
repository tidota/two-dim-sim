// simci.cpp

// this is a class implementing the covariance intersection.

#ifndef SIMCI_HPP
#define SIMCI_HPP

#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simparam.hpp>

using namespace std;
using namespace Eigen;

// helper function
void getOmega(MatrixXd &C1, MatrixXd &C2, double &omega, const int &mode);

class SimCi: public SimParam
{
  private: const bool mode4_original;
  private: const int mode4_omega_mode;
  private: const double mode4_rateQ;
  private: const double mode4_original_omega;
  private: const bool mode4_original_force;

  public: SimCi(const YAML::Node& doc);
  public: virtual ~SimCi(){}

  protected: virtual void mutualLocModeImpl(
    const VectorXd& z, const std::pair<int,int>& edge,
    const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
    const VectorXd& z_diff) override;
};

#endif
