// simpfcons.cpp

// this is a class implementing the conservative data exchange in a
// non-parametric way.

#ifndef SIMPFCONS_HPP
#define SIMPFCONS_HPP

#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simparam.hpp>

using namespace std;
using namespace Eigen;

class SimPfCons: public SimNonParam
{
  private: const double mode3_omega_variable;
  private: double mode3_omega;
  private: const double mode3_rateQ;

  public: SimPfCons(const YAML::Node& doc);
  public: virtual ~SimPfCons(){}

  protected: virtual void mutualLocModeImpl(
    const VectorXd& z, const std::pair<int,int>& edge,
    const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
    const VectorXd& z_diff) override;
};

#endif
