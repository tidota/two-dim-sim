// simcons.cpp

// this is a class implementing the conservative data exchange.

#ifndef SIMCONS_HPP
#define SIMCONS_HPP

#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simparam.hpp>

using namespace std;
using namespace Eigen;

class SimCons: public SimParam
{
  private: const double mode3_omega_variable;
  private: double mode3_omega;
  private: const double mode3_rateQ;

  public: SimCons(const YAML::Node& doc);
  public: virtual ~SimCons(){}

  protected: virtual void mutualLocModeImpl(
    const VectorXd& z, const std::pair<int,int>& edge,
    const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
    const VectorXd& z_diff) override;
};

// this one is the original design where a square root is applied. That is a
// special case of omega = 0.5
class SimConsOrig: public SimParam
{
  private: const double mode2_rate1;
  private: const double mode2_rate2;
  private: const double mode2_rateQ;

  public: SimConsOrig(const YAML::Node& doc);
  public: virtual ~SimConsOrig(){}

  protected: virtual void mutualLocModeImpl(
    const VectorXd& z, const std::pair<int,int>& edge,
    const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
    const VectorXd& z_diff) override;
};

#endif
