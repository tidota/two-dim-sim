// simnaive.cpp

// this is a class implementing the naive data exchange.

#ifndef SIMNAIVE_HPP
#define SIMNAIVE_HPP

#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simparam.hpp>

using namespace std;
using namespace Eigen;

class SimNaive: public SimParam
{
  protected: const double mode1_rateQ;

  public: SimNaive(const YAML::Node& doc);
  public: ~SimNaive(){}

  protected: virtual void mutualLocModeImpl(
    const VectorXd& z, const std::pair<int,int>& edge,
    const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
    const VectorXd& z_diff) override;
};

// This is the original version of the naive method
class SimNaiveOrig: public SimParam
{
  public: SimNaiveOrig(const YAML::Node& doc): SimParam(doc){}
  public: ~SimNaiveOrig(){}

  protected: virtual void mutualLocModeImpl(
    const VectorXd& z, const std::pair<int,int>& edge,
    const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
    const VectorXd& z_diff) override;
};

#endif
