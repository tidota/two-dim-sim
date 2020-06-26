// simcent.hpp

// this is a class implementing the centralized localization.

#ifndef SIMCENT_HPP
#define SIMCENT_HPP

#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simparam.hpp>

using namespace std;
using namespace Eigen;

class SimCent: public SimParam
{
  // covariance matrix
  private: MatrixXd cent_vars;

  public: SimCent(const YAML::Node& doc);
  public: virtual ~SimCent(){}

  public: virtual void predict() override;
  protected: virtual void globalLocModeImpl(
    const MatrixXd& H, const MatrixXd& Q, const VectorXd& z_diff) override;
  protected: virtual void mutualLocModeImpl(
    const VectorXd& z, const std::pair<int,int>& edge,
    const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
    const VectorXd& z_diff) override;
};

#endif
