// simpfcons.cpp

// this is a class implementing the conservative data exchange in a
// non-parametric way.

#ifndef SIMPFCONS_HPP
#define SIMPFCONS_HPP

#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simnonparam.hpp>

using namespace std;
using namespace Eigen;

class SimPfCons: public SimNonParam
{
  private: double mode6_omega;
  private: double mode6_sigma;

  public: SimPfCons(const YAML::Node& doc);
  public: virtual ~SimPfCons(){}

  private: void evalByOmega(
    const std::vector<VectorXd>& est, std::vector<double>& cumul_weights);
  private: void evalByZ(
    const std::vector<VectorXd>& est_target,
    const std::vector<double>& cumul_weights_target,
    const std::vector<VectorXd>& est_ref,
    const std::vector<double>& cumul_weights_ref,
    std::vector<double>& cumul_weights, const VectorXd& z);
  private: void resample(
    std::vector<VectorXd>& est, const std::vector<double>& cumul_weights);

  protected: virtual void mutualLocImpl(
    const VectorXd& z, const std::pair<int,int>& edge) override;
};

#endif
