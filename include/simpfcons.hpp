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
  private: double mode6_omega;

  public: SimPfCons(const YAML::Node& doc);
  public: virtual ~SimPfCons(){}

  protected: virtual void mutualLocImpl(
    const VectorXd& z, const std::pair<int,int>& edge) override;
};

#endif
