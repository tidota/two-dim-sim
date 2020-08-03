// simnonparam.hpp

// this is a middle layer class providing abstract model of nonparametric
// filter.

#ifndef SIMNONPARAM_HPP
#define SIMNONPARAM_HPP

#include <fstream>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simbase.hpp>

using namespace std;
using namespace Eigen;

class SimNonParam: public SimBase
{
  // weight of motion model
  private: const double alpha1M;
  private: const double alpha2M;
  // constant value of motion model
  private: const double betaM;

  // means on robot's locations
  protected: std::vector< std::vector<VectorXd> > ests;
  protected: int n_particles;

  // buffers for covariances to plot
  private: std::vector<VectorXd> last_loc;
  private: std::vector< std::vector<VectorXd> > last_est;

  // file for plotting (trajectories)
  private: std::ofstream fout;
  // file for plotting (particles at the end)
  private: std::ofstream fout_pf;

  public: SimNonParam(const YAML::Node& doc);
  public: virtual ~SimNonParam(){}

  // print sim infor
  public: virtual void printSimInfo() override;
  // start logging
  public: virtual bool startLog(const std::string& fname) override final;
  // end logging
  public: virtual void endLog() override final;
  // give one line for the time step to plot
  public: virtual void plotImpl() override final;

  public: virtual void predict() override;
  protected: virtual void globalLocImpl(const VectorXd& z) override final;
  protected: virtual void mutualLocImpl(
    const VectorXd& z, const std::pair<int,int>& edge) = 0;
  public: void calcErrors() override final;
};

#endif
