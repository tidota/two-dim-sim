// simparam.hpp

// this is a middle layer class providing abstract model of parametric filter.

#ifndef SIMPARAM_HPP
#define SIMPARAM_HPP

#include <fstream>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simbase.hpp>

using namespace std;
using namespace Eigen;

class SimParam: public SimBase
{
  // weight of motion model
  private: const double alpha1M;
  private: const double alpha2M;
  // constant value of motion model
  private: const double betaM;

  // means on robot's locations
  protected: std::vector<VectorXd> means;
  // variance
  protected: std::vector<MatrixXd> vars;

  // buffers for covariances to plot
  private: std::vector<VectorXd> last_loc;
  private: std::vector<VectorXd> last_mean;
  private: std::vector<MatrixXd> last_var;

  // file for plotting
  private: std::ofstream fout;

  public: SimParam(const YAML::Node& doc);
  public: virtual ~SimParam(){}

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
  protected: virtual void globalLocModeImpl(
    const MatrixXd& H, const MatrixXd& Q, const VectorXd& z_diff);
  protected: virtual void mutualLocImpl(
    const VectorXd& z, const std::pair<int,int>& edge) override final;
  protected: virtual void mutualLocModeImpl(
    const VectorXd& z, const std::pair<int,int>& edge,
    const MatrixXd& H1, const MatrixXd& H2, const MatrixXd& Q,
    const VectorXd& z_diff) = 0;
  public: void calcErrors() override final;
};

#endif
