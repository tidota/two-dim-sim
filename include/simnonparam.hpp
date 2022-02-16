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
  protected: int n_particles;
  protected: std::vector< std::vector<VectorXd> > ests;
  protected: std::vector< std::vector<double> > est_oris;

  // buffers for covariances to plot
  private: std::vector<VectorXd> last_loc;
  private: std::vector< std::vector<VectorXd> > last_est;
  private: std::vector< std::vector<double> > last_est_ori;

  // buffer for calculation of average errors
  private: std::vector<double> cumul_errors;
  // buffer for calculation of average errors for each phase
  // (a new phase is added when scaleControl is called)
  private: std::vector< std::vector<double> > cumul_errors_phase;
  // time steps for each phase
  // (a new phase is added when scaleControl is called)
  private: std::vector<int> time_steps_phase;

  // buffer for calculation of average errors in orientation
  // It is calculated if use_orientation is set to true.
  private: std::vector<double> cumul_errors_ori;
  // buffer for calculation of average errors in orientation for each phase
  // (a new phase is added when scaleControl is called)
  private: std::vector< std::vector<double> > cumul_errors_ori_phase;

  // random value generator for sampling
  protected: std::mt19937 gen_pf;
  // using random seed?
  protected: const bool use_random_seed_pf;
  protected: const unsigned int random_seed_pf;
  // constant param for the evaluation process
  private: const double gl_eval_cons;
  // plot particles?
  private: const bool plot_particles;

  // file for plotting (trajectories)
  private: std::ofstream fout;
  // file for plotting (particles at the end)
  private: std::ofstream fout_pf;

  public: SimNonParam(const YAML::Node& doc);
  public: virtual ~SimNonParam(){}

  // add items for the new phase
  // (this function is called by scaleControl in SimBase)
  protected: virtual void startNewPhase() override
  {
    cumul_errors_phase.push_back(std::vector<double>(n_robots, 0));
    cumul_errors_ori_phase.push_back(std::vector<double>(n_robots, 0));
    time_steps_phase.push_back(0);
  }

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

  protected: int drawRandIndx(const std::vector<double>& cumul_weights);
};

#endif
