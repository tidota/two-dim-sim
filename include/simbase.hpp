// simbase.hpp
//
// This is a base class providing a fundamental framework.

#ifndef SIMBASE_HPP
#define SIMBASE_HPP

#include <random>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class SimBase
{
  // robot's locations
  protected: std::vector<VectorXd> robots;

  // number of robots
  protected: const int n_robots;
  // velocities
  protected: std::vector<VectorXd> vels;
  // accelerations
  protected: std::vector<VectorXd> accs;
  // errors
  protected: std::vector<double> errors;
  // mode of estimation
  protected: const int mode;
  // use orientation
  protected: const bool use_orientation;
  // number of dimensions
  protected: const int n_dim;

  // max time
  protected: const double max_time;
  // simulation update frequendy
  protected: const double sim_freq;
  // phase, i.e., deta T
  protected: const double deltaT; // = 1.0/sim_freq;
  // time step
  protected: int t_step;
  // time
  protected: double t; // t_step * deltaT
  public: double getProgress()
  {
    return (double)t/max_time;
  }

  // using random seed?
  protected: const bool use_random_seed;
  protected: const unsigned int random_seed;

  // weight of motion noise
  protected: const double sigmaM;
  // sigma for floating effect
  protected: const double sigmaF;
  // sigma for sensor model
  protected: const double sigmaS_R;
  protected: const double sigmaS_T;
  // sigma for global localizaiton
  protected: const double sigmaGlobalLocR;
  protected: const double sigmaGlobalLocT;

  // thresholds and rates to generate accelerations
  protected: double repul_thresh;
  protected: const double repul_rate;
  protected: double attr_thresh;
  protected: const double attr_rate;
  protected: const double fric_rate;
  public: virtual void changeControl(
    const double repul_th, const double attr_th) final
  {
    repul_thresh = repul_th;
    attr_thresh = attr_th;
  }
  public: virtual void scaleControl(
    const double repul_rate, const double attr_rate) final
  {
    repul_thresh *= repul_rate;
    attr_thresh *= attr_rate;
  }

  // control rate
  protected: const double ctrl_rate;
  // control max
  protected: const double ctrl_max;

  // communication radius
  protected: const double comm_radius;
  // topology/connectivity
  protected: const std::string topology;
  protected: std::vector< std::pair<int,int> > edges;

  // random order to update estimates?
  protected: const bool enable_random_order;
  // probabilistically update?
  protected: const bool enable_prob_update;
  protected: const double prob_update_p;

  // enable_update_step
  protected: const bool enable_update_step;
  // global localization
  protected: const bool enable_global_loc;
  // enable_bidirectional
  protected: const bool enable_bidirectional;
  protected: std::mt19937 gen;
  protected: const int global_loc_steps;
  protected: const int plot_interval;

  // buffer for showing covariances
  protected: std::vector< std::vector< std::vector<double> > > cov_buff;
  protected: bool show_covs;

  public: SimBase(const YAML::Node& doc);
  public: virtual ~SimBase(){}

  // print sim infor
  public: virtual void printSimInfo();
  // start logging
  public: virtual bool startLog(const std::string& fname) = 0;
  // end logging
  public: virtual void endLog() = 0;
  // give one line for the time step to plot
  public: virtual void plot() final;
  protected: virtual void plotImpl() = 0;

  // increment the time step
  public: virtual int stepForward() final;
  // check if the simulation is done
  public: virtual bool isDone() final;
  // update the simulation environment
  public: virtual void updateSim() final;

  public: virtual void predict() = 0;
  public: virtual void globalLoc() final;
  protected: virtual void globalLocImpl(const VectorXd& z) = 0;
  public: virtual void mutualLoc() final;
  protected: virtual void mutualLocImpl(
    const VectorXd& z, const std::pair<int,int>& edge) = 0;
  public: virtual void calcErrors() = 0;
};

#endif
