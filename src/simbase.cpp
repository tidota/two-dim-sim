// simbase.cpp

#include<cstdlib>

#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simbase.hpp>

// =============================================================================
SimBase::SimBase(const YAML::Node& doc):
  n_robots(doc["robots"].size()),
  errors(n_robots, 0),
  n_dim(doc["robots"][0].size()),
  max_time(doc["max_time"].as<double>()),
  sim_freq(doc["sim_freq"].as<double>()),
  deltaT(1.0/doc["sim_freq"].as<double>()),
  use_random_seed(doc["use_random_seed"].as<bool>()),
  random_seed(doc["random_seed"].as<unsigned int>()),
  sigmaM(doc["sigmaM"].as<double>()),
  sigmaF(doc["sigmaF"].as<double>()),
  sigmaS_R(doc["sigmaS_R"].as<double>()),
  sigmaS_T(doc["sigmaS_T"].as<double>()),
  sigmaGlobalLocR(doc["sigmaGlobalLocR"].as<double>()),
  sigmaGlobalLocT(doc["sigmaGlobalLocT"].as<double>()),
  repul_thresh(doc["repul_thresh"].as<double>()),
  repul_rate(doc["repul_rate"].as<double>()),
  attr_thresh(doc["attr_thresh"].as<double>()),
  attr_rate(doc["attr_rate"].as<double>()),
  fric_rate(doc["fric_rate"].as<double>()),
  ctrl_rate(doc["ctrl_rate"].as<double>()),
  ctrl_max(doc["ctrl_max"].as<double>()),
  comm_radius(doc["comm_radius"].as<double>()),
  topology(doc["topology"].as<std::string>()),
  enable_random_order(doc["enable_random_order"].as<bool>()),
  enable_prob_update(doc["enable_prob_update"].as<bool>()),
  prob_update_p(doc["prob_update_p"].as<double>()),
  enable_update_step(doc["enable_update_step"].as<bool>()),
  enable_global_loc(doc["enable_global_loc"].as<bool>()),
  enable_bidirectional(doc["enable_bidirectional"].as<bool>()),
  global_loc_steps(doc["global_loc_steps"].as<int>()),
  plot_interval(doc["plot_interval"].as<int>())
{
  for(int i = 0; i < n_robots; ++i)
  {
    VectorXd buff(n_dim);
    for (int j = 0; j < n_dim; ++j)
    {
      buff(j) = doc["robots"][i][j].as<double>();
    }
    robots.push_back(buff);
    vels.push_back(VectorXd::Zero(n_dim));
    accs.push_back(VectorXd::Zero(n_dim));
  }

  if (topology == "star")
  {
    for (int i = 1; i < n_robots; ++i)
    {
      edges.push_back(std::pair<int, int>(0, i));
    }
  }
  else if (topology == "loop")
  {
    for (int i = 0; i < n_robots; ++i)
    {
      edges.push_back(std::pair<int, int>(i, (i + 1)%n_robots));
    }
  }
  else if (topology == "complete")
  {
    for (int i = 0; i < n_robots - 1; ++i)
    {
      for (int j = i + 1; j < n_robots; ++j)
      {
        edges.push_back(std::pair<int, int>(i, j));
      }
    }
  }
  else if (topology == "dynamic")
  {
    for (int i = 0; i < n_robots - 1; ++i)
    {
      for (int j = i + 1; j < n_robots; ++j)
      {
        VectorXd diff = robots[j] - robots[i];
        if (diff.norm() <= comm_radius)
        {
          edges.push_back(std::pair<int, int>(i, j));
        }
      }
    }
  }
  else // unknown
  {
    std::cout << "Error: unknown topology type [" << topology << "]" << std::endl;
    exit(-1);
  }

  unsigned int seed;
  if (use_random_seed)
  {
    seed = random_seed;
  }
  else
  {
    std::random_device rd{};
    seed = rd();
  }
  gen.seed(seed);

  t_step = 0;
  t = 0;
}

// =============================================================================
void SimBase::printSimInfo()
{
  std::cout << "network topology: " << topology << std::endl;
  std::cout << "max time: " << max_time << std::endl;
  std::cout << "simulation update frequency (Hz): " << sim_freq << std::endl;
  std::cout << "sigmaM: " << sigmaM << std::endl;
  std::cout << "sigmaF: " << sigmaF << std::endl;
  std::cout << "sigmaS_R: " << sigmaS_R << std::endl;
  std::cout << "sigmaS_T: " << sigmaS_T << std::endl;
  std::cout << "sigmaGlobalLocR: " << sigmaGlobalLocR << std::endl;
  std::cout << "sigmaGlobalLocT: " << sigmaGlobalLocT << std::endl;
  std::cout << "ctrl_rate: " << ctrl_rate << std::endl;
  std::cout << "ctrl_max: " << ctrl_max << std::endl;
  std::cout << "# of robots: " << n_robots << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "robot[" << i << "]: ";
    for (int idim = 0; idim < 2; ++idim)
      std::cout << robots[i](idim) << ((idim == (2 - 1))? "": ", ");
    std::cout << std::endl;
  }
  std::cout << "deltaT: " << deltaT << std::endl;
  std::cout << "network topology (" << topology << "):";
  for (auto edge: edges)
  {
    std::cout << "(" << edge.first << ", " << edge.second << "), ";
  }
  std::cout << std::endl;
}

// =============================================================================
void SimBase::plot()
{
  if (t_step % plot_interval == 0)
    plotImpl();
}

// =============================================================================
int SimBase::stepForward()
{
  t = ++t_step * deltaT;
  return t_step;
}

// =============================================================================
bool SimBase::isDone()
{
  return t_step * deltaT > max_time;
}

// =============================================================================
void SimBase::updateSim()
{
  // === update simulation env. ===
  // NOTE: robots are forming a circle.

  // accelerations for robots:
  for (int i = 0; i < n_robots; ++i)
  {
    accs[i] = VectorXd::Zero(n_dim);
  }
  // Repulsive force
  for (int i = 0; i < n_robots - 1; ++i)
  {
    for (int j = i + 1; j < n_robots; ++j)
    {
      // if it is too close, apply a repulsive acceleration.
      VectorXd diff = robots[i] - robots[j];
      double dist = repul_thresh - diff.norm();
      if (0 < dist && dist < repul_thresh)
      {
        diff = diff / diff.norm();
        accs[i] += dist * dist * diff * repul_rate;
        accs[j] -= dist * dist * diff * repul_rate;
      }
    }
  }
  // Attractive force
  for (int i = 0; i < n_robots; ++i)
  {
    int j = (i + 1) % n_robots;

    // if it is too far, apply an attractive acceleration.
    VectorXd diff = robots[i] - robots[j];
    double dist = diff.norm() - attr_thresh;
    if (dist > 0)
    {
      diff = diff / diff.norm();
      accs[i] -= dist * dist * diff * attr_rate;
      accs[j] += dist * dist * diff * attr_rate;
    }
  }
  // Friction
  for (int i = 0; i < n_robots; ++i)
  {
    accs[i] -= fric_rate * vels[i];
  }

  // for all robots, apply the resulted accelerations
  for (int i = 0; i < n_robots; ++i)
  {
    std::normal_distribution<> motion_noise(0, sigmaM*std::fabs(vels[i].norm()));

    // update the position by the current velocity.
    robots[i] += deltaT * (vels[i] + motion_noise(gen) * VectorXd::Ones(n_dim));

    // update the velocity by the resulted acceleration.
    vels[i] += deltaT * accs[i];
  }

  // update connectivity based on the topology mode
  if (topology == "dynamic")
  {
    edges.clear();
    for (int i = 0; i < n_robots - 1; ++i)
    {
      for (int j = i + 1; j < n_robots; ++j)
      {
        VectorXd diff = robots[j] - robots[i];
        if (diff.norm() <= comm_radius)
        {
          edges.push_back(std::pair<int, int>(i, j));
        }
      }
    }
  }
}

// =============================================================================
void SimBase::globalLoc()
{
  if (t_step % global_loc_steps == 0)
  {
    // take measurement (communication)
    std::normal_distribution<> global_locR_noise(0, sigmaGlobalLocR);
    std::normal_distribution<> global_locT_noise(0, sigmaGlobalLocT);
    VectorXd z(2);
    z(0) = robots[0].norm() + global_locR_noise(gen);
    z(1) = std::atan2(robots[0](1), robots[0](0)) + global_locT_noise(gen);

    // it is supposed to be followed by the inherited function.
    if (enable_global_loc)
      globalLocImpl(z);
  }
}

// =============================================================================
void SimBase::mutualLoc()
{
  // === estimation update ===
  // for all edges of network
  std::vector<int> indx_list(edges.size());
  for (int i = 0; i < static_cast<int>(edges.size()); ++i)
    indx_list[i] = i;
  if (enable_random_order) // to add a flag to shuffle them
    std::random_shuffle(indx_list.begin(), indx_list.end());
  // mutual measurement
  for (auto indx: indx_list)
  {
    // throw a die and decide to update it or not.
    std::uniform_real_distribution<> dis(0, 1.0);
    double val = dis(gen);
    if (enable_prob_update && val > prob_update_p)
    {
      continue;
    }

    // get indexes
    auto edge = edges[indx];

    // difference between the two robots
    VectorXd diff = robots[edge.second] - robots[edge.first];
    // take measurement (communication)
    VectorXd z(2);
    std::normal_distribution<> sensor_noiseR(0, sigmaS_R);
    std::normal_distribution<> sensor_noiseT(0, sigmaS_T);
    z(0) = diff.norm() + sensor_noiseR(gen);
    z(1) = std::atan2(diff(1), diff(0)) + sensor_noiseT(gen);

    // if the update step is not enabled, just skip the rest. This just
    // generates random numbers which are the same sequence for the enabled
    // case so that the grand-truth trajectories are the same in both cases.
    if (!enable_update_step)
    {
      continue;
    }

    // it is supposed to be followed by the inherited function.
    mutualLocImpl(z, edge);
  }
}