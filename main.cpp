/*

The main part of the simulation.

*/

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main (int argc, char** argv)
{
  std::ifstream fin("settings.yaml");
  YAML::Node doc = YAML::Load(fin);

  // get max time
  const double max_time = doc["max_time"].as<double>();
  // get simulation update frequendy
  const double sim_freq = doc["sim_freq"].as<double>();
  // phase, i.e., deta T
  const double deltaT = 1.0/sim_freq;

  // weight of motion noise
  const double sigmaM = doc["sigmaM"].as<double>();
  // get sigma for floating effect
  const double sigmaF = doc["sigmaF"].as<double>();
  // get sigma for sensor model
  const double sigmaS = doc["sigmaS"].as<double>();
  // get sigma for global localizaiton
  const double sigmaGlobalLocR = doc["sigmaGlobalLocR"].as<double>();
  const double sigmaGlobalLocT = doc["sigmaGlobalLocT"].as<double>();

  // control rate
  const double ctrl_rate = doc["ctrl_rate"].as<double>();
  // control max
  const double ctrl_max = doc["ctrl_max"].as<double>();

  // weight of motion model
  const double alphaM = doc["alphaM"].as<double>();
  // constant value of motion model
  const double betaM = doc["betaM"].as<double>();

  // mode
  const int mode = doc["mode"].as<int>();

  // enable_update_step
  const bool enable_update_step = doc["enable_update_step"].as<bool>();

  // global localization
  const bool enable_global_loc = doc["enable_global_loc"].as<bool>();
  // global localization at every specified steps.
  const int global_loc_steps = doc["global_loc_steps"].as<int>();

  // plotting interval
  const int plot_interval = doc["plot_interval"].as<int>();

  // robot's locations
  std::vector<VectorXd> robots;
  // means on robot's locations
  std::vector<VectorXd> means;
  // variance
  std::vector<MatrixXd> vars;
  // number of robots
  int n_robots = doc["robots"].size();
  // velocities
  std::vector<VectorXd> vels;
  // accelerations
  std::vector<VectorXd> accs;
  // errors
  std::vector<double> errors(n_robots, 0);
  // number of dimensions
  int n_dim = doc["robots"][0].size();
  // for all robots
  for(int i = 0; i < n_robots; ++i)
  {
    // init location x
    VectorXd buff(n_dim);
    for (int j = 0; j < n_dim; ++j)
    {
      buff(j) = doc["robots"][i][j].as<double>();
    }
    robots.push_back(buff);
    means.push_back(buff);
    vars.push_back(MatrixXd::Zero(n_dim, n_dim));
    vels.push_back(VectorXd::Zero(n_dim));
    accs.push_back(VectorXd::Zero(n_dim));
  }

  // topology/connectivity
  std::string topology = doc["topology"].as<std::string>();
  std::vector< std::pair<int,int> > edges;

  if (topology == "star")
  {
    for (int i = 1; i < n_robots; ++i)
    {
      edges.push_back(std::pair<int, int>(0, i));
    }
  }
  else if (topology == "loop")
  {
    for (int i = 0; i < n_robots - 1; ++i)
    {
      edges.push_back(std::pair<int, int>(i, i + 1));
    }
  }
  else if (topology == "complete")
  {
    for (int i = 1; i < n_robots - 1; ++i)
    {
      for (int j = i + 1; j < n_robots; ++j)
      {
        edges.push_back(std::pair<int, int>(i, j));
      }
    }
  }
  else if (topology == "dynamic")
  {
    // TODO
  }
  else // unknown
  {
    std::cout << "Error: unknown topology type [" << topology << "]" << std::endl;
    return -1;
  }

  // random order to update estimates?
  bool enable_random_order = doc["enable_random_order"].as<bool>();
  // probabilistically update?
  bool enable_prob_update = doc["enable_prob_update"].as<bool>();
  double prob_update_p = doc["prob_update_p"].as<double>();
  // allow loopy updates (the last robot can communicate with the first one)
  bool enable_loop = doc["enable_loop"].as<bool>();
  fin.close();

  std::cout << "network topology: " << topology << std::endl;
  std::cout << "max time: " << max_time << std::endl;
  std::cout << "simulation update frequency (Hz): " << sim_freq << std::endl;
  std::cout << "sigmaM: " << sigmaM << std::endl;
  std::cout << "sigmaF: " << sigmaF << std::endl;
  std::cout << "sigmaS: " << sigmaS << std::endl;
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
    std::cout << "| ";
    std::cout << "means[" << i << "]: ";
    for (int idim = 0; idim < 2; ++idim)
      std::cout << means[i](idim) << ((idim == (2 - 1))? "": ", ");
    std::cout << "| ";
    std::cout << "vars[" << i << "]: " << vars[i].determinant() << std::endl;
  }
  std::cout << "deltaT: " << deltaT << std::endl;
  std::cout << "network topology (" << topology << "):";
  for (auto edge: edges)
  {
    std::cout << "(" << edge.first << ", " << edge.second << "), ";
  }
  std::cout << std::endl;

  // random numbers
  std::random_device rd{};
  std::mt19937 gen{rd()};

  // print header
  std::cout << "    t   |";
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "R[" << setw(7 * n_dim - 7) << i << "].x |";
    std::cout << "R[" << setw(7 * n_dim - 7) << i << "].m |";
    std::cout << " det(v) |";
  }
  std::cout << std::endl;

  // prep output file
  ofstream fout ("output.dat");

  // main part of simulation
  for (int t_step = 0; t_step * deltaT <= max_time; ++t_step)
  {
    double t = t_step * deltaT;
    // print the current states
    if ((int)(t * sim_freq) % plot_interval == 0)
    {
      std::cout << std::fixed << std::setprecision(2);
      std::cout << std::right << std::setw(8) << t << "|";
      fout << std::fixed << std::setprecision(3);
      fout << std::right << std::setw(8) << t << " ";

      for (int i = 0; i < n_robots; ++i)
      {
        for (int j = 0; j < n_dim; ++j)
        {
          std::cout << std::right << std::setw(6) << robots[i](j);
          if (j < n_dim - 1)
            std::cout << ",";
          else
            std::cout << "|";
        }
        for (int j = 0; j < n_dim; ++j)
        {
          std::cout << std::right << std::setw(6) << means[i](j);
          if (j < n_dim - 1)
            std::cout << ",";
          else
            std::cout << "|";
        }
        std::cout << std::right << std::setw(8) << std::sqrt(vars[i].determinant()) << "|";

        for (int j = 0; j < n_dim; ++j)
        {
          fout << std::right << std::setw(9) << robots[i](j);
          fout << " ";
        }
        for (int j = 0; j < n_dim; ++j)
        {
          fout << std::right << std::setw(9) << means[i](j);
          fout << " ";
        }
        fout << std::right << std::setw(9) << std::sqrt(vars[i].determinant()) << " ";
      }

      std::cout << std::endl;
      fout << std::endl;
    }

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
        const double thresh = 10.0;
        const double repul_rate = 0.1;
        VectorXd diff = robots[i] - robots[j];
        double dist = thresh - diff.norm();
        if (0 < dist && dist < thresh)
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
      const double thresh = 3.0;
      const double attr_rate = 0.1;
      VectorXd diff = robots[i] - robots[j];
      double dist = diff.norm() - thresh;
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
      const double fric_rate = 0.6;
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
      // TODO
    }

    // === prediction ===
    // for all robots
    for (int i = 0; i < n_robots; ++i)
    {
      // motion model
      means[i] = means[i] + MatrixXd::Identity(n_dim, n_dim) * vels[i] * deltaT;
      MatrixXd V = MatrixXd::Identity(n_dim, n_dim) * deltaT;
      MatrixXd M = (alphaM * alphaM) * MatrixXd::Identity(n_dim, n_dim).cwiseProduct(vels[i] * vels[i].transpose());
      vars[i] = vars[i] + V.transpose() * M * V;
    }

    // === estimation update ===
    std::vector<VectorXd> means_buff(means);
    std::vector<MatrixXd> vars_buff(vars);

    // global localization.
    if (enable_global_loc && t_step % global_loc_steps == 0)
    {
      // take measurement (communication)
      VectorXd z(2);
      std::normal_distribution<> global_locR_noise(0, sigmaGlobalLocR);
      std::normal_distribution<> global_locT_noise(0, sigmaGlobalLocT);
      z(0) = robots[0].norm() + global_locR_noise(gen);
      z(1) = std::atan2(robots[0](1), robots[0](0)) + global_locT_noise(gen);

      // update its location estimate
      VectorXd z_hat(2);// = means_buff[0];
      z_hat(0) = means[0].norm();
      z_hat(1) = std::atan2(means[0](1), means[0](0));
      double q = means[0].squaredNorm();
      MatrixXd H(2,n_dim);
      H(0, 0) = means[0](0)/std::sqrt(q);
      H(1, 0) = -means[0](1)/q;
      H(0, 1) = means[0](1)/std::sqrt(q);
      H(1, 1) = means[0](0)/q;
      MatrixXd Q = MatrixXd::Zero(2,2);
      Q(0,0) = sigmaGlobalLocR * sigmaGlobalLocR;
      Q(1,1) = sigmaGlobalLocT * sigmaGlobalLocT;
      MatrixXd St = H * vars_buff[0] * H.transpose() + Q;
      MatrixXd K = vars_buff[0] * H.transpose() * St.inverse();
      means_buff[0] += K * (z - z_hat);
      vars_buff[0] = (MatrixXd::Identity(n_dim, n_dim) - K * H) * vars_buff[0];
    }

    // for all edges of network
    std::vector<int> indx_list(edges.size());
    if (true) // to add a flag to shuffle them
      std::random_shuffle(indx_list.begin(), indx_list.end());
    for (auto indx: indx_list)
    {
      // mutual measurement
      auto edge = edges[indx];

      // take measurement (communication)
      VectorXd z(2);
      std::normal_distribution<> global_locR_noise(0, sigmaGlobalLocR);
      std::normal_distribution<> global_locT_noise(0, sigmaGlobalLocT);
      z(0) = robots[0].norm() + global_locR_noise(gen);
      z(1) = std::atan2(robots[0](1), robots[0](0)) + global_locT_noise(gen);

      // update its location estimate
      VectorXd z_hat(2);// = means_buff[0];
      z_hat(0) = means[0].norm();
      z_hat(1) = std::atan2(means[0](1), means[0](0));
      double q = means[0].squaredNorm();
      MatrixXd H(2,n_dim);
      H(0, 0) = means[0](0)/std::sqrt(q);
      H(1, 0) = -means[0](1)/q;
      H(0, 1) = means[0](1)/std::sqrt(q);
      H(1, 1) = means[0](0)/q;
      MatrixXd Q = MatrixXd::Zero(2,2);
      Q(0,0) = sigmaGlobalLocR * sigmaGlobalLocR;
      Q(1,1) = sigmaGlobalLocT * sigmaGlobalLocT;
      MatrixXd St = H * vars_buff[0] * H.transpose() + Q;
      MatrixXd K = vars_buff[0] * H.transpose() * St.inverse();
      means_buff[0] += K * (z - z_hat);
      vars_buff[0] = (MatrixXd::Identity(n_dim, n_dim) - K * H) * vars_buff[0];
    }

    // apply the updated estimations
    if (enable_update_step)
    {
      for (int i = 0; i < n_robots; ++i)
      {
        means[i] = means_buff[i];
        vars[i] = vars_buff[i];
      }
    }

    // calculate errors
    for (int i = 0; i < n_robots; ++i)
    {
      // TODO
      //errors[i] += std::fabs(means[i] - robots[i]);
    }
  }

  fout.close();

  // output of gnuplot command
  // std::cout << std::endl;
  // std::cout << "gnuplot command" << std::endl;
  // for (int i = 0; i < n_robots; ++i)
  // {
  //   if (i == 0)
  //     std::cout << "plot ";
  //   else
  //     std::cout << "     ";
  //   std::cout << "\"output.dat\" u 1:" << std::to_string(2+i*3)
  //             << " title \"x" << std::to_string(1+i) << "\", \\"
  //             << std::endl;
  // }
  // for (int i = 0; i < n_robots; ++i)
  // {
  //   std::cout << "     ";
  //   std::cout << "\"output.dat\" u 1:"
  //             << std::to_string(3+i*3) << ":" << std::to_string(4+i*3)
  //             << " with errorbars title \"m" << std::to_string(1+i) << "\"";
  //   if (i < n_robots - 1)
  //     std::cout << ", \\";
  //   std::cout << std::endl;
  // }
  // std::cout << std::endl;

  // display errors
  double total_error = 0;
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "robot[" << i << "]'s average error:" << (errors[i]/(max_time*sim_freq)) << std::endl;
    total_error += errors[i];
  }
  std::cout << "overall average error: " << (total_error/(max_time*sim_freq)/n_robots) << std::endl;


  // just for test of eigen
  Matrix3d A = Matrix3d::Identity();
  std::cout << "A = I: " << std::endl << A << std::endl;

  A(1,1) = 5;
  std::cout << "A(1,1) = 5: " << std::endl << A << std::endl;

  return 0;
}
