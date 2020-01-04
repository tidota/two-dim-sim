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
  const double sigmaS_R = doc["sigmaS_R"].as<double>();
  const double sigmaS_T = doc["sigmaS_T"].as<double>();
  // get sigma for global localizaiton
  const double sigmaGlobalLocR = doc["sigmaGlobalLocR"].as<double>();
  const double sigmaGlobalLocT = doc["sigmaGlobalLocT"].as<double>();

  // thresholds and rates to generate accelerations
  const double repul_thresh = doc["repul_thresh"].as<double>();
  const double repul_rate = doc["repul_rate"].as<double>();
  const double attr_thresh = doc["attr_thresh"].as<double>();
  const double attr_rate = doc["attr_rate"].as<double>();
  const double fric_rate = doc["fric_rate"].as<double>();

  // control rate
  const double ctrl_rate = doc["ctrl_rate"].as<double>();
  // control max
  const double ctrl_max = doc["ctrl_max"].as<double>();

  // weight of motion model
  const double alpha1M = doc["alpha1M"].as<double>();
  const double alpha2M = doc["alpha2M"].as<double>();
  // constant value of motion model
  const double betaM = doc["betaM"].as<double>();

  // mode
  const int mode = doc["mode"].as<int>();
  // params for mode2
  const double mode2_rate1 = doc["mode2_rate1"].as<double>();
  const double mode2_rate2 = doc["mode2_rate2"].as<double>();
  // params for mode3
  const double mode3_omega = doc["mode3_omega"].as<double>();

  // communication radius
  const double comm_radius = doc["comm_radius"].as<double>();

  // enable_update_step
  const bool enable_update_step = doc["enable_update_step"].as<bool>();

  // enable_bidirectional
  const bool enable_bidirectional = doc["enable_bidirectional"].as<bool>();

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
    return -1;
  }

  // random order to update estimates?
  bool enable_random_order = doc["enable_random_order"].as<bool>();
  // probabilistically update?
  bool enable_prob_update = doc["enable_prob_update"].as<bool>();
  double prob_update_p = doc["prob_update_p"].as<double>();
  fin.close();

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
    std::cout << " err |";
  }
  std::cout << std::endl;

  // prep output file
  ofstream fout ("output.dat");

  // buffers for covariances to plot
  std::vector<VectorXd> last_loc(n_robots);
  std::vector<VectorXd> last_mean(n_robots);
  std::vector<MatrixXd> last_var(n_robots);

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
        std::cout << std::right << std::setw(5) << (robots[i] - means[i]).norm() << "|";

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
        for (int j = 0; j < n_dim; ++j)
        {
          for (int k = 0; k < n_dim; ++k)
          {
            fout << std::right << std::setw(9) << vars[i](j,k);
            fout << " ";
          }
        }
        fout << std::right << std::setw(8)
             << (robots[i] - means[i]).norm() << " ";

        // store data to the buffer to plot
        last_loc[i] = robots[i];
        last_mean[i] = means[i];
        last_var[i] = vars[i];
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

    // === prediction ===
    // for all robots
    for (int i = 0; i < n_robots; ++i)
    {
      // motion model
      means[i] = means[i] + MatrixXd::Identity(n_dim, n_dim) * vels[i] * deltaT;

      MatrixXd FloatEffect = MatrixXd::Identity(2, 2) * betaM;
      double v = vels[i].norm();
      MatrixXd EigenVecs(2, 2);
      EigenVecs(0, 0) = vels[i](0) / v;
      EigenVecs(1, 0) = vels[i](1) / v;
      EigenVecs(0, 1) = -vels[i](1) / v;
      EigenVecs(1, 1) = vels[i](0) / v;
      MatrixXd EigenVals = MatrixXd::Identity(2, 2);
      EigenVals(0, 0) = alpha1M * v * v;
      EigenVals(1, 1) = alpha2M * v * v;
      // transpose == inverse because eigen vectors are perpendicular to each other
      MatrixXd M = EigenVecs * EigenVals * EigenVecs.transpose() + FloatEffect;
      MatrixXd V = MatrixXd::Identity(n_dim, n_dim) * deltaT;
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
      VectorXd z_diff = z - z_hat;
      if (z_diff(1) > M_PI)
      {
        z_diff(1) -= 2*M_PI;
      }
      else if (z_diff(1) <= -M_PI)
      {
        z_diff(1) += 2*M_PI;
      }
      means_buff[0] += K * z_diff;
      vars_buff[0] = (MatrixXd::Identity(n_dim, n_dim) - K * H) * vars_buff[0];
      // apply the updated estimations
      if (enable_update_step)
      {
        means[0] = means_buff[0];
        vars[0] = vars_buff[0];
      }
    }

    // for all edges of network
    std::vector<int> indx_list(edges.size());
    for (int i = 0; i < static_cast<int>(edges.size()); ++i)
      indx_list[i] = i;
    if (enable_random_order) // to add a flag to shuffle them
      std::random_shuffle(indx_list.begin(), indx_list.end());
    // mutual measurement
    for (auto indx: indx_list)
    {
      // if probabilistic update is enabled
      if (enable_prob_update)
      {
        // throw a die and decide to update it or not.
        std::uniform_real_distribution<> dis(0, 1.0);
        double val = dis(gen);
        if (val > prob_update_p)
          continue;
      }

      // get indexes
      auto edge = edges[indx];
      VectorXd diff = robots[edge.second] - robots[edge.first];
      VectorXd diff_hat = means[edge.second] - means[edge.first];

      // take measurement (communication)
      VectorXd z(2);
      std::normal_distribution<> sensor_noiseR(0, sigmaS_R);
      std::normal_distribution<> sensor_noiseT(0, sigmaS_T);
      z(0) = diff.norm() + sensor_noiseR(gen);
      z(1) = std::atan2(diff(1), diff(0)) + sensor_noiseT(gen);

      // update its location estimate
      VectorXd z_hat(2);// = means_buff[0];
      z_hat(0) = diff_hat.norm();
      z_hat(1) = std::atan2(diff_hat(1), diff_hat(0));
      double q = diff_hat.squaredNorm();

      MatrixXd Q = MatrixXd::Zero(2,2);
      Q(0,0) = sigmaS_R * sigmaS_R;
      Q(1,1) = sigmaS_T * sigmaS_T;

      MatrixXd H1(2,n_dim);
      H1(0, 0) = -diff_hat(0)/std::sqrt(q);
      H1(1, 0) = diff_hat(1)/q;
      H1(0, 1) = -diff_hat(1)/std::sqrt(q);
      H1(1, 1) = -diff_hat(0)/q;

      MatrixXd H2(2,n_dim);
      H2(0, 0) = diff_hat(0)/std::sqrt(q);
      H2(1, 0) = -diff_hat(1)/q;
      H2(0, 1) = diff_hat(1)/std::sqrt(q);
      H2(1, 1) = diff_hat(0)/q;

      MatrixXd St1 = H1 * vars_buff[edge.first] * H1.transpose() + Q;
      MatrixXd St2 = H2 * vars_buff[edge.second] * H2.transpose() + Q;

      if (mode == 0)
      {
        St1 = H1 * vars_buff[edge.first] * H1.transpose() + Q;
        St2 = H2 * vars_buff[edge.second] * H2.transpose() + Q;
      }
      else if (mode == 1)
      {
        St1 = H1 * vars_buff[edge.first] * H1.transpose()
            + H2 * vars_buff[edge.second] * H2.transpose() + Q;
        St2 = St1;
      }
      else if (mode == 2)
      {
        vars_buff[edge.first] *= mode2_rate1;
        vars_buff[edge.second] *= mode2_rate2;
        St1 = H1 * vars_buff[edge.first] * H1.transpose()
            + H2 * vars_buff[edge.second] * H2.transpose() + Q;
        St2 = St1;
      }
      else if (mode == 3)
      {
        St1 = H1 * (vars_buff[edge.first] / mode3_omega) * H1.transpose()
            + H2 * (vars_buff[edge.second] / (1 - mode3_omega)) * H2.transpose()
            + Q;
        St2 = H1 * (vars_buff[edge.first] / (1 - mode3_omega)) * H1.transpose()
            + H2 * (vars_buff[edge.second] / mode3_omega) * H2.transpose()
            + Q;
        vars_buff[edge.first] /= mode3_omega;
        vars_buff[edge.second] /= mode3_omega;
      }

      VectorXd z_diff = z - z_hat;
      if (z_diff(1) > M_PI)
      {
        z_diff(1) -= 2*M_PI;
      }
      else if (z_diff(1) <= -M_PI)
      {
        z_diff(1) += 2*M_PI;
      }
      if (enable_bidirectional)
      {
        MatrixXd K1 = vars_buff[edge.first] * H1.transpose() * St1.inverse();
        means_buff[edge.first] += K1 * z_diff;
        vars_buff[edge.first] = (MatrixXd::Identity(n_dim, n_dim) - K1 * H1) * vars_buff[edge.first];
      }

      MatrixXd K2 = vars_buff[edge.second] * H2.transpose() * St2.inverse();
      means_buff[edge.second] += K2 * z_diff;
      vars_buff[edge.second] = (MatrixXd::Identity(n_dim, n_dim) - K2 * H2) * vars_buff[edge.second];

      // apply the updated estimations
      if (enable_update_step)
      {
        means[edge.first] = means_buff[edge.first];
        vars[edge.first] = vars_buff[edge.first];
        means[edge.second] = means_buff[edge.second];
        vars[edge.second] = vars_buff[edge.second];
      }
    }

    // calculate errors
    for (int i = 0; i < n_robots; ++i)
    {
      errors[i] += (means[i] - robots[i]).norm();
    }
  }

  fout.close();

  // output of gnuplot command
  const int off_next_robot = n_dim + n_dim + n_dim*n_dim + 1;
  std::cout << std::endl;
  std::cout << "~~~ gnuplot command (errors)~~~" << std::endl;
  std::cout << "clear" << std::endl;
  std::cout << "unset object" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    if (i == 0)
      std::cout << "plot ";
    else
      std::cout << "     ";
    std::cout << "\"output.dat\" u 1:"
              << std::to_string(2+(i+1)*off_next_robot-1)
              << " title \"err" << std::to_string(1+i) << "\" with line";
    if (i < n_robots - 1)
      std::cout << ", \\";
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "~~~ gnuplot command (locations) ~~~" << std::endl;
  std::cout << "clear" << std::endl;
  std::cout << "unset object" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    Eigen::EigenSolver<MatrixXd> s(vars[i]);
    auto eigen_val = s.eigenvalues();
    auto eigen_vec = s.eigenvectors();
    double var_ang
      = std::atan2(eigen_vec.col(0)[1].real(), eigen_vec.col(0)[0].real())
        / M_PI*180.0;
    std::cout << "set object " << std::to_string(i + 1)
              << " ellipse center "
              << std::to_string(last_mean[i](0)) << ","
              << std::to_string(last_mean[i](1))
              << " size "
              << std::to_string(std::sqrt(eigen_val[0].real()) * 2) << ","
              << std::to_string(std::sqrt(eigen_val[1].real()) * 2)
              << " angle "
              << std::to_string(var_ang)
              << " front fillstyle empty border -1" << std::endl;
  }
  std::cout << "h1 = 227/360.0" << std::endl;
  std::cout << "h2 = 40/360.0" << std::endl;
  std::cout << "set palette model HSV functions (1-gray)*(h2-h1)+h1,1,0.68"
            << std::endl;
  std::cout << "set size ratio -1" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    if (i == 0)
      std::cout << "plot   ";
    else
      std::cout << "replot ";
    std::cout << "\"output.dat\" u ";
    for (int j = 0; j < n_dim; ++j)
    {
      std::cout << std::to_string(2+i*off_next_robot+j);
      if (j < n_dim - 1)
        std::cout << ":";
      else
        std::cout << ":1";
    }
    std::cout << " title \"R" << std::to_string(1+i) << "\""
              << " with linespoints lt -1 lw 1.0 ps 3.0"
              << " pt " << std::to_string(i + 1)
              << " lc palette" << std::endl;
  }
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "replot ";
    std::cout << "\"output.dat\" u ";
    for (int j = 0; j < n_dim; ++j)
    {
      std::cout << std::to_string(2+n_dim+i*off_next_robot+j);
      if (j < n_dim - 1)
        std::cout << ":";
      else
        std::cout << ":1";
    }
    std::cout << " title \"Est" << std::to_string(1+i) << "\""
              << " with linespoints lt 1 lw 3.0 ps 3.0"
              << " pt " << std::to_string(i + 1)
              << " lc palette" << std::endl;
  }
  std::cout << std::endl;

  // display errors
  double total_error = 0;
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "robot[" << i << "]'s average error:" << (errors[i]/(max_time*sim_freq)) << std::endl;
    total_error += errors[i];
  }
  std::cout << "overall average error: " << (total_error/(max_time*sim_freq)/n_robots) << std::endl;

  return 0;
}
