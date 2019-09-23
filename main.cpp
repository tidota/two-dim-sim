/*

The main part of the simulation.
The parameter list:

- # of robots
- max range of sensory data
- size of the world

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

  // field size
  double field_size = doc["field_size"].as<double>();
  // get max time steps
  int max_time_steps = doc["max_time_steps"].as<int>();
  // get delta T
  double deltaT = doc["deltaT"].as<double>();
  // weight of motion noise
  double sigmaM = doc["sigmaM"].as<double>();
  // get sigma for floating effect
  double sigmaF = doc["sigmaF"].as<double>();
  // get sigma for sensor model
  double sigmaS = doc["sigmaS"].as<double>();
  // get sigma for global localizaiton
  double sigmaGlobalLoc = doc["sigmaGlobalLoc"].as<double>();
  // get # of robots
  int n_robots = doc["n_robots"].as<int>();
  // control rate
  double ctrl_rate = doc["ctrl_rate"].as<double>();
  // control max
  double ctrl_max = doc["ctrl_max"].as<double>();
  // robot's locations
  std::vector<double> robots(n_robots);
  // means on robot's locations
  std::vector<double> means(n_robots);
  // variance
  std::vector<double> vars(n_robots);
  // errors
  std::vector<double> errors(n_robots, 0);
  // mode
  // 0: default (each robot assumes other robot's locations are accurate)
  // 1: other robot's uncertainty is considered
  // 2: reduced uncertainty
  int mode = doc["mode"].as<int>();
  // enable_update_step
  bool enable_update_step = doc["enable_update_step"].as<bool>();
  // global localization
  bool global_loc = doc["global_loc"].as<bool>();
  // global localization at every specified steps.
  int global_loc_steps = doc["global_loc_steps"].as<int>();
  // for all robots
  for(unsigned int i = 0; i < doc["robots"].size(); ++i)
  {
    // init location x
    robots[i] = doc["robots"][i]["x"].as<double>();
    means[i] = robots[i];
    vars[i] = 0;
  }
  // random order to update estimates?
  bool random_order = doc["random_order"].as<bool>();
  // probabilistically update?
  bool prob_update = doc["prob_update"].as<bool>();
  double prob_update_p = doc["prob_update_p"].as<double>();
  // allow loopy updates (the last robot can communicate with the first one)
  bool loopy_updates = doc["enable_loop"].as<bool>();
  fin.close();

  std::cout << "field size: " << field_size << std::endl;
  std::cout << "max time steps: " << max_time_steps << std::endl;
  std::cout << "delta T: " << deltaT << std::endl;
  std::cout << "sigmaF: " << sigmaF << std::endl;
  std::cout << "sigmaS: " << sigmaS << std::endl;
  std::cout << "n_robots: " << n_robots << std::endl;
  std::cout << "ctrl_rate: " << ctrl_rate << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "robot[" << i << "]: " << robots[i] << std::endl;
    std::cout << "means[" << i << "]: " << means[i] << std::endl;
    std::cout << "vars[" << i << "]: " << vars[i] << std::endl;
  }

  // random numbers
  std::random_device rd{};
  std::mt19937 gen{rd()};

  // print header
  std::cout << "   t  |";
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "R[" << setw(3) << i << "].x |";
    std::cout << "R[" << setw(3) << i << "].m |";
    std::cout << "R[" << setw(3) << i << "].V |";
  }
  std::cout << std::endl;

  // prep output file
  ofstream fout ("output.dat");

  // each step
  for (int t = 0; t < max_time_steps; ++t)
  {
    // each robot
    for (int i = 0; i < n_robots; ++i)
    {
      // decide control
      int i_next = (n_robots + i + 1) % n_robots;
      double dist_next = ((i == n_robots-1)? field_size: robots[i_next]) - robots[i];
      if (dist_next < 0)
        dist_next += field_size;
      int i_prev = (n_robots + i - 1) % n_robots;
      double dist_prev = robots[i] - ((i == 0)? 0: robots[i_prev]);
      if (dist_prev < 0)
        dist_prev += field_size;
      double u = ctrl_rate * (dist_next - dist_prev) / 2;
      if (u > ctrl_max)
        u = ctrl_max;

      // execute the control
      // motion noise
      std::normal_distribution<> motion_noise(0, sigmaM*std::fabs(u));
      // floating effect
      std::normal_distribution<> floating_noise(0, sigmaF * deltaT);
      // predict the next locaiton
      robots[i] += (u + motion_noise(gen)) * deltaT + floating_noise(gen);

      if (robots[i] < 0)
        robots[i] += field_size;
      else if (robots[i] >= field_size)
        robots[i] -= field_size;

      // predict the location
      // Note:
      //   mu = mu + u * DeltaT
      //   sigma = sigma + sigmaM * (u * DeltaT) ^ 2 + (sigF * DeltaT) ^ 2
      means[i] += u * deltaT;
      if (means[i] < 0)
        means[i] += field_size;
      else if (means[i] >= field_size)
        means[i] -= field_size;

      vars[i]
        += ((sigmaM * u) * (sigmaM * u) * deltaT * deltaT
         + sigmaF * sigmaF) * deltaT * deltaT;
    }

    std::vector<double> means_buff(means);
    std::vector<double> vars_buff(vars);

    std::vector<int> indx_list(n_robots);
    for (int i = 0; i < n_robots; ++i)
      indx_list[i] = i;
    if (random_order)
      std::random_shuffle(indx_list.begin(), indx_list.end());

    // each robot
    for (int indx = 0; indx < n_robots; ++indx)
    {
      if (prob_update)
      {
        // the estimates are updated at a specified probability.
        std::uniform_real_distribution<> dis(0, 1.0);
        double val = dis(gen);
        if (val > prob_update_p)
          continue;
      }

      // get the index for the robot to update.
      int i = indx_list[indx];

      // get the index of the next robot to communicate
      int j = (n_robots + i + 1) % n_robots;

      // allow the communication between the last and first robots only if
      // loopy updates are enabled.
      if (i == n_robots - 1 && loopy_updates == false)
        continue;

      // take measurement (communication)
      double z = robots[j] - robots[i];
      if (z < 0)
        z += field_size;
      std::normal_distribution<> sensor_noise{0,sigmaS};
      z += sensor_noise(gen);

      // update its location estimate
      double z_hat = means_buff[j] - means_buff[i];
      if (z_hat < 0)
        z_hat += field_size;
      double Ht1 = -1;
      double Ht2 = 1;
      double St1, St2;
      if (mode == 0)
      {
        St1 = Ht1*vars_buff[i]*Ht1 + sigmaS*sigmaS;
        St2 = Ht2*vars_buff[j]*Ht2 + sigmaS*sigmaS;
      }
      else if (mode == 1)
      {
        St1 = Ht1*vars_buff[i]*Ht1 + Ht2*vars_buff[j]*Ht2 + sigmaS*sigmaS;
        St2 = St1;
      }
      else if (mode == 2)
      {
        vars_buff[i] *= 2;
        vars_buff[j] *= 2;
        St1 = Ht1*vars_buff[i]*Ht1 + Ht2*vars_buff[j]*Ht2 + sigmaS*sigmaS;
        St2 = St1;
      }
      double K1 = vars_buff[i]*Ht1/St1;
      double K2 = vars_buff[j]*Ht2/St2;
      if (loopy_updates)
      {
        means_buff[i] += K1*(z - z_hat);
        if (means_buff[i] < 0)
          means_buff[i] += field_size;
        else if (means_buff[i] >= field_size)
          means_buff[i] -= field_size;
        vars_buff[i] *= (1 - K1*Ht1);
      }
      means_buff[j] += K2*(z - z_hat);
      if (means_buff[j] < 0)
        means_buff[j] += field_size;
      else if (means_buff[j] >= field_size)
        means_buff[j] -= field_size;
      vars_buff[j] *= (1 - K2*Ht2);
    }

    // if it is close to the origin, perform global localization
    if (global_loc && t % global_loc_steps == 0)
    {
      // take measurement (communication)
      double z = robots[0];
      std::normal_distribution<> global_loc_noise(0, sigmaGlobalLoc);
      z += global_loc_noise(gen);

      // update its location estimate
      double z_hat = means_buff[0];
      double Ht = 1;
      double St;
      St = Ht*vars_buff[0]*Ht + sigmaGlobalLoc*sigmaGlobalLoc;
      double K = vars_buff[0]*Ht/St;
      means_buff[0] += K*(z - z_hat);
      if (means_buff[0] < 0)
        means_buff[0] += field_size;
      else if (means_buff[0] >= field_size)
        means_buff[0] -= field_size;
      vars_buff[0] *= (1 - K*Ht);
    }

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
      errors[i] += std::fabs(means[i] - robots[i]);
    }

    // print the estimates
    std::cout << setw(6) << t << "|";
    fout << setw(6) << t << " ";
    for (int i = 0; i < n_robots; ++i)
    {
      std::cout << setw(9) << setprecision(3) << robots[i] << "|";
      std::cout << setw(9) << setprecision(3) << means[i] << "|";
      std::cout << setw(9) << setprecision(3) << std::sqrt(vars[i]) << "|";

      fout << setw(9) << setprecision(3) << robots[i] << " ";
      fout << setw(9) << setprecision(3) << means[i] << " ";
      fout << setw(9) << setprecision(3) << std::sqrt(vars[i]) << " ";
    }
    std::cout << std::endl;
    fout << std::endl;
  }
  fout.close();

  // output of gnuplot command
  std::cout << std::endl;
  std::cout << "gnuplot command" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    if (i == 0)
      std::cout << "plot ";
    else
      std::cout << "     ";
    std::cout << "\"output.dat\" u 1:" << std::to_string(2+i*3)
              << " title \"x" << std::to_string(1+i) << "\", \\"
              << std::endl;
  }
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "     ";
    std::cout << "\"output.dat\" u 1:"
              << std::to_string(3+i*3) << ":" << std::to_string(4+i*3)
              << " with errorbars title \"m" << std::to_string(1+i) << "\"";
    if (i < n_robots - 1)
      std::cout << ", \\";
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // display errors
  double total_error = 0;
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "robot[" << i << "]'s average error:" << (errors[i]/max_time_steps) << std::endl;
    total_error += errors[i];
  }
  std::cout << "overall average error: " << (total_error/max_time_steps/n_robots) << std::endl;


  // just for test of eigen
  Matrix3d A = Matrix3d::Identity();
  std::cout << "A = I: " << std::endl << A << std::endl;

  A(1,1) = 5;
  std::cout << "A(1,1) = 5: " << std::endl << A << std::endl;

  return 0;
}
