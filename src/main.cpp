/*

The main part of the simulation.

*/

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simnaive.hpp>
#include <simcons.hpp>
#include <simci.hpp>
#include <simcent.hpp>
#include <simpfcons.hpp>

using namespace std;
using namespace Eigen;

void getOmega(MatrixXd &, MatrixXd &, double &, const int &);

int main (int argc, char** argv)
{
  std::ifstream fin("settings.yaml");
  YAML::Node doc = YAML::Load(fin);

  //=========================================================================================================
  // in the main function

  // mode
  const int mode = doc["mode"].as<int>();

  // second destinations (if it is true, the robot intervals will be doubled in
  // the middle of operation)
  const bool second_dest = doc["second_dest"].as<bool>();
  const double second_dest_rate = doc["second_dest_rate"].as<double>();
  bool second_dest_udpated = false;

  // instantiate the simulation object
  std::shared_ptr<SimBase> sim;

  if (mode == 0)
    sim = std::make_shared<SimNaiveOrig>(doc);
  else if (mode == 1)
    sim = std::make_shared<SimNaive>(doc);
  else if (mode == 2)
    sim = std::make_shared<SimConsOrig>(doc);
  else if (mode == 3)
    sim = std::make_shared<SimCons>(doc);
  else if (mode == 4)
    sim = std::make_shared<SimCi>(doc);
  else if (mode == 5)
    sim = std::make_shared<SimCent>(doc);
  else if (mode == 6)
    sim = std::make_shared<SimPfCons>(doc);
  else
  {
    std::cout << "Invalid mode number: " << mode << std::endl;
    return -1;
  }

  fin.close();

  // print sim infor
  sim->printSimInfo();

  // start log from parametric
  sim->startLog("output.dat");

  // main part of simulation
  for (; !sim->isDone(); sim->stepForward())
  {
    if (second_dest && !second_dest_udpated && sim->getProgress() > 0.5)
    {
      std::cout << "SECOND DESTINATIONS: thresholds are being doubled."
                << std::endl;
      sim->scaleControl(second_dest_rate, second_dest_rate);
      second_dest_udpated = true;
    }

    // print the current status
    sim->plot();

    // update sim
    sim->updateSim();

    // predict by the motion model
    sim->predict();

    // global localization
    sim->globalLoc();

    // mutual localization
    sim->mutualLoc();

    // calcError from parametric
    sim->calcErrors();
  }

  // endLog from parametric
  sim->endLog();

  return 0;
}
