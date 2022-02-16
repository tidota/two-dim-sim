// simnonparam.cpp

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include <simnonparam.hpp>

// =============================================================================
SimNonParam::SimNonParam(const YAML::Node& doc): SimBase(doc),
  alpha1M(doc["alpha1M"].as<double>()),
  alpha2M(doc["alpha2M"].as<double>()),
  betaM(doc["betaM"].as<double>()),
  n_particles(doc["n_particles"].as<int>()),
  ests(n_robots),
  est_oris(n_robots),
  last_loc(n_robots),
  last_est(n_robots),
  last_est_ori(n_robots),
  cumul_errors(n_robots, 0),
  cumul_errors_phase(1, std::vector<double>(n_robots, 0)),
  // NOTE: time_steps_phase[0] is set to -1 so that step 0 is not counted.
  //(step 0 is the initial state and there is no error)
  time_steps_phase(1, -1),
  cumul_errors_ori(n_robots, 0),
  cumul_errors_ori_phase(1, std::vector<double>(n_robots, 0)),
  use_random_seed_pf(doc["use_random_seed_pf"].as<bool>()),
  random_seed_pf(doc["random_seed_pf"].as<unsigned int>()),
  gl_eval_cons(doc["gl_eval_cons"].as<double>()),
  plot_particles(doc["plot_particles"].as<bool>())
{
  // for all robots
  for(int i = 0; i < n_robots; ++i)
  {
    // init location x
    VectorXd buff(n_dim);
    for (int j = 0; j < n_dim; ++j)
    {
      buff(j) = doc["robots"][i][j].as<double>();
    }
    for (int j = 0; j < n_particles; ++j)
    {
      ests[i].push_back(buff);
      last_est[i].push_back(buff);
      if (use_orientation)
      {
        est_oris[i].push_back(doc["robots"][i][n_dim].as<double>());
        last_est_ori[i].push_back(doc["robots"][i][n_dim].as<double>());
      }
    }
  }

  // random value generator
  unsigned int seed;
  if (use_random_seed_pf)
  {
    seed = random_seed_pf;
  }
  else
  {
    std::random_device rd{};
    seed = rd();
  }
  gen_pf.seed(seed);
}

// =============================================================================
void SimNonParam::printSimInfo()
{
  SimBase::printSimInfo();

  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "average[" << i << "]: ";
    VectorXd mean = VectorXd::Zero(n_dim);
    for (int ip = 0; ip < n_particles; ++ip)
    {
      mean = mean + ests[i][ip];
    }
    mean = mean / n_particles;
    for (int idim = 0; idim < n_dim; ++idim)
    {
      std::cout << mean(idim) << ((idim == (n_dim - 1))? "": ", ");
    }
    std::cout << "| ";
    VectorXd var = VectorXd::Zero(n_dim);
    for (int ip = 0; ip < n_particles; ++ip)
    {
      VectorXd diff = ests[i][ip] - mean;
      for (int idim = 0; idim < n_dim; ++idim)
      {
        diff(idim) = diff(idim) * diff(idim);
      }
      var = var + diff;
    }
    var = var / n_particles;
    std::cout << "variance[" << i << "]: ";
    for (int idim = 0; idim < n_dim; ++idim)
    {
      std::cout << var(idim) << ((idim == (n_dim - 1))? "": ", ");
    }
    std::cout << std::endl;
  }
}

// =============================================================================
bool SimNonParam::startLog(const std::string& fname)
{
  // start log for nonparametric
  std::cout << "  t  |";
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "     R["
              << setw(2) << i
              << ((use_orientation)? "].x       |": "].x |");
    std::cout << "   Est["
              << setw(2) << i
              << ((use_orientation)? "].x       |": "].x |");
    // std::cout << "R[" << setw(7 * n_dim - 7) << i << "].m |";
    // std::cout << " det(v) |";
    std::cout << " err |";
  }
  std::cout << std::endl;

  // prep output file
  if (!fout)
    fout.close();
  fout.open(fname);

  return fout.good();
}

// =============================================================================
void SimNonParam::endLog()
{
  fout.close();

  // write the particles into a file at the end
  if (!fout_pf)
    fout_pf.close();
  fout_pf.open("particles.dat");
  for (int ip = 0; ip < n_particles; ++ip)
  {
    for (int irobot = 0; irobot < n_robots; ++irobot)
    {
      for (int idim = 0; idim < n_dim; ++idim)
      {
        fout_pf << this->last_est[irobot][ip](idim);
        fout_pf << " ";
      }
      if (use_orientation)
      {
        fout_pf << this->last_est_ori[irobot][ip];
        fout_pf << " ";
      }
    }
    fout_pf << std::endl;
  }
  fout_pf.close();

  // output of gnuplot command
  // NOTE: # of columns is different if using orientation.
  //       If not using orientation,
  //        actual position (x, y), estimated (x, y), error, p95.
  //       If using orientation,
  //        actual position (x, y, theta), estimated (x, y, theta),
  //        error of location, p95, error of orientation, reserved
  const int off_next_robot
    = (use_orientation)? n_dim + n_dim + 6: n_dim + n_dim + 2;
  const int off_next_column = (use_orientation)? n_dim + 1: n_dim;

  std::cout << std::endl;
  std::cout << "~~~ gnuplot command (errors vs uncertrainty) ~~~" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "--- ROBOT " << (i + 1) << " ---" << std::endl;
    std::cout << "clear" << std::endl;
    std::cout << "unset object" << std::endl;
    std::cout << "plot \"output.dat\" u 1:"
              << std::to_string(2 + i*off_next_robot + 2*off_next_column)
              << " title \"err of robot" << std::to_string(1+i) << "\" with line";
    std::cout << std::endl;
    std::cout << "replot \"output.dat\" u 1:"
              << std::to_string(2 + i*off_next_robot + 2*off_next_column + 1)
              << " title \"Uncertainty of robot" << std::to_string(1+i) << "\" with line";
    std::cout << std::endl;
    if (use_orientation)
    {
      std::cout << "--- orientation ---" << std::endl;
      std::cout << "clear" << std::endl;
      std::cout << "unset object" << std::endl;
      std::cout << "set yrange [-180:180]" << std::endl;
      std::cout << "plot \"output.dat\" u 1:"
                << std::to_string(2 + i*off_next_robot + 2*off_next_column + 2)
                << " title \"Error in Orientation (deg)" << std::to_string(1+i) << "\" with line";
      std::cout << std::endl;
    }
  }
  std::cout << std::endl;

  std::cout << "~~~ gnuplot command (errors as a whole) ~~~" << std::endl;
  std::cout << "clear" << std::endl;
  std::cout << "unset object" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    if (i == 0)
      std::cout << "plot ";
    else
      std::cout << "     ";
    std::cout << "\"output.dat\" u 1:"
              << std::to_string(2 + i*off_next_robot + 2*off_next_column)
              << " title \"err" << std::to_string(1+i) << "\" with line";
    if (i < n_robots - 1)
      std::cout << ", \\";
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "~~~ (orientation errors) ~~~" << std::endl;
  std::cout << "clear" << std::endl;
  std::cout << "unset object" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    if (i == 0)
      std::cout << "plot ";
    else
      std::cout << "     ";
    std::cout << "\"output.dat\" u 1:"
              << std::to_string(2 + i*off_next_robot + 2*off_next_column + 2)
              << " title \"orientation err" << std::to_string(1+i) << " (deg)\" with line";
    if (i < n_robots - 1)
      std::cout << ", \\";
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "~~~ trajectories of the robots (ground-truth vs estimation) ~~~"
            << std::endl;
  std::cout << "writing a gnuplot script..." << std::endl;
  std::fstream f_gnuplot("gnuplot_traj.plt", std::fstream::out);
  f_gnuplot << "clear" << std::endl;

  // --- particles at the end --- //
  f_gnuplot << "unset object" << std::endl;
  if (show_covs)
  {
    int count_cov = 0;
    for (unsigned int i = 0; i < cov_buff.size(); ++i)
    {
      // if (i % 10 > 0)
      //   continue;
      for (auto data: cov_buff[i])
      {
        double s1 = std::sqrt(data[2]) * 2 * std::sqrt(5.991);
        double s2 = std::sqrt(data[3]) * 2 * std::sqrt(5.991);
        f_gnuplot << "set object " << std::to_string(count_cov + 1)
                  << " ellipse center "
                  << std::to_string(data[0]) << ","
                  << std::to_string(data[1])
                  << " size "
                  // https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
                  << ((std::isnan(s1))? "0.0": std::to_string(s1)) << ","
                  << ((std::isnan(s2))? "0.0": std::to_string(s2))
                  << " angle "
                  << std::to_string(data[4])
                  << " front fillstyle empty border -1" << std::endl;
        ++count_cov;
      }
    }
  }
  f_gnuplot << "h1 = 227/360.0" << std::endl;
  f_gnuplot << "h2 = 40/360.0" << std::endl;
  f_gnuplot << "set palette model HSV functions (1-gray)*(h2-h1)+h1,1,0.68"
            << std::endl;
  f_gnuplot << "set size ratio -1" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    if (i == 0)
      f_gnuplot << "plot   ";
    else
      f_gnuplot << "replot ";
    f_gnuplot << "\"output.dat\" u ";
    for (int j = 0; j < n_dim; ++j)
    {
      f_gnuplot << std::to_string(2+i*off_next_robot+j);
      if (j < n_dim - 1)
        f_gnuplot << ":";
      else
        f_gnuplot << ":1";
    }
    f_gnuplot << " title \"R" << std::to_string(1+i) << "\""
              << " with linespoints dt 2 lt -1 lw 2.0 ps 3.0"
              << " pt " << std::to_string(i + 1)
              << " lc palette" << std::endl;
  }
  for (int i = 0; i < n_robots; ++i)
  {
    f_gnuplot << "replot ";
    f_gnuplot << "\"output.dat\" u ";
    for (int j = 0; j < n_dim; ++j)
    {
      f_gnuplot << std::to_string(2+off_next_column+i*off_next_robot+j);
      if (j < n_dim - 1)
        f_gnuplot << ":";
      else
        f_gnuplot << ":1";
    }
    f_gnuplot << " title \"Est" << std::to_string(1+i) << "\""
              << " with linespoints lt 1 lw 3.0 ps 3.0"
              << " pt " << std::to_string(i + 1)
              << " lc palette" << std::endl;
  }
  for (int i = 0; i < n_robots; ++i)
  {
    f_gnuplot << "replot ";
    f_gnuplot << "\"particles.dat\" u ";
    for (int j = 0; j < n_dim; ++j)
    {
      f_gnuplot << std::to_string(1+i*off_next_column+j);
      if (j < n_dim - 1)
        f_gnuplot << ":";
    }
    f_gnuplot << " notitle ps 1.0 "
              << " pt " << std::to_string(i + 1)
              << " lc black" << std::endl;
  }
  if (use_orientation)
  {
    f_gnuplot << "xf(phi)=0.2*cos(phi)" << std::endl;
    f_gnuplot << "yf(phi)=0.2*sin(phi)" << std::endl;
    for (int i = 0; i < n_robots; ++i)
    {
      f_gnuplot << "replot ";
      f_gnuplot << "\"particles.dat\" u ";
      for (int j = 0; j < n_dim; ++j)
      {
        f_gnuplot << std::to_string(1+i*off_next_column+j);
        if (j < n_dim - 1)
          f_gnuplot << ":";
      }
      f_gnuplot << ":(xf($" << std::to_string(1+i*off_next_column+n_dim) << "))"
                << ":(yf($" << std::to_string(1+i*off_next_column+n_dim) << "))"
                << " notitle with vectors head size 0.1,20,60 filled"
                << " lc black" << std::endl;
    }
  }
  f_gnuplot << "pause -1 \"Hit any key to continue\"" << std::endl;
  f_gnuplot.close();
  std::cout << "To run the gnuplot script:" << std::endl << std::endl;
  std::cout << "        gnuplot gnuplot_traj.plt" << std::endl << std::endl;

  std::cout << "writing a matplotlib.pyplot script..." << std::endl;
  std::fstream f_pyplot("pyplot_traj.py", std::fstream::out);
  f_pyplot <<
    "import numpy as np" << std::endl <<
    "import matplotlib.pyplot as plt" << std::endl <<
    "from matplotlib.lines import Line2D" << std::endl <<
    "import colorsys" << std::endl <<
    "import matplotlib as mpl" << std::endl <<
    "from mpl_toolkits.axes_grid1 import make_axes_locatable" << std::endl <<
    "import copy" << std::endl <<
    "import matplotlib.patches as mpatches" << std::endl <<
    "import math" << std::endl <<

    "mpl.rcParams['text.usetex'] = True" << std::endl <<

    "cm_new = True" << std::endl <<

    "# load data" << std::endl <<
    "f = open('output.dat', 'r')" << std::endl <<
    "data = np.array([[float(v) for v in line.split()] for line in f])" << std::endl <<
    "f.close()" << std::endl <<

    "maxt = data[-1][0]" << std::endl <<

    "def f(t):" << std::endl <<
    "    return np.exp(-t) * np.cos(2*np.pi*t)" << std::endl <<
    "def g(t):" << std::endl <<
    "    return np.exp(-2*t) * np.sin(2*np.pi*t)" << std::endl <<

    "def cmf(val):" << std::endl <<
    "    if cm_new:" << std::endl <<
    "        val /= maxt" << std::endl <<
    "        return plt.cm.jet(0.97-val*0.82)" << std::endl <<
    "    else:" << std::endl <<
    "        val /= maxt" << std::endl <<
    "        h1 = 227/360" << std::endl <<
    "        h2 = 40/360" << std::endl <<
    "        h = (1-val)*(h2-h1)+h1" << std::endl <<
    "        return np.array(list(colorsys.hsv_to_rgb(h,1,0.68)) + [1])" << std::endl <<

    "subtle_alpha = 0.90" << std::endl <<
    "def cmf_s(val):" << std::endl <<
    "    if cm_new:" << std::endl <<
    "        val /= maxt" << std::endl <<
    "        c = plt.cm.jet(0.97-val*0.82)" << std::endl <<
    //"        c[3] = subtle_alpha" << std::endl <<
    "        return c" << std::endl <<
    "    else:" << std::endl <<
    "        val /= maxt" << std::endl <<
    "        h1 = 227/360" << std::endl <<
    "        h2 = 40/360" << std::endl <<
    "        h = (1-val)*(h2-h1)+h1" << std::endl <<
    "        return np.array(list(colorsys.hsv_to_rgb(h,1,0.68)) + [subtle_alpha])" << std::endl <<

    "# def cm1f(val): # original" << std::endl <<
    "#     val /= maxt" << std::endl <<
    "#     h1 = 227/360" << std::endl <<
    "#     h2 = 40/360" << std::endl <<
    "#     h = (1-val)*(h2-h1)+h1" << std::endl <<
    "#     return np.array(list(colorsys.hsv_to_rgb(h,1,0.68)) + [1])" << std::endl <<

    "# def cm2f(val): # kind of new" << std::endl <<
    "#     val /= maxt" << std::endl <<
    "#     return plt.cm.jet(0.97-val*0.82)" << std::endl <<

    "def plotPath(points, style, cm, lw=1, ms=5):" << std::endl <<
    "    # p[0]: x_1" << std::endl <<
    "    # p[1]: x_2" << std::endl <<
    "    # p[2]: t" << std::endl <<
    "    for i in range(1,len(points)):" << std::endl <<
    "        p_prev = points[i-1]" << std::endl <<
    "        p = points[i]" << std::endl <<
    "        plt.plot(" << std::endl <<
    "            [p_prev[0], p[0]]," << std::endl <<
    "            [p_prev[1], p[1]]," << std::endl <<
    "            style if '-' in style else '-' + style," << std::endl <<
    "            color=cm(p[2])," << std::endl <<
    "            linewidth=lw," << std::endl <<
    "            markersize=ms)" << std::endl <<

    "# custom colormap" << std::endl <<
    "cm = copy.deepcopy(plt.cm.get_cmap('jet'))" << std::endl <<
    "cm._init()" << std::endl <<
    "cm._lut = np.array([cmf(x) for x in np.linspace(0, maxt, cm.N)])" << std::endl <<
    "mpl.cm.register_cmap(name='cm', cmap=cm)" << std::endl <<

    "cm_s = copy.deepcopy(plt.cm.get_cmap('jet'))" << std::endl <<
    "cm_s._init()" << std::endl <<
    "cm_s._lut = np.array([cmf_s(x) for x in np.linspace(0, maxt, cm_s.N)])" << std::endl <<
    "mpl.cm.register_cmap(name='cm_s', cmap=cm_s)" << std::endl <<

    "# prep figure" << std::endl <<
    "fig, ax = plt.subplots()" << std::endl <<
    "#fig.set_size_inches(6,5)" << std::endl <<
    "fig.set_dpi(200)" << std::endl <<
    "# ax.set_aspect('equal')" << std::endl <<
    "#ax.set_title('Trajectories')" << std::endl <<

    "# markers" << std::endl <<
    "markers = [" << std::endl <<
    "    'o', '*', '^', 's', '1', 'x', 'v', 'D', '1', '+'" << std::endl <<
    "]" << std::endl <<

    "# plot ground-truth" << std::endl <<
    "for i in range(8):" << std::endl <<
    "    print('plotting ground-truth of robot ' + str(i+1) + '...')" << std::endl <<
    "    plotPath(" << std::endl <<
    "        np.concatenate(" << std::endl <<
    "            (data[:,1+" << off_next_robot << "*i:2+" << off_next_robot << "*i]," << std::endl <<
    "             data[:,2+" << off_next_robot << "*i:3+" << off_next_robot << "*i]," << std::endl <<
    "             data[:,0:1]),axis=1).tolist()," << std::endl <<
    "        markers[i]+'--', cmf_s, lw=1, ms=4)" << std::endl <<
    "# plot estimated trajectories" << std::endl <<
    "for i in range(8):" << std::endl <<
    "    print('plotting estimation of robot ' + str(i+1) + '...')" << std::endl <<
    "    plotPath(" << std::endl <<
    "        np.concatenate(" << std::endl <<
    "            (data[:,1+" << off_next_column << "+" << off_next_robot << "*i:2+"
                             << off_next_column << "+" << off_next_robot << "*i]," << std::endl <<
    "             data[:,2+" << off_next_column << "+" << off_next_robot << "*i:3+"
                             << off_next_column << "+" << off_next_robot << "*i]," << std::endl <<
    "             data[:,0:1]),axis=1).tolist()," << std::endl <<
    "        markers[i], cmf, lw=1, ms=4)" << std::endl;

  if (plot_particles)
  {
    f_pyplot <<
      "#plot particles" << std::endl <<
      "f = open('particles.dat', 'r')" << std::endl <<
      "pdata = np.array([[float(v) for v in line.split()] for line in f])" << std::endl <<
      "f.close()" << std::endl <<
      "for i in range(" << n_robots << "):" << std::endl <<
      "    print('plotting particles of robot ' + str(i+1) + '...')" << std::endl <<
      "    plt.plot(" << std::endl <<
      "        pdata[:,0+" << off_next_column << "*i:1+" << off_next_column << "*i].tolist(), pdata[:,1+" << off_next_column << "*i:2+" << off_next_column << "*i].tolist()," << std::endl <<
      "        markers[i], color='black', linewidth=0.5, markersize=1, zorder=100)" << std::endl;

    if(use_orientation)
    {
      f_pyplot <<
        "n_particles = len(pdata[:,0:1].tolist())" << std::endl <<
        "for i in range(" << n_robots << "):" << std::endl <<
        "    print('plotting arrows of robot ' + str(i+1) + '...')" << std::endl <<
        "    for j in range(n_particles):" << std::endl <<
        "        plt.arrow(" << std::endl <<
        "            pdata[j," << off_next_column << "*i], pdata[j,1+" << off_next_column << "*i]," << std::endl <<
        "            0.1*math.cos(pdata[j,2+" << off_next_column << "*i]), 0.1*math.sin(pdata[j,2+" << off_next_column << "*i])," <<
        "            head_width=0.05, head_length=0.05, fill=True, color='black', zorder=100)" << std::endl;
    }
  }

  f_pyplot <<
    "# fix the ratio" << std::endl <<
    "xleft, xright = ax.get_xlim()" << std::endl <<
    "ybottom, ytop = ax.get_ylim()" << std::endl <<
    "ratio = abs((xright-xleft)/(ybottom-ytop))" << std::endl <<
    "ax.set_aspect(1)" << std::endl <<

    "# labels" << std::endl <<
    "plt.ylabel(r'$x_2$ \\small{(meters)}', fontsize=18)" << std::endl <<
    "plt.xlabel(r'$x_1$ \\small{(meters)}', fontsize=18)" << std::endl <<

    "# colorbar" << std::endl <<
    "print('adding colorbar...')" << std::endl <<
    "divider = make_axes_locatable(plt.gca())" << std::endl <<
    "norm = mpl.colors.Normalize(vmin = 0, vmax = maxt)" << std::endl <<
    "ax_cb = divider.new_horizontal(size=\"5%\", pad=0.05)" << std::endl <<
    "ax_cb.set_title('$t$ \\small{(sec)}', fontsize=18)" << std::endl <<
    "cb = mpl.colorbar.ColorbarBase(ax_cb, cmap=cm, norm=norm, orientation='vertical')" << std::endl <<
    "plt.gcf().add_axes(ax_cb)" << std::endl <<

    "# legends" << std::endl <<
    "print('adding legends...')" << std::endl <<
    "# legend_elements = [" << std::endl <<
    "#     Line2D([0], [0]," << std::endl <<
    "#         color='k', lw=1," << std::endl <<
    "#         linestyle='--'," << std::endl <<
    "#         marker=markers[i]," << std::endl <<
    "#         markersize=5, label='R'+str(i+1)) for i in range(8)]" << std::endl <<
    "# legend_elements += [" << std::endl <<
    "#     Line2D([0], [0]," << std::endl <<
    "#         color='k', lw=1," << std::endl <<
    "#         marker=markers[i]," << std::endl <<
    "#         markersize=8, label='Est'+str(i+1)) for i in range(8)]" << std::endl <<
    "legend_elements = [" << std::endl <<
    "    Line2D([0], [0]," << std::endl <<
    "        color='k', lw=1," << std::endl <<
    "        linestyle=''," << std::endl <<
    "        marker=markers[i]," << std::endl <<
    "        markersize=5, label='Robot'+str(i+1)) for i in range(8)]" << std::endl <<
    "ax.legend(" << std::endl <<
    "    handles=legend_elements," << std::endl <<
    "    loc='best'," << std::endl <<
    "    fontsize=5)" << std::endl;

  f_pyplot <<
    "# show the figure" << std::endl <<
    "print('showing the figure...')" << std::endl <<
    "plt.show()" << std::endl <<

    "# save the figure" << std::endl <<
    "#fig.savefig('result.png')" << std::endl;
  f_pyplot.close();
  std::cout << "To run the matplot.pyplot script: " << std::endl << std::endl;
  std::cout << "        python3 pyplot_traj.py" << std::endl << std::endl;

  std::cout << "~~~ average errors ~~~" << std::endl;
  // display errors
  {
    std::cout << "time steps: " << (max_time*sim_freq) << std::endl;
    double total_error = 0;
    for (int i = 0; i < n_robots; ++i)
    {
      std::cout << "robot[" << i << "]'s average error:" << (cumul_errors[i]/(max_time*sim_freq)) << std::endl;
      total_error += cumul_errors[i];
    }
    std::cout << "overall average error: " << (total_error/(max_time*sim_freq)/n_robots) << std::endl;
  }
  for (unsigned int iphase = 0; iphase < cumul_errors_phase.size(); ++iphase)
  {
    std::cout << "~~~ average errors (phase " << (iphase+1) << ") ~~~" << std::endl;
    // display errors
    {
      std::cout << "time steps: " << time_steps_phase[iphase] << std::endl;
      double total_error_phase = 0;
      for (int i = 0; i < n_robots; ++i)
      {
        std::cout << "robot[" << i << "]'s average error:"
                  << (cumul_errors_phase[iphase][i]/time_steps_phase[iphase]) << std::endl;
        total_error_phase += cumul_errors_phase[iphase][i];
      }
      std::cout << "overall average error: "
                << (total_error_phase/time_steps_phase[iphase]/n_robots) << std::endl;
    }
  }
  std::cout << std::endl;

  if (use_orientation)
  {
    std::cout << "~~~ average errors (ori) ~~~" << std::endl;
    // display errors in orientation
    {
      std::cout << "time steps: " << (max_time*sim_freq) << std::endl;
      double total_error_ori = 0;
      for (int i = 0; i < n_robots; ++i)
      {
        std::cout << "robot[" << i << "]'s average error (ori in deg):"
                  << (cumul_errors_ori[i]/(max_time*sim_freq) / M_PI * 180.0) << std::endl;
        total_error_ori += cumul_errors_ori[i];
      }
      std::cout << "overall average error (ori in deg): "
                << (total_error_ori/(max_time*sim_freq)/n_robots / M_PI * 180.0) << std::endl;
    }
    for (unsigned int iphase = 0; iphase < cumul_errors_ori_phase.size(); ++iphase)
    {
      std::cout << "~~~ average errors (ori, phase " << (iphase+1) << ") ~~~" << std::endl;
      // display errors
      {
        std::cout << "time steps: " << time_steps_phase[iphase] << std::endl;
        double total_error_ori_phase = 0;
        for (int i = 0; i < n_robots; ++i)
        {
          std::cout << "robot[" << i << "]'s average error (ori in deg):"
                    << (cumul_errors_ori_phase[iphase][i]/time_steps_phase[iphase] / M_PI * 180.0) << std::endl;
          total_error_ori_phase += cumul_errors_ori_phase[iphase][i];
        }
        std::cout << "overall average error (ori in deg): "
                  << (total_error_ori_phase/time_steps_phase[iphase]/n_robots / M_PI * 180.0) << std::endl;
      }
    }
  }
}

// =============================================================================
void SimNonParam::plotImpl()
{
  // A line contains each robot's information.
  // - the actual location
  // - the average estimated location
  // - errors

  // Time
  std::cout << std::fixed << std::setprecision(2);
  std::cout << std::right << std::setw(5) << t << "|";
  fout << std::fixed << std::setprecision(3);
  fout << std::right << std::setw(8) << t << " ";

  std::vector< std::vector <double> > covs_at_time;
  // For each robot
  for (int i = 0; i < n_robots; ++i)
  {
    // calculate average estimated location
    VectorXd average = VectorXd::Zero(n_dim);
    // NOTE: calculating the average orientation by using vectors
    //double ave_ori = 0;
    double sum_sin = 0;
    double sum_cos = 0;
    VectorXd squared_average = VectorXd::Zero(n_dim);
    double mix_average = 0;
    for (int ip = 0; ip < n_particles; ++ip)
    {
      average += this->ests[i][ip];
      if (use_orientation)
      {
        sum_sin += std::sin(this->est_oris[i][ip]);
        sum_cos += std::cos(this->est_oris[i][ip]);
        // if (ave_ori >= this->est_oris[i][ip])
        // {
        //   if (ave_ori - this->est_oris[i][ip] <= M_PI)
        //     ave_ori = (ip * ave_ori + this->est_oris[i][ip])/(ip + 1);
        //   else
        //   {
        //     ave_ori
        //       = (ip * (ave_ori - 2*M_PI) + this->est_oris[i][ip])/(ip + 1);
        //     if (ave_ori <= -M_PI)
        //       ave_ori += 2*M_PI;
        //   }
        // }
        // else
        // {
        //   if (this->est_oris[i][ip] - ave_ori <= M_PI)
        //     ave_ori = (this->est_oris[i][ip] + ip * ave_ori)/(ip + 1);
        //   else
        //   {
        //     ave_ori
        //       = ((this->est_oris[i][ip] - 2*M_PI) + ip * ave_ori)/(ip + 1);
        //     if (ave_ori <= -M_PI)
        //       ave_ori += 2*M_PI;
        //   }
        // }
      }
      if (show_covs)
      {
        VectorXd sqrd(2);
        sqrd << this->ests[i][ip][0] * this->ests[i][ip][0],
                this->ests[i][ip][1] * this->ests[i][ip][1];
        squared_average += sqrd;
        mix_average += this->ests[i][ip][0] * this->ests[i][ip][1];
      }
    }
    // for covariance
    average /= n_particles;
    double p95 = 0;
    if (show_covs)
    {
      squared_average /= n_particles;
      mix_average /= n_particles;
      MatrixXd vars(2,2);
      vars << squared_average[0] - average[0]*average[0],
              mix_average - average[0]*average[1],
              mix_average - average[0]*average[1],
              squared_average[1] - average[1]*average[1];
      Eigen::EigenSolver<MatrixXd> s(vars);
      auto eigen_val = s.eigenvalues();
      auto eigen_vec = s.eigenvectors();
      double var_ang
        = std::atan2(eigen_vec.col(0)[1].real(), eigen_vec.col(0)[0].real())
          / M_PI*180.0;
      std::vector<double> cov
        = {average[0], average[1],
           eigen_val[0].real(), eigen_val[1].real(),
           var_ang};
      covs_at_time.push_back(std::move(cov));

      VectorXd v = (robots[i] - average);
      Rotation2D<double> Rot(
          -std::atan2(eigen_vec.col(0)[1].real(), eigen_vec.col(0)[0].real()));
      v = Rot * v;
      double a = v(0);
      double b = v(1);
      double sx = std::sqrt(eigen_val[0].real());
      double sy = std::sqrt(eigen_val[1].real());
      if (vars.determinant() > 0 && v.norm() > 0)
      {
        p95 = std::sqrt(5.991) * sx * sy
            * std::sqrt((a*a+b*b)/(a*a*sy*sy + b*b*sx*sx));
      }
    }
    // calc average orientation
    double ave_ori;
    if (use_orientation)
    {
      assert(sum_sin + sum_cos != 0);
      ave_ori = std::atan2(sum_sin, sum_cos);
    }

    // --- Terminal Output --- //
    // current location
    for (int j = 0; j < n_dim; ++j)
    {
      std::cout << std::right << std::setw(6) << robots[i](j);
      if (j < n_dim - 1 || use_orientation)
        std::cout << ",";
      else
        std::cout << "|";
    }
    if (use_orientation)
    {
      std::cout << std::right << std::setw(5) << oris[i] << "|";
    }
    // average estimation
    for (int j = 0; j < n_dim; ++j)
    {
      std::cout << std::right << std::setw(6) << average(j);
      if (j < n_dim - 1 || use_orientation)
        std::cout << ",";
      else
        std::cout << "|";
    }
    if (use_orientation)
    {
      std::cout << std::right << std::setw(5) << ave_ori << "|";
    }
    // error
    std::cout << std::right << std::setw(5) << errors[i] << "|";

    // --- File Output --- //
    // current location
    for (int j = 0; j < n_dim; ++j)
    {
      fout << std::right << std::setw(9) << robots[i](j);
      fout << " ";
    }
    if (use_orientation)
    {
      fout << std::right << std::setw(9) << oris[i];
      fout << " ";
    }
    // average estimation
    for (int j = 0; j < n_dim; ++j)
    {
      fout << std::right << std::setw(9) << average(j);
      fout << " ";
    }
    if (use_orientation)
    {
      fout << std::right << std::setw(9) << ave_ori;
      fout << " ";
    }
    // error
    fout << std::right << std::setw(8)
         << (robots[i] - average).norm() << " ";
    // length at the ellipse counter crosses the lin to the ground-truth point.
    fout << std::right << std::setw(8)
         << p95 << " ";
    if (use_orientation)
    {
      double diff_ori = ave_ori - oris[i];
      if (diff_ori > M_PI)
        diff_ori -= 2*M_PI;
      if (diff_ori <= -M_PI)
        diff_ori += 2*M_PI;
      diff_ori = diff_ori / M_PI * 180.0;
      fout << std::right << std::setw(8)
           << diff_ori << " ";
      fout << "0 ";
    }

    // store data to the buffer to plot
    last_loc[i] = robots[i];
    for (int ip = 0; ip < n_particles; ++ip)
    {
      last_est[i][ip] = this->ests[i][ip];
      if (use_orientation)
        last_est_ori[i][ip] = this->est_oris[i][ip];
    }
  }
  cov_buff.push_back(std::move(covs_at_time));

  std::cout << std::endl;
  fout << std::endl;
}

// =============================================================================
void SimNonParam::predict()
{
  // Prediction based on the motion model.
  // === prediction ===
  // for all robots
  for (int i = 0; i < n_robots; ++i)
  {
    if (use_orientation && use_vel_ctrl)
    {
      std::normal_distribution<> rnd1(0, std::sqrt(alpha1M) * vels_rot[i](0));
      std::normal_distribution<> rnd2(0, std::sqrt(alpha2M) * vels_rot[i](0));
      std::normal_distribution<> rndO(0, vels_rot[i](1)/5);
      std::normal_distribution<> rndf(0, std::sqrt(betaM));
      //vels_rot[i];
      for (int ip = 0; ip < n_particles; ++ip)
      {
        // convert it into Cartesian coordinates by the particle's theta
        VectorXd vel_cart = VectorXd::Zero(2);
        vel_cart(0) = vels_rot[i](0) * std::cos(est_oris[i][ip]);
        vel_cart(1) = vels_rot[i](0) * std::sin(est_oris[i][ip]);

        // add noise
        VectorXd noise1_x1(n_dim);
        noise1_x1 << std::cos(est_oris[i][ip]), std::sin(est_oris[i][ip]);
        noise1_x1 *= rnd1(gen_pf);
        VectorXd noise1_x2(n_dim);
        noise1_x2 << std::cos(est_oris[i][ip] + M_PI/2),
                     std::sin(est_oris[i][ip] + M_PI/2);
        noise1_x2 *= rnd2(gen_pf);
        VectorXd noise2(n_dim);
        noise2 << rndf(gen_pf), rndf(gen_pf);

        // calcualte deltaX
        VectorXd deltaX = deltaT * (vel_cart + noise1_x1 + noise1_x2 + noise2);
        // update the particle
        ests[i][ip] += deltaX;

        // calculate deltaTheta
        double deltaTheta = deltaT * (vels_rot[i](1) + rndO(gen_pf));
        est_oris[i][ip] += deltaTheta;
        if (est_oris[i][ip] > M_PI)
          est_oris[i][ip] -= 2*M_PI;
        else if (est_oris[i][ip] <= -M_PI)
          est_oris[i][ip] += 2*M_PI;
      }
    }
    else
    {
      // velocity based on velocity control input
      VectorXd deltaX = MatrixXd::Identity(n_dim, n_dim) * vels[i] * deltaT;

      // === motion model === //
      // calculate eigen vectors and values
      double v = vels[i].norm();
      MatrixXd EigenVecs(2, 2);
      EigenVecs(0, 0) = vels[i](0) / v;
      EigenVecs(1, 0) = vels[i](1) / v;
      EigenVecs(0, 1) = -vels[i](1) / v;
      EigenVecs(1, 1) = vels[i](0) / v;

      std::normal_distribution<> rnd1(0, std::sqrt(alpha1M) * v);
      std::normal_distribution<> rnd2(0, std::sqrt(alpha2M) * v);
      std::normal_distribution<> rndf(0, std::sqrt(betaM));

      for (int ip = 0; ip < n_particles; ++ip)
      {
        // get a random value based on the first eigen value
        double x1 = rnd1(gen_pf);

        // get a random value based on the second eigen value
        double x2 = rnd2(gen_pf);

        // make a vector of these random values
        VectorXd rnd2D(n_dim);
        rnd2D << x1, x2;

        // transform it by eigen vectors
        //EigenVecs * vec
        VectorXd noise1 = EigenVecs * rnd2D;

        // get a random vector based on the float effect
        // VectorXd FloatEffect = VectorXd::Ones(2) * betaM;
        double x1f = rndf(gen_pf);
        double x2f = rndf(gen_pf);
        VectorXd noise2(n_dim);
        noise2 << x1f, x2f;

        // calculate velocity with the noises
        // add it to the estimation
        ests[i][ip] += deltaX + (noise1 + noise2) * deltaT;
      }
    }
  }
}

// =============================================================================
void SimNonParam::globalLocImpl(const VectorXd& z)
{
  // z[0]: distance from the origin to the first robot
  // z[1]: driection toward the first robot

  // create cumulative weights for resampling
  std::vector<double> cumul_weights(n_particles, 0);

  // for all particles of the first robot
  for (int i = 0; i < n_particles; ++i)
  {
    VectorXd z_hat(2);
    z_hat(0) = ests[0][i].norm();

    if (use_beacon_sensor)
    {
      z_hat(1) = std::atan2(-ests[0][i](1), -ests[0][i](0)) - est_oris[0][i];
      if (z_hat(1) > M_PI)
      {
        z_hat(1) -= 2*M_PI;
      }
      else if (z_hat(1) <= -M_PI)
      {
        z_hat(1) += 2*M_PI;
      }
    }
    else
    {
      z_hat(1) = std::atan2(ests[0][i](1), ests[0][i](0));
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

    // evaluate
    cumul_weights[i]
      = exp(
          -z_diff(0)*z_diff(0)/sigmaGlobalLocR/sigmaGlobalLocR/gl_eval_cons
          -z_diff(1)*z_diff(1)/sigmaGlobalLocT/sigmaGlobalLocT/gl_eval_cons);

    if (i > 0)
      cumul_weights[i] += cumul_weights[i - 1];
  }

  // new  population
  std::vector<VectorXd> new_est(n_particles, VectorXd(n_dim));
  std::vector<double> new_est_ori(n_particles);

  // resample
  for (int i = 0; i < n_particles; ++i)
  {
    // decide the index to pick up
    int indx = drawRandIndx(cumul_weights);

    // add the picked one
    new_est[i] = this->ests[0][indx];
    if (use_orientation)
      new_est_ori[i] = this->est_oris[0][indx];
  }

  // swap
  std::swap(this->ests[0], new_est);
  if (use_orientation)
    std::swap(this->est_oris[0], new_est_ori);
}

// =============================================================================
void SimNonParam::calcErrors()
{
  // calculate errors
  for (int i = 0; i < n_robots; ++i)
  {
    VectorXd mean = VectorXd::Zero(n_dim);
    for (int ip = 0; ip < n_particles; ++ip)
    {
      mean += ests[i][ip];
    }
    mean /= n_particles;
    errors[i] = (mean - robots[i]).norm();
    cumul_errors[i] += errors[i];
    cumul_errors_phase[cumul_errors_phase.size() - 1][i] += errors[i];
  }
  // NOTE: calculate errors in orientation
  if (use_orientation)
  {
    for (int i = 0; i < n_robots; ++i)
    {
      double sum_sin = 0;
      double sum_cos = 0;
      for (int ip = 0; ip < n_particles; ++ip)
      {
        sum_sin += std::sin(this->est_oris[i][ip]);
        sum_cos += std::cos(this->est_oris[i][ip]);
      }
      double ave_ori;
      assert(sum_sin + sum_cos != 0);
      ave_ori = std::atan2(sum_sin, sum_cos);

      double diff_ori = ave_ori - oris[i];
      if (diff_ori > M_PI)
        diff_ori -= 2*M_PI;
      if (diff_ori <= -M_PI)
        diff_ori += 2*M_PI;
      double error_ori = std::fabs(diff_ori);
      cumul_errors_ori[i] += error_ori;
      cumul_errors_ori_phase[cumul_errors_ori_phase.size() - 1][i] += error_ori;
    }
  }
  time_steps_phase[time_steps_phase.size() - 1]++;
}

// =============================================================================
int SimNonParam::drawRandIndx(const std::vector<double>& cumul_weights)
{
  std::uniform_real_distribution<>
    dist(0, cumul_weights[cumul_weights.size() - 1]);

  double val = dist(gen_pf);
  int lo = 0;
  int hi = cumul_weights.size() - 1;

  while (lo < hi - 1)
  {
    int mid = (hi + lo)/2;

    if (val > cumul_weights[mid])
    {
      lo = mid;
    }
    else
    {
      hi = mid;
    }
  }
  int ans = hi;
  if (val <= cumul_weights[lo])
    ans = lo;
  return ans;
}
