// simnonparam.cpp

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
  last_loc(n_robots),
  last_est(n_robots)
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
    std::vector<VectorXd> est;
    for (int j = 0; j < n_particles; ++j)
    {
      est.push_back(buff);
    }
    ests.push_back(est);
  }
}

// =============================================================================
void SimNonParam::printSimInfo()
{
  SimBase::printSimInfo();

  // TODO
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
      std::cout << mean(idim) << ((idim == (2 - 1))? "": ", ");
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
      std::cout << var(idim) << ((idim == (2 - 1))? "": ", ");
    }
    std::cout << std::endl;
  }
}

// =============================================================================
bool SimNonParam::startLog(const std::string& fname)
{
  // start log from nonparametric
  std::cout << "    t   |";
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "R[" << setw(7 * n_dim - 7) << i << "].x |";
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

  // TODO
  // write in the particles at the end


/*
  // output of gnuplot command
  const int off_next_robot = n_dim + n_dim + n_dim*n_dim + 2;
  std::cout << std::endl;
  std::cout << "~~~ gnuplot command (errors vs determinants) ~~~" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "--- ROBOT " << (i + 1) << " ---" << std::endl;
    std::cout << "clear" << std::endl;
    std::cout << "unset object" << std::endl;
    std::cout << "plot \"output.dat\" u 1:"
              << std::to_string(2+(i+1)*off_next_robot-2)
              << " title \"|Sigma|^0.5 of robot" << std::to_string(1+i) << "\" with line, ";
    std::cout << "\"output.dat\" u 1:"
              << std::to_string(2+(i+1)*off_next_robot-1)
              << " title \"err of robot" << std::to_string(1+i) << "\" with line";
    std::cout << std::endl;
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
              << std::to_string(2+(i+1)*off_next_robot-1)
              << " title \"err" << std::to_string(1+i) << "\" with line";
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
  f_gnuplot << "unset object" << std::endl;
  for (int i = 0; i < n_robots; ++i)
  {
    Eigen::EigenSolver<MatrixXd> s(vars[i]);
    auto eigen_val = s.eigenvalues();
    auto eigen_vec = s.eigenvectors();
    double var_ang
      = std::atan2(eigen_vec.col(0)[1].real(), eigen_vec.col(0)[0].real())
        / M_PI*180.0;
    f_gnuplot << "set object " << std::to_string(i + 1)
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
              << " with linespoints lt -1 lw 1.0 ps 3.0"
              << " pt " << std::to_string(i + 1)
              << " lc palette" << std::endl;
  }
  for (int i = 0; i < n_robots; ++i)
  {
    f_gnuplot << "replot ";
    f_gnuplot << "\"output.dat\" u ";
    for (int j = 0; j < n_dim; ++j)
    {
      f_gnuplot << std::to_string(2+n_dim+i*off_next_robot+j);
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

    "mpl.rcParams['text.usetex'] = True" << std::endl <<

    "cm_new = False" << std::endl <<

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

    "subtle_alpha = 0.25" << std::endl <<
    "def cmf_s(val):" << std::endl <<
    "    if cm_new:" << std::endl <<
    "        val /= maxt" << std::endl <<
    "        c = plt.cm.jet(0.97-val*0.82)" << std::endl <<
    "        c[3] = subtle_alpha" << std::endl <<
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
    "mpl.cm.register_cmap(name='cm_s', cmap=cm)" << std::endl <<

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
    "            (data[:,1+10*i:2+10*i]," << std::endl <<
    "             data[:,2+10*i:3+10*i]," << std::endl <<
    "             data[:,0:1]),axis=1).tolist()," << std::endl <<
    "        markers[i]+'--', cmf_s, lw=1, ms=4)" << std::endl <<
    "# plot estimated trajectories" << std::endl <<
    "for i in range(8):" << std::endl <<
    "    print('plotting estimation of robot ' + str(i+1) + '...')" << std::endl <<
    "    plotPath(" << std::endl <<
    "        np.concatenate(" << std::endl <<
    "            (data[:,3+10*i:4+10*i]," << std::endl <<
    "             data[:,4+10*i:5+10*i]," << std::endl <<
    "             data[:,0:1]),axis=1).tolist()," << std::endl <<
    "        markers[i], cmf, lw=1, ms=4)" << std::endl <<

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
    "cb = mpl.colorbar.ColorbarBase(ax_cb, cmap='cm', norm=norm, orientation='vertical')" << std::endl <<
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
    "    fontsize=5)" << std::endl <<

    "# add covariances" << std::endl <<
    "print('adding covariances...')" << std::endl <<
    "def addEllipse(center, r1, r2, deg):" << std::endl <<
    "    center = np.array(center)" << std::endl <<
    "    major_ax = r1" << std::endl <<
    "    minor_ax = r2" << std::endl <<
    "    angle_deg = deg" << std::endl <<
    "    patch = mpatches.Ellipse(center, major_ax, minor_ax, angle_deg, fc='none', ls='solid', ec='k', lw='1', zorder=100)" << std::endl <<
    "    ax.add_patch(patch)" << std::endl;

  for (int i = 0; i < n_robots; ++i)
  {
    Eigen::EigenSolver<MatrixXd> s(vars[i]);
    auto eigen_val = s.eigenvalues();
    auto eigen_vec = s.eigenvectors();
    double var_ang
      = std::atan2(eigen_vec.col(0)[1].real(), eigen_vec.col(0)[0].real())
        / M_PI*180.0;
    f_pyplot << "addEllipse(["
             << std::to_string(last_mean[i](0)) << ", "
             << std::to_string(last_mean[i](1)) << "], "
             << std::to_string(std::sqrt(eigen_val[0].real()) * 2) << ", "
             << std::to_string(std::sqrt(eigen_val[1].real()) * 2) << ", "
             << std::to_string(var_ang) << ")" << std::endl;
  }

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
  double total_error = 0;
  for (int i = 0; i < n_robots; ++i)
  {
    std::cout << "robot[" << i << "]'s average error:" << (errors[i]/(max_time*sim_freq)) << std::endl;
    total_error += errors[i];
  }
  std::cout << "overall average error: " << (total_error/(max_time*sim_freq)/n_robots) << std::endl;
*/
}

// =============================================================================
void SimNonParam::plotImpl()
{
  // for the terminal,
  // just print the actual locations and errors

  // for the file
  // print the actual location, the average estimated locations, and errors

  std::cout << std::fixed << std::setprecision(2);
  std::cout << std::right << std::setw(8) << t << "|";
  fout << std::fixed << std::setprecision(3);
  fout << std::right << std::setw(8) << t << " ";

  for (int i = 0; i < n_robots; ++i)
  {
    // calculate average estimated location
    VectorXd average = VectorXd::Zero(n_dim);
    for (int ip = 0; ip < n_particles; ++ip)
    {
      average += this->ests[i][ip];
    }
    average = average /= n_particles;

    // current location
    for (int j = 0; j < n_dim; ++j)
    {
      std::cout << std::right << std::setw(6) << robots[i](j);
      if (j < n_dim - 1)
        std::cout << ",";
      else
        std::cout << "|";
    }
    // average estimation
    for (int j = 0; j < n_dim; ++j)
    {
      std::cout << std::right << std::setw(6) << average(j);
      if (j < n_dim - 1)
        std::cout << ",";
      else
        std::cout << "|";
    }

    // error
    std::cout << std::right << std::setw(5) << (robots[i] - average).norm() << "|";

    // current location
    for (int j = 0; j < n_dim; ++j)
    {
      fout << std::right << std::setw(9) << robots[i](j);
      fout << " ";
    }
    // average estimation
    for (int j = 0; j < n_dim; ++j)
    {
      fout << std::right << std::setw(9) << average(j);
      fout << " ";
    }

    // error
    fout << std::right << std::setw(8)
         << (robots[i] - average).norm() << " ";

    // store data to the buffer to plot
    last_loc[i] = robots[i];
    for (int ip = 0; ip < n_particles; ++ip)
    {
      last_est[i][ip] = this->ests[i][ip];
    }
  }

  std::cout << std::endl;
  fout << std::endl;
}

// =============================================================================
void SimNonParam::predict()
{
/*
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
  */
}

// =============================================================================
void SimNonParam::globalLocImpl(const VectorXd& z)
{
  /*
  // update its location estimate
  VectorXd z_hat(2);
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
  VectorXd z_diff = z - z_hat;
  if (z_diff(1) > M_PI)
  {
    z_diff(1) -= 2*M_PI;
  }
  else if (z_diff(1) <= -M_PI)
  {
    z_diff(1) += 2*M_PI;
  }

  globalLocModeImpl(H, Q, z_diff);
  */
}

// =============================================================================
void SimNonParam::calcErrors()
{
  // TODO
  // // calculate errors
  // for (int i = 0; i < n_robots; ++i)
  // {
  //   errors[i] += (means[i] - robots[i]).norm();
  // }
}
