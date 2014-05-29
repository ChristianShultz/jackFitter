/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : three_point_fit_bias_functions.cc

* Purpose :

* Creation Date : 28-05-2014

* Last Modified : Wed 28 May 2014 09:43:47 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "three_point_fit_bias_functions.h"


// an example of the a bias function
double noBiasFunction(const std::vector<double> &p)
{
  return 1.;
}

// an example of a bias function -- exp(-p[0]*(p[1] -p[2])^2)
double gaussianBiasFunction3(const std::vector<double> &p)
{
  if(p.size() != 3)
  {
    std::cout << __func__ << ": wrong number of parameters" << std::endl;
    exit(1);
  }

  return exp(-p[0]*(p[1] -p[2])*(p[1] -p[2]));
}

// that didn't work too well, lets try sin(x)/x to make it fall less rapidly
// -- sin(x)/x , x = (p[0]-p[1])*pi/p[2]
double sinXdivXBiasFunction3(const std::vector<double> &p)
{
  if(p.size() != 3)
  {
    std::cout << __func__ << ": wrong number of parameters" << std::endl;
    exit(1);
  }

  if(p[0] - p[1] == 0.)
    return 1.;

  double x = (p[0] - p[1])*acos(-1.)/p[2];

  return sin(x)/x;
}


