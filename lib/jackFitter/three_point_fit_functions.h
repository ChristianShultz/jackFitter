#ifndef THREE_POINT_FIT_FUNCTIONS_H
#define THREE_POINT_FIT_FUNCTIONS_H 


#include <string>
#include <vector>
#include "jackknife_fitter.h"



////////////////////////////////////////////////////////////////
///////////////////     FUNCTIONS    ///////////////////////////
////////////////////////////////////////////////////////////////


struct ThreePointConstant : public FitFunction
{
  ThreePointConstant(void)
    : FitFunction(1) 
  {
    setParName(0,"C");
  }

  std::string getFitType(void) const {return std::string("ThreePointConstant");}

  double operator()(const vector<double> &pars, double t) const
  {
    return pars[0];
  }

};



// A1exp(-E1(tf-t)) + A2exp(-E2(t-ti)) + C -- C is the formfactor

struct ThreePointDoubleExpPlusConst : public FitFunction
{
  ThreePointDoubleExpPlusConst(const double tf, const double ti)
    : FitFunction(5) , m_tf(tf) , m_ti(ti)
  { 
    setParName(0,"C");
    setParName(1,"Ar");
    setParName(2,"Er");
    setParName(3,"Al");
    setParName(4,"El");

    // CJS -- tuning this for 743 
    // setParamLowerLimit("E1",0.0); 
    // setParamLowerLimit("E2",0.0); 
  }

  std::string getFitType(void) const {return std::string("ThreePointDoubleExpPlusConst");}

  double operator()(const vector<double> &pars, double t) const
  {
    return pars[0] + pars[1]*exp(-pars[2]*(m_tf - t)) + pars[3]*exp(-pars[4]*(t - m_ti));
  }

  double m_tf;
  double m_ti;

};



// A1*( exp (-E1(tf-t) ) + exp(-E2(t-ti))) + C -- C is the formfactor

struct ThreePointSymmetricExpPlusConst : public FitFunction
{
  ThreePointSymmetricExpPlusConst(const double tf, const double ti)
    : FitFunction(3) , m_tf(tf) , m_ti(ti)
  { 
    setParName(0,"C");
    setParName(1,"A");
    setParName(2,"E");

    // CJS -- tuning this for 743 
    // setParamLowerLimit("E",0.0); 
  }

  std::string getFitType(void) const {return std::string("ThreePointSymmetricExpPlusConst");}

  double operator()(const vector<double> &pars, double t) const
  {
    return pars[0] + pars[1]*(exp(-pars[2]*(m_tf - t)) + exp(-pars[2]*(t - m_ti)) );
  }

  double m_tf;
  double m_ti;

};


// A*exp (-E(t-ti) )  + C -- C is the formfactor

struct ThreePointLeftExpPlusConst : public FitFunction
{
  ThreePointLeftExpPlusConst(const double tf, const double ti)
    : FitFunction(3) , m_tf(tf) , m_ti(ti)
  { 
    setParName(0,"C");
    setParName(1,"A");
    setParName(2,"E");
    // setParamLowerLimit("E",0.0); 
  }

  std::string getFitType(void) const {return std::string("ThreePointLeftExpPlusConst");}

  double operator()(const vector<double> &pars, double t) const
  {
    return pars[0] + pars[1]*exp(-pars[2]*(t - m_ti)) ;
  }

  double m_tf;
  double m_ti;

};

// A*exp (-E(tf-t) )  + C -- C is the formfactor

struct ThreePointRightExpPlusConst : public FitFunction
{
  ThreePointRightExpPlusConst(const double tf, const double ti)
    : FitFunction(3) , m_tf(tf) , m_ti(ti)
  { 
    setParName(0,"C");
    setParName(1,"A");
    setParName(2,"E");
    // setParamLowerLimit("E",0.00); 
  }

  std::string getFitType(void) const {return std::string("ThreePointRightExpPlusConst");}

  double operator()(const vector<double> &pars, double t) const
  {
    return pars[0] + pars[1]*exp(-pars[2]*(m_tf - t)) ;
  }

  double m_tf;
  double m_ti;

};

template<int N> 
struct LinearSumExp : public FitFunction
{
  LinearSumExp(void) 
    : FitFunction( 3*N )
  {
    for(int i = 0; i < N; ++i)
    {
      std::stringstream ss; 
      std::string idx; 
      ss << i; 
      idx = ss.str(); 
      int j = 3*i; 
      setParName(j,"E" + idx); 
      setParName(j + 1, "Z" + idx + "l");
      setParName(j + 2, "Z" + idx + "r"); 
    }
  }

  std::string getFitType(void) const {return std::string("LinearSumExp");}

  double operator()(const std::vector<double> &pars, double t) const 
  {
    double sum = 0; 
    for(int i =0; i < N; ++i)
    {
      int j = 3*i; 
      sum += exp(-t * pars[j]) * pars[j+1] *pars[j+2]; 
    }
    return sum; 
  }
};


#endif /* THREE_POINT_FIT_FUNCTIONS_H */
