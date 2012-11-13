#ifndef THREE_POINT_FIT_FORMS_H
#define THREE_POINT_FIT_FORMS_H


#include <iostream>
#include <iomanip>
#include <string>
#include "jackknife_fitter.h"
#include "ensem_data.h"


// A1exp(-E1(tf-t)) + A2exp(-E2(t-ti)) + C -- C is the formfactor

struct DoubleExpPlusConst : public FitFunction
{
  DoubleExpPlusConst(const double tf, const double ti)
    : FitFunction(5) , m_tf(tf) , m_ti(ti)
  { 
    setParName(0,"C");
    setParName(1,"A1");
    setParName(2,"E1");
    setParName(3,"A2");
    setParName(4,"E2");
  }

  std::string getFitType(void) const {return std::string("DoubleExpPlusConst");}

  double operator()(const vector<double> &pars, double t) const
  {
    return pars[0] + pars[1]*exp(-pars[2]*(m_tf - t)) + pars[3]*exp(-pars[4]*(t - m_ti));
  }

  double m_tf;
  double m_ti;

};


struct FitComparatorBiasFunction
{
  FitComparatorBiasFunction(double (*fPtr)(const std::vector<double> &))
    : m_fPtr(fPtr)  {  }


  FitComparatorBiasFunction(const FitComparatorBiasFunction &o)
    : m_fPtr(o.m_fPtr) {  }

  virtual double operator()(const std::vector<double> &bias) const
  {
    if(!!!m_fPtr)
    {
      std::cout << __func__ << ": need to give a function pointer" << std::endl;
      exit(1);
    }

    return m_fPtr(bias);
  }

  virtual std::string biasFunctionName(void)  const {return std::string("empty");}

  double (*m_fPtr)(const std::vector<double> &);
};



// an example of a bias function -- exp(-p[0]*(p[1] -p[2])^2)
double gaussianBiasFunction3(const std::vector<double> &p); 

struct FitComparatorGaussianBias : public FitComparatorBiasFunction
{

  FitComparatorGaussianBias(void)
    : FitComparatorBiasFunction(&gaussianBiasFunction3)  { }

  FitComparatorGaussianBias(const FitComparatorGaussianBias &o)
    : FitComparatorBiasFunction(o)  {}


  std::string biasFunctionName(void) const {return std::string("gaussianBiasFunction3");}
};



// that didn't work too well, lets try sin(x)/x to make it fall less rapidly
// -- sin(x)/x , x = (p[0]-p[1])*pi/p[3]
double sinXdivXBiasFunction3(const std::vector<double> &p);


struct FitComparatorSinXdivXBias : public FitComparatorBiasFunction
{
  FitComparatorSinXdivXBias(void)
    : FitComparatorBiasFunction(&sinXdivXBiasFunction3) {}

  FitComparatorSinXdivXBias(const FitComparatorSinXdivXBias &o)
    : FitComparatorBiasFunction(o)  {}

  std::string biasFunctionName(void) const {return std::string("sinXdivXBiasFunction3");}
};

#if 0
// this lets us apply an extra bias by 
struct BiasedFitComparator : public FitComparator
{
  BiasedFitComparator(const std::string &fitCompType, const FitComparatorBiasFunction &biasFunction)
    :  FitComparator(), m_fitCompType(fitCompType) , m_biasFunction(biasFunction)
  { 
    if(m_fitCompType == std::string("Chisq"))
      m_fitComp  = ADAT::Handle<FitComparator>(new CompareFitsByChisqPerNDoF); 

    m_fitComp  = ADAT::Handle<FitComparator>(new CompareFitsByChisqPerNDoF); 
  }

  double operator()(const FitDescriptor &fitDesc, const JackFit& fit) const
  {
    return (m_fitComp->operator()(fitDesc,fit) * m_biasFunction(fit.getBiasParameters()));
  }

  virtual std::string biasFunctionName(void) const {return m_biasFunction->biasFunctionName();}

  std::string m_fitCompType;
  ADAT::Handle<FitComparator> m_fitComp;
  ADAT::Handle<FitComparatorBiasFunction> m_biasFunction; 

};

#endif

struct BiasedFitComparator : public FitComparator
{
  BiasedFitComparator(ADAT::Handle<FitComparator> &fitComp, ADAT::Handle<FitComparatorBiasFunction> &biasFunction)
    :  FitComparator(), m_fitComp(fitComp) , m_biasFunction(biasFunction)
  {   }

  double operator()(const FitDescriptor &fitDesc, const JackFit& fit) const
  {
    return (m_fitComp->operator()(fitDesc,fit) * m_biasFunction->operator()(fit.getBiasParameters()));
  }

  virtual std::string biasFunctionName(void) const {return m_biasFunction->biasFunctionName();}

  ADAT::Handle<FitComparator> m_fitComp;
  ADAT::Handle<FitComparatorBiasFunction> m_biasFunction; 
};



// the three point fitter for radmat
struct FitThreePoint 
{

  // constructor
  FitThreePoint(EnsemData data, int t_f, int t_i, ADAT::Handle<FitComparator> fitComp, int minTSlice);

  // the fit
  void saveFitPlot(const std::string &filename) const;

  // self explanatory 
  std::string getFitPlotString(void) const {return m_axis_plot;}
  std::string getFitSummary(void) const {return m_fit_summary;}
  std::string getFitType(void) const {return m_fit_type;}
  std::string getFitName(void) const {return m_best_fit_name;}

  // FF is the form factor, everything else is just fit parameters
  ENSEM::EnsemReal getFF(void) const {return m_FF;}
  ENSEM::EnsemReal getA1(void) const {return m_A1;}
  ENSEM::EnsemReal getE1(void) const {return m_E1;}
  ENSEM::EnsemReal getA2(void) const {return m_A2;}
  ENSEM::EnsemReal getE2(void) const {return m_E2;}

  // self explanatory 
  double getChisq(void) const {return m_chisq;}
  double getNDoF(void) const {return m_nDoF;}

  // leaving this public
  ADAT::Handle<FitComparator> m_fitComp;
  JackFitLog m_fits;
  std::string m_fit_summary;
  std::string m_axis_plot;
  std::string m_fit_type;
  std::string m_best_fit_name;  
  double m_chisq;
  double m_nDoF;
  ENSEM::EnsemReal m_FF;
  ENSEM::EnsemReal m_A1;
  ENSEM::EnsemReal m_E1; 
  ENSEM::EnsemReal m_A2;
  ENSEM::EnsemReal m_E2;

};

#endif /* THREE_POINT_FIT_FORMS_H */
