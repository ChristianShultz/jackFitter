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


struct FitThreePoint 
{

  FitThreePoint(EnsemData data, int t_f, int t_i, ADAT::Handle<FitComparator> fitComp, int minTSlice);

  void saveFitPlot(const std::string &filename) const;

  std::string getFitPlotString(void) const {return m_axis_plot;}
  std::string getFitSummary(void) const {return m_fit_summary;}
  std::string getFitType(void) const {return m_fit_type;}
  std::string getFitName(void) const {return m_best_fit_name;}

  ENSEM::EnsemReal getFF(void) const {return m_FF;}
  ENSEM::EnsemReal getA1(void) const {return m_A1;}
  ENSEM::EnsemReal getE1(void) const {return m_E1;}
  ENSEM::EnsemReal getA2(void) const {return m_A2;}
  ENSEM::EnsemReal getE2(void) const {return m_E2;}

  double getChisq(void) const {return m_chisq;}
  double getNDoF(void) const {return m_nDoF;}

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
