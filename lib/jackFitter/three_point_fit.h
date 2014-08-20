#ifndef THREE_POINT_FIT_H
#define THREE_POINT_FIT_H 

#include "three_point_fit_compatators.h"
#include "three_point_fit_functions.h"
#include "three_point_xml.h"

////////////////////////////////////////////////////////////////
/////////////////     THREE-PT FITTER    ///////////////////////
////////////////////////////////////////////////////////////////




// the three point fitter for radmat
struct FitThreePoint 
{

  // constructor
  FitThreePoint(EnsemData data, 
      const int tsnk,                        // sink position 
      const int tsrc,                        // source position
      const int t_f,                         // high end of fit range
      const int t_i,                         // low end of fit range 
      ADAT::Handle<FitComparator> fitComp,
      const int minTSlice,                   // min length of fit range
      const std::string &fit_type, 
      const bool useVals = false,
      const FitParValue=FitParValue());



  // the fit
  void saveFitPlot(const std::string &filename) const;

  // self explanatory 
  std::string getFitPlotString(void) const {return m_axis_plot;}
  std::string getFitSummary(void) const {return m_fit_summary;}
  std::string getFitType(void) const {return m_fit_type;}
  std::string getFitName(void) const {return m_best_fit_name;}

  // basically write out the old axis plot but also plot the components that 
  // went into the fit
  std::string getFitPlotStringWithComponents(void) const {return m_axis_plot_component;}

  // FF is the form factor, everything else is just a fit parameter
  ENSEM::EnsemReal getFF(void) const {return m_FF;}
  ENSEM::EnsemReal getAl(void) const {return m_Al;}
  ENSEM::EnsemReal getEl(void) const {return m_El;}
  ENSEM::EnsemReal getAr(void) const {return m_Ar;}
  ENSEM::EnsemReal getEr(void) const {return m_Er;}

  // self explanatory 
  double getChisq(void) const {return m_chisq;}
  double getNDoF(void) const {return m_nDoF;}

  int tlow(void) {return t_low;}
  int thigh(void) {return t_high;}

  int t_low, t_high;
  
  // leaving this public
  ADAT::Handle<FitComparator> m_fitComp;
  JackFitLog m_fits;
  std::string m_fit_summary;
  std::string m_axis_plot;
  std::string m_axis_plot_component;
  std::string m_fit_type;
  std::string m_best_fit_name;  
  double m_chisq;
  double m_nDoF;
  double dlow,dhigh; 
  ENSEM::EnsemReal m_FF;
  ENSEM::EnsemReal m_Al;
  ENSEM::EnsemReal m_El; 
  ENSEM::EnsemReal m_Ar;
  ENSEM::EnsemReal m_Er;

};

#endif /* THREE_POINT_FIT_H */
