#ifndef THREE_POINT_FIT_FORMS_H
#define THREE_POINT_FIT_FORMS_H


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <numeric> 
#include "jackknife_fitter.h"
#include "ensem_data.h"
#include "adat/handle.h"
#include "io/adat_xmlio.h"




////////////////////////////////////////////////////////////////
///////////////////     XML STUFF    ///////////////////////////
////////////////////////////////////////////////////////////////

struct ThreePointComparatorProps_t 
{
  std::string fit_type;
  std::string baseProp; 
  std::string biasProp;
  ADATXML::Array<std::string> extraProps;
  int tlow;
  int thigh; // the fit range
  int minTSlice; 
};

std::string toString(const ThreePointComparatorProps_t &);
std::ostream& operator<<(std::ostream&, const ThreePointComparatorProps_t &);
void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointComparatorProps_t &); 


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
    setParName(1,"A1");
    setParName(2,"E1");
    setParName(3,"A2");
    setParName(4,"E2");
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
  ThreePointDoubleExpPlusConst(const double tf, const double ti)
    : FitFunction(5) , m_tf(tf) , m_ti(ti)
  { 
    setParName(0,"C");
    setParName(1,"A");
    setParName(2,"E");
  }

  std::string getFitType(void) const {return std::string("ThreePointSymmetricExpPlusConst");}

  double operator()(const vector<double> &pars, double t) const
  {
    return pars[0] + pars[1]*(exp(-pars[2]*(m_tf - t)) + exp(-pars[2]*(t - m_ti)) );
  }

  double m_tf;
  double m_ti;

};




////////////////////////////////////////////////////////////////
//////////////////     COMPARATORS   ///////////////////////////
////////////////////////////////////////////////////////////////

// allow for easy stacking on of extra fit comparators
struct ThreePointMultiComparator : public FitComparator
{
  ThreePointMultiComparator(const ADAT::Handle<FitComparator> &base_comparator, const std::vector<ADAT::Handle<FitComparator> > &pieces) 
    : FitComparator() , m_base_comparator(base_comparator) , m_comps(pieces) 
  {  }

  inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const
  {

    double ret = (*m_base_comparator)(fitDesc,fit);

    if(!!!m_comps.empty())
    {
      std::vector<ADAT::Handle<FitComparator> >::const_iterator it;
      for(it = m_comps.begin(); it != m_comps.end(); ++it)
        ret *= (*it)->operator()(fitDesc,fit);
    }
    return ret; 
  }

  ADAT::Handle<FitComparator> m_base_comparator; 
  std::vector<ADAT::Handle<FitComparator> > m_comps; 
};


// base comparators
// should use one of the following -- most are stolen from old code, they were scattered so collect here under new names 

// maximize NDoF *statQ
class ThreePointCompareFitsByQN : public FitComparator
{
  public:
    ThreePointCompareFitsByQN() : FitComparator(){};

    inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const
    {
      return fit.getNDoF()* statQ( fit.getAvgChisq(), fit.getNDoF() );
    }
};

//this is an example where the inverse of the chisq / nDoF is returned (to be maximised)
class ThreePointCompareFitsByChisqPerNDoF : public FitComparator
{
  public:
    ThreePointCompareFitsByChisqPerNDoF() : FitComparator(){};

    inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const
    {
      return double(fit.getNDoF()) / fit.getAvgChisq();

    }
};


//this is an example where the Q of the fit is returned (to be maximised)
class ThreePointCompareFitsByQ : public FitComparator
{
  public:
    ThreePointCompareFitsByQ() : FitComparator(){};

    inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const
    {
      return statQ( fit.getAvgChisq(), fit.getNDoF() );

    }
};



// can use as many of thes as like

// move toward solutions with low noise -- still have flat directions
// where the noise in one parameter goes to zero and the noise in another to +inf
struct ThreePointLeanParameters : public FitComparator
{
  ThreePointLeanParameters() : FitComparator() {}

  inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const
  {
    vector<double> errs = fit.getAvgFitParErrors();
    double dum(1.);
    std::vector<double>::const_iterator start,stop;
    start = errs.begin(); 
    stop = errs.end();
    while(start != stop)
      dum /= (*start++);
    return dum;
  }
};


// a stab at biasing towards lean amplitudes
struct ThreePointLeanAmplitudes : public FitComparator
{
  ThreePointLeanAmplitudes() : FitComparator() {}

  inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const
  {

    vector<double> errs = fit.getAvgFitParErrors();

    double dum(1.); 

    if( (*(fitDesc.ff)).getFitType() == "ThreePointDoubleExpPlusConst" )
    {
      dum  = 1./(errs[1]*errs[3]);
    }

    return dum;
  }

};


// should also add in a lean constant
struct ThreePointLeanConstant: public FitComparator
{
  ThreePointLeanConstant() : FitComparator() {}

  inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const
  {

    vector<double> errs = fit.getAvgFitParErrors();

    return 1./errs[0];
  }
};

// a stab at biasing towards higher masses
struct ThreePointBigE : public FitComparator
{
  ThreePointBigE() : FitComparator() {}

  inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const
  {

    vector<double> pars = fit.getAvgFitParValues();
    vector<double> errs = fit.getAvgFitParErrors();

    double dum(1.); 

    if( (*(fitDesc.ff)).getFitType() == "ThreePointDoubleExpPlusConst" )
    {
      dum  *= fabs(pars[2]/errs[2]) * fabs(pars[4]/errs[4]);
    }

    return dum;
  }

};


struct ThreePointSmallAmplitudes : public FitComparator
{
  ThreePointSmallAmplitudes() : FitComparator() {}

  inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const
  {
    vector<double> pars = fit.getAvgFitParValues();
    vector<double> errs = fit.getAvgFitParErrors();

    double dum(1.);

    if( (*(fitDesc.ff)).getFitType() == "DoubleExpPlusConst" )
    {
      dum  /= (pars[1]*pars[3]);
    }

    return dum;
  }

};



////////////////////////////////////////////////////////////////
/////////////////   BIAS  FUNCTIONS    /////////////////////////
////////////////////////////////////////////////////////////////


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


// no bi
double noBiasFunction(const std::vector<double> &p);

struct FitComparatorNoBias : public FitComparatorBiasFunction
{
  FitComparatorNoBias(void)
    : FitComparatorBiasFunction(&noBiasFunction) {}

  FitComparatorNoBias(const FitComparatorNoBias &o)
    : FitComparatorBiasFunction(o) {}

  std::string biasFunctionName(void) const {return std::string("noBiasFunction");}

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



// quick way to make a comparator from xml
ADAT::Handle<FitComparator> constructThreePointFitComparator(const ThreePointComparatorProps_t &);



////////////////////////////////////////////////////////////////
/////////////////     THREE-PT FITTER    ///////////////////////
////////////////////////////////////////////////////////////////




// the three point fitter for radmat
struct FitThreePoint 
{

  // constructor
  FitThreePoint(EnsemData data, int t_f, int t_i, ADAT::Handle<FitComparator> fitComp, int minTSlice, const std::string &fit_type);

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
  std::string m_axis_plot_component;
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
