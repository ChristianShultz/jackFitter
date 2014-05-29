#ifndef THREE_POINT_FIT_COMPATATORS_H
#define THREE_POINT_FIT_COMPATATORS_H 


#include <string>
#include <vector>
#include "jackknife_fitter.h"
#include "ensem_data.h"
#include "adat/handle.h"
#include "three_point_fit_bias_functions.h"
#include "three_point_xml.h"

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


// quick way to make a comparator from xml
ADAT::Handle<FitComparator> constructThreePointFitComparator(const ThreePointComparatorProps_t &);


#endif /* THREE_POINT_FIT_COMPATATORS_H */
