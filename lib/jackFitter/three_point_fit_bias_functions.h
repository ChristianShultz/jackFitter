#ifndef THREE_POINT_FIT_BIAS_FUNCTIONS_H
#define THREE_POINT_FIT_BIAS_FUNCTIONS_H 

#include <vector>
#include "jackknife_fitter.h"

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




#endif /* THREE_POINT_FIT_BIAS_FUNCTIONS_H */
