/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : three_point_fit_compatators.cc

* Purpose :

* Creation Date : 28-05-2014

* Last Modified : Wed 28 May 2014 09:58:23 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "three_point_fit_compatators.h"


namespace
{

  //////////////////////////////////////////////////////////////////////
  //
  // return a FitComparatorBiasFunction
  ADAT::Handle<FitComparatorBiasFunction> getBiasComparator(const std::string &name)
  {

    typedef std::map<std::string, ADAT::Handle<FitComparatorBiasFunction> > type_t; 
    typedef type_t::value_type value_type; 
    type_t m_map; 

    m_map.insert(value_type( "gaussianBiasFunction3", ADAT::Handle<FitComparatorBiasFunction>(new FitComparatorGaussianBias())));
    m_map.insert(value_type(  "sinXdivXBiasFunction3", ADAT::Handle<FitComparatorBiasFunction>(new FitComparatorSinXdivXBias())));
    m_map.insert(value_type( "noBiasFunction", ADAT::Handle<FitComparatorBiasFunction>(new FitComparatorNoBias())));
    m_map.insert(value_type( "none", ADAT::Handle<FitComparatorBiasFunction>(new FitComparatorNoBias())));

    type_t::const_iterator it; 
    it = m_map.find(name); 
    if(it != m_map.end())
      return it->second; 

    std::cout << __func__ << ": error: " << name << " is unsupported, options were "<< std::endl;
    for(it = m_map.begin(); it != m_map.end(); ++it)
      std::cout << it->first << std::endl;

    std::cout << "exiting.." << std::endl;
    exit(1); 
  }



  //////////////////////////////////////////////////////////////////////
  //
  // return an unmodified fit comparator
  ADAT::Handle<FitComparator> getBasicComparator(const std::string &name)
  {

    typedef std::map<std::string, ADAT::Handle<FitComparator> > type_t; 
    typedef type_t::value_type value_type; 
    type_t m_map; 

    m_map.insert(value_type("CompareFitsByQN", ADAT::Handle<FitComparator>(new ThreePointCompareFitsByQN())));
    m_map.insert(value_type("CompareFitsByChisqPerNDoF", ADAT::Handle<FitComparator>(new ThreePointCompareFitsByChisqPerNDoF())));
    m_map.insert(value_type("CompareFitsByQ", ADAT::Handle<FitComparator>(new ThreePointCompareFitsByQ())));

    type_t::const_iterator it; 
    it = m_map.find(name); 
    if(it != m_map.end())
      return it->second; 

    std::cout << __func__ << ": error: " << name << " is unsupported, options were "<< std::endl;
    for(it = m_map.begin(); it != m_map.end(); ++it)
      std::cout << it->first << std::endl;

    std::cout << "exiting.." << std::endl;
    exit(1); 
  }



  //////////////////////////////////////////////////////////////////////
  //
  // return a biased comparator
  ADAT::Handle<FitComparator> getExtraComparator(const std::string &name)
  {
    typedef std::map<std::string, ADAT::Handle<FitComparator> > type_t; 
    typedef type_t::value_type value_type; 
    type_t m_map; 

    m_map.insert(value_type("LeanParameters",ADAT::Handle<FitComparator>(new ThreePointLeanParameters()))); 
    m_map.insert(value_type("LeanAmplitudes",ADAT::Handle<FitComparator>(new ThreePointLeanAmplitudes()))); 
    m_map.insert(value_type("LeanConstant",ADAT::Handle<FitComparator>(new ThreePointLeanConstant()))); 
    m_map.insert(value_type("BigE",ADAT::Handle<FitComparator>(new ThreePointBigE()))); 
    m_map.insert(value_type("SmallAmplitudes",ADAT::Handle<FitComparator>(new ThreePointSmallAmplitudes()))); 


    type_t::const_iterator it; 
    it = m_map.find(name); 
    if(it != m_map.end())
      return it->second; 

    std::cout << __func__ << ": error: " << name << " is unsupported, options were "<< std::endl;
    for(it = m_map.begin(); it != m_map.end(); ++it)
      std::cout << it->first << std::endl;

    std::cout << "exiting.." << std::endl;
    exit(1); 

  }

  //////////////////////////////////////////////////////////////////////
  //
  // get the vector of comparators
  std::vector<ADAT::Handle<FitComparator> > getExtraComparators(const ThreePointComparatorProps_t &prop)
  {
    std::vector<ADAT::Handle<FitComparator> > ret; 
    for(int i = 0; i < prop.extraProps.size(); ++i)
      ret.push_back(getExtraComparator(prop.extraProps[i]));
    return ret; 
  }

} // namespace anonomyous


//////////////////////////////////////////////////////////////////////
//
// the actual comparator we will actually be using
ADAT::Handle<FitComparator> constructThreePointFitComparator(const ThreePointComparatorProps_t &prop)
{
  ADAT::Handle<FitComparator> multiComp(new
      ThreePointMultiComparator(getBasicComparator(prop.baseProp),getExtraComparators(prop)));
  ADAT::Handle<FitComparatorBiasFunction> biasFunction(getBiasComparator(prop.biasProp));

  return ADAT::Handle<FitComparator>(new BiasedFitComparator(multiComp,biasFunction));
}

