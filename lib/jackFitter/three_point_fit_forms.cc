/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : three_point_fit_forms.cc

 * Purpose :

 * Creation Date : 09-11-2012

 * Last Modified : Mon 14 Oct 2013 06:37:42 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "three_point_fit_forms.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <algorithm> 



// xml junk


std::string toString(const ThreePointComparatorProps_t &prop)
{
  std::stringstream ss;
  ss << "fit_type = " << prop.fit_type;  
  ss << " baseProp = " << prop.baseProp << "   biasProp = " << prop.biasProp << "   extraProps = ";
  for(int i = 0; i < prop.extraProps.size(); ++i)
    ss << prop.extraProps[i] << " ";
  ss << "   tlow = " << prop.tlow << "   thigh = " << prop.thigh << "   minTSlice = " << prop.minTSlice; 
  return ss.str();
}

std::ostream& operator<<(std::ostream &o, const ThreePointComparatorProps_t &prop)
{
  o << toString(prop);
  return o;
}

namespace
{

  template<typename T>
    void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
    {
      if(ptop.count(path) > 0)
        read(ptop,path,place);
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " << f << " trying to read path, " << path
          << ", path was empty, exiting" << std::endl;
        exit(1);
      }
    }

} // namespace anonomyous 



void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointComparatorProps_t &prop)
{
  ADATXML::XMLReader ptop(xml,path);
  doXMLRead(ptop,"fit_type",prop.fit_type,__PRETTY_FUNCTION__); 
  doXMLRead(ptop,"baseProp",prop.baseProp,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"biasProp",prop.biasProp,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"extraProps",prop.extraProps,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"tlow",prop.tlow,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"thigh",prop.thigh,__PRETTY_FUNCTION__); 
  doXMLRead(ptop,"minTSlice",prop.minTSlice,__PRETTY_FUNCTION__);
}

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



namespace
{
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

  std::vector<ADAT::Handle<FitComparator> > getExtraComparators(const ThreePointComparatorProps_t &prop)
  {
    std::vector<ADAT::Handle<FitComparator> > ret; 
    for(int i = 0; i < prop.extraProps.size(); ++i)
      ret.push_back(getExtraComparator(prop.extraProps[i]));
    return ret; 
  }

} // namespace anonomyous


ADAT::Handle<FitComparator> constructThreePointFitComparator(const ThreePointComparatorProps_t &prop)
{
  ADAT::Handle<FitComparator> multiComp(new
      ThreePointMultiComparator(getBasicComparator(prop.baseProp),getExtraComparators(prop)));
  ADAT::Handle<FitComparatorBiasFunction> biasFunction(getBiasComparator(prop.biasProp));

  return ADAT::Handle<FitComparator>(new BiasedFitComparator(multiComp,biasFunction));
}





////////////////////////////////////////////////////////////////
/////////////////     THREE-PT FITTER    ///////////////////////
////////////////////////////////////////////////////////////////



namespace
{

  std::vector<double> constructBiasParameters(const std::string &fname, 
      double tmax, 
      double tmin,
      double thigh, 
      double tlow)
  {

    std::vector<double> biasParameters; 

    if(fname == std::string("gaussianBiasFunction3"))
    {
      // try to bias the fits toward the center
      double half_range = double(tmax - tmin)/2.;
      double sigma = 2.*log(20)/(half_range*half_range);             // 1/20 supression on the edges varrying smoothly to centered      
      double midpt = double(tmax - tmin)/2. + tmin;                  // it looked reasonable in mathematica..
      double fit_midpt = double(thigh - tlow)/2. + tlow; 
      biasParameters.push_back(sigma);
      biasParameters.push_back(midpt);
      biasParameters.push_back(fit_midpt);
    }
    else if(fname == std::string("sinXdivXBiasFunction3"))
    {
      double half = double(tmax - tmin)/2.;
      double midpt = double(tmax - tmin)/2. + tmin; 
      double fit_midpt = double(thigh - tlow)/2. + tlow; 
      biasParameters.push_back(midpt);
      biasParameters.push_back(fit_midpt);
      biasParameters.push_back(half);
      
    }
    else if((fname == std::string("none")) || (fname == std::string("noBiasFunction")))
    {
      // nothing     
    }
    else
    {
      std::cout << __func__ << ": fname (" << fname << ") is not supported" << std::endl;
      exit(1);
    }

    return biasParameters;
  }



  std::string makeFancyPlot(const AxisPlot &pt,
      const FitThreePoint &tp,
      const std::string &mode,
      const int tf,
      const int ti)
  {

    if((mode != std::string("dExpPC") && (mode != std::string("Constant") && (mode != std::string("symExpPC")) )))
    {
      std::cout << __func__ << ": error: mode " << mode << "not supported" << std::endl;
      __builtin_trap(); 
    }

    // A1exp(-E1(tf-t)) + A2exp(-E2(t-ti)) + C -- C is the formfactor
    AxisPlot plot(pt); 
    ENSEM::EnsemReal C, A1,A2,E1,E2,dum; 
    C = tp.getFF();
    A1 = tp.getA1();
    A2 = tp.getA2();
    E1 = tp.getE1();
    E2 = tp.getE2(); 


    double dC,dA1,dA2,dE1,dE2;
    dC = ENSEM::toDouble(ENSEM::mean(C));
    dA1 = ENSEM::toDouble(ENSEM::mean(A1)); 
    dA2 = ENSEM::toDouble(ENSEM::mean(A2));
    dE1 = ENSEM::toDouble(ENSEM::mean(E1));
    dE2 = ENSEM::toDouble(ENSEM::mean(E2)); 

    std::vector<double> m,mp,mm, time; 
    double min, max, mean,var, tmp;

    // add the form factor
    mean = ENSEM::toDouble(ENSEM::mean(C)); 
    var = ENSEM::toDouble(ENSEM::variance(C)); 
    min = mean - var; 
    max = mean + var;
    for(int t = 0; t < tf - ti; ++t)
    {
      time.push_back(t+ti); 
      m.push_back(mean);
      mp.push_back(mean + var);
      mm.push_back(mean - var); 
    }

    plot.addLineData(time,m,3);
    plot.addLineData(time,mm,3);
    plot.addLineData(time,mp,3);


    m.clear();
    mm.clear();
    mp.clear();


    // add the sink exp
    for(int t = ti; t <= tf; ++t)
    {
      dum = A1*ENSEM::exp(-E1*ENSEM::Real(double(tf-t)));
      mean = ENSEM::toDouble(ENSEM::mean(dum));
      var = ENSEM::toDouble(ENSEM::variance(dum));
      m.push_back(mean);
      mm.push_back(mean - var);
      mp.push_back(mean + var);
    }

    plot.addLineData(time,m,4);
    plot.addLineData(time,mm,4);
    plot.addLineData(time,mp,4);


    // these things actually sum so keep track of absolute ranges
    std::transform(mp.begin(),mp.end(),mp.begin(),std::bind1st(std::plus<double>(),max)); 

    tmp = *std::max_element(mp.begin(),mp.end());
    max = std::max(max,tmp);
    tmp = *std::min_element(mm.begin(),mm.end());
    min = std::min(min,tmp);

    m.clear();
    mm.clear();
    mp.clear();

    for(int t = ti; t <= tf; ++t)
    {
      dum = A2*ENSEM::exp(-E2*ENSEM::Real(t-ti));
      mean = ENSEM::toDouble(ENSEM::mean(dum));
      var = ENSEM::toDouble(ENSEM::variance(dum));
      m.push_back(mean);
      mm.push_back(mean -var);
      mp.push_back(mean + var);
    } 

    plot.addLineData(time,m,5);
    plot.addLineData(time,mm,5);
    plot.addLineData(time,mp,5);

    // these things actually sum so keep track of absolute ranges
    std::transform(mp.begin(),mp.end(),mp.begin(),std::bind1st(std::plus<double>(),max)); 


    tmp = *std::max_element(mp.begin(),mp.end());
    max = std::max(max,tmp);
    tmp = *std::min_element(mm.begin(),mm.end());
    min = std::min(min,tmp);


    plot.setYRange(min,max);
    plot.setXRange(ti - 2, tf + 7);


    double height = max - min;
    double unit = height/5.;

    std::stringstream lab;
    lab << "FF = " << std::setprecision(4) << dC;
    plot.addLabel(tf+1,0.5*unit + min,lab.str(),3,0.8);
    lab.str("");
    lab << "A1 = " << std::setprecision(4) << dA1;
    plot.addLabel(tf+1,1.5*unit + min,lab.str(),4,0.8);
    lab.str("");
    lab << "E1 = " << std::setprecision(4) << dE1;
    plot.addLabel(tf+1,2.5*unit + min,lab.str(),4,0.8);
    lab.str("");
    lab << "A2 = " << std::setprecision(4) << dA2;
    plot.addLabel(tf+1,3.5*unit + min,lab.str(),5,0.8);
    lab.str("");
    lab << "E2 = " << std::setprecision(4) << dE2;
    plot.addLabel(tf+1,4.5*unit + min,lab.str(),5,0.8);




    return plot.getAxisPlotString(); 

  }



  // do a fit using double exp plus const fit function
  ADAT::Handle<FitFunction> trydExpPCFit(const int t_f,
      const int t_i, 
      const int minTSlice,
      EnsemData data,
      ADAT::Handle<FitComparator> &m_fitComp,
      JackFitLog &m_fits)
  {


    if(minTSlice < 5) 
    {
      std::cout << __func__ << " need at least 5 data points to do a fit" << std::endl;
      exit(1);
    }


    // try the double exp plus const fit 
    ADAT::Handle<FitFunction> dExpPC (new ThreePointDoubleExpPlusConst(t_f,t_i));

    dExpPC->setDefaultParValue("C",1.);     // a terrible guess at the form factor 
    dExpPC->setDefaultParValue("E1",0.18); // this is like delta m so around a pion mass in lattice units?
    dExpPC->setDefaultParValue("E2",0.18); // this is like delta m so around a pion mass in lattice units?

    // guess at the sign by computing the slope -- also avoid contact terms
    ENSEM::Real xp,xm;
    double tp , tm; 
    tp = double(t_f -1);
    tm = double(t_f - 3); 
    xp = ENSEM::mean(data.getYUsingNearestX(tp));
    xm = ENSEM::mean(data.getYUsingNearestX(tm));
    double sign,slope;
    slope = ENSEM::toDouble(xp-xm);
    if(slope> 0)
      sign = 1.;
    else
      sign = -1.;

    dExpPC->setDefaultParValue("A1",sign*1.);

    // guess at the sign by computing the slope -- also avoid contact terms
    tp = double(t_i + 3);
    tm = double(t_i + 1);
    xp = ENSEM::mean(data.getYUsingNearestX(tp));
    xm = ENSEM::mean(data.getYUsingNearestX(tm));
    slope = ENSEM::toDouble(xp-xm);
    if(slope> 0)
      sign = -1.;
    else
      sign = 1.;


    // so imagine not very noisy correlators -- we could try
    // to set lower/upper limits using the slope, not doing it now 
    // but might be a good idea to include the option for 
    // optimized three-point functions


    dExpPC->setDefaultParValue("A2",sign*1.);


    // do a loop over all possible fit separations.. this is probably going to be slow
    for(int t_low = t_i + 1; t_low < t_f; ++t_low)
      for(int t_high = t_f -1; t_high > t_low; --t_high)
      {
        // skip the impossible ones
        if(t_high - t_low < minTSlice)
          continue;

        m_fits.getEnsemData().showAll();
        m_fits.getEnsemData().hideDataAboveX(t_high + 0.1);
        m_fits.getEnsemData().hideDataBelowX(t_low - 0.1); 

        // this class is so stupidly broken -- why is it this hard to add something
        std::vector<std::pair<std::string,int> > named_ints; 
        named_ints.push_back(std::pair<std::string,int>("t_low",t_low));
        named_ints.push_back(std::pair<std::string,int>("t_high",t_high));

        std::stringstream ss; 
        ss << "DoubleExpPlusConst: t_low = " << t_low << " t_high = " << t_high; 
        // sanity 2 

        if(m_fits.getEnsemData().getNData() >= minTSlice)
          m_fits.addFit(ss.str(),dExpPC,constructBiasParameters(m_fitComp->biasFunctionName(),t_f, t_i, t_high,t_low),named_ints);   
      }

    return dExpPC; 
  }

  // do a fit using symm double exp plus const fit function
  ADAT::Handle<FitFunction> trySymExpPCFit(const int t_f,
      const int t_i, 
      const int minTSlice,
      EnsemData data,
      ADAT::Handle<FitComparator> &m_fitComp,
      JackFitLog &m_fits)
  {


    if(minTSlice < 3) 
    {
      std::cout << __func__ << " need at least 5 data points to do a fit" << std::endl;
      exit(1);
    }


    // try the double exp plus const fit 
    ADAT::Handle<FitFunction> dExpPC (new ThreePointSymmetricExpPlusConst(t_f,t_i));

    dExpPC->setDefaultParValue("C",1.);     // a terrible guess at the form factor 
    dExpPC->setDefaultParValue("E",0.18); // this is like delta m so around a pion mass in lattice units?

    // guess at the sign by computing the slope -- also avoid contact terms
    ENSEM::Real xp,xm;
    double tp , tm; 
    tp = double(t_f -1);
    tm = double(t_f - 3); 
    xp = ENSEM::mean(data.getYUsingNearestX(tp));
    xm = ENSEM::mean(data.getYUsingNearestX(tm));
    double sign,slope;
    slope = ENSEM::toDouble(xp-xm);
    if(slope> 0)
      sign = 1.;
    else
      sign = -1.;

    dExpPC->setDefaultParValue("A",sign*1.);



    // do a loop over all possible fit separations.. this is probably going to be slow
    for(int t_low = t_i + 1; t_low < t_f; ++t_low)
      for(int t_high = t_f -1; t_high > t_low; --t_high)
      {
        // skip the impossible ones
        if(t_high - t_low < minTSlice)
          continue;

        m_fits.getEnsemData().showAll();
        m_fits.getEnsemData().hideDataAboveX(t_high + 0.1);
        m_fits.getEnsemData().hideDataBelowX(t_low - 0.1); 

        std::stringstream ss; 
        ss << "DoubleExpPlusConst: t_low = " << t_low << " t_high = " << t_high; 
        // sanity 2 
        
        // this class is so stupidly broken -- why is it this hard to add something
        std::vector<std::pair<std::string,int> > named_ints; 
        named_ints.push_back(std::pair<std::string,int>("t_low",t_low));
        named_ints.push_back(std::pair<std::string,int>("t_high",t_high));

        if(m_fits.getEnsemData().getNData() >= minTSlice)
          m_fits.addFit(ss.str(),dExpPC,constructBiasParameters(m_fitComp->biasFunctionName(),t_f, t_i, t_high,t_low),named_ints);   
      }

    return dExpPC; 
  }
  struct stupidDerivative
  {
    stupidDerivative(const double xx, const double ff)
      : x(xx) , fprime(ff)
    {}

    double x;
    double fprime; 
  };

  bool isCompatibleWithConstant(const std::vector<double> x, ENSEM::EnsemVectorReal y)
  {
    const int sz = x.size(); 
    std::vector<stupidDerivative> derivative;
    std::vector<double> yy(x.size(),0.); 
    double mean(0.);

    for(int i =0; i < sz; ++i)
      yy[i] = ENSEM::toDouble(ENSEM::mean(ENSEM::peekObs(y,i))); 

    for(int i = 1; i < sz - 1; ++i)
    {
      derivative.push_back(stupidDerivative(x[i],(yy[i+1] - yy[i-1])/(x[i+1] - x[i-1]))); 
      mean += fabs(yy[i]); 
    }

    std::vector<stupidDerivative>::const_iterator it; 
    double mean_deriv(0.);
    double abs_mean_deriv(0.);

    for(it = derivative.begin(); it != derivative.end(); ++it)
    {
      mean_deriv += it->fprime; 
      abs_mean_deriv += fabs(it->fprime); 
    }

    // try to divide out the scale
    mean_deriv /= (mean*double(sz - 2)); 
    abs_mean_deriv /= (mean*double(sz-2)); 

    if((fabs(mean_deriv) < 0.1) && (abs_mean_deriv < 0.2))
      return true; 

    return false; 
  }


  // do a fit using const fit function
  ADAT::Handle<FitFunction> tryConstantFit(const int t_f,
      const int t_i, 
      const int minTSlice,
      EnsemData data,
      ADAT::Handle<FitComparator> &m_fitComp,
      JackFitLog &m_fits)
  {
    ADAT::Handle<FitFunction> tpConst  (new ThreePointConstant());

    /*
       if(!!!isCompatibleWithConstant(data.getAllXData(), data.getAllYData()))
       return tpConst; 
     */

    std::vector<double> x = data.getAllXData();

    // can't be const as EnsemData doesn't want a const x...
    // still this is a one off vector so it shouldn't be harmful
    std::vector<double>::iterator it;
    double guess(0.);
    for(it = x.begin(); it != x.end(); ++it)
        guess +=  ENSEM::toDouble(ENSEM::mean(data.getYUsingNearestX(*it))); 
    
    guess /= double(x.size()); 

    // global average
    tpConst->setDefaultParValue("C",guess); 

    // do a loop over all possible fit separations.. this is probably going to be slow
    for(int t_low = t_i + 1; t_low < t_f; ++t_low)
      for(int t_high = t_f -1; t_high > t_low; --t_high)
      {
        // skip the impossible ones
        if(t_high - t_low < minTSlice)
          continue;

        m_fits.getEnsemData().showAll();
        m_fits.getEnsemData().hideDataAboveX(t_high + 0.1);
        m_fits.getEnsemData().hideDataBelowX(t_low - 0.1); 

        std::stringstream ss; 
        ss << "Constant: t_low = " << t_low << " t_high = " << t_high; 
        // sanity 2 

        // this class is so stupidly broken -- why is it this hard to add something
        std::vector<std::pair<std::string,int> > named_ints; 
        named_ints.push_back(std::pair<std::string,int>("t_low",t_low));
        named_ints.push_back(std::pair<std::string,int>("t_high",t_high));

        if(m_fits.getEnsemData().getNData() >= minTSlice)
          m_fits.addFit(ss.str(),tpConst,constructBiasParameters(m_fitComp->biasFunctionName(),t_f, t_i, t_high,t_low),named_ints);   
      }

    return tpConst; 
  }


} // namespace anonomyous 






  FitThreePoint::FitThreePoint(EnsemData data, int t_f, int t_i, ADAT::Handle<FitComparator> fitComp, int minTSlice, const std::string &fit_type) 
: m_fits(data) , m_fitComp(fitComp)
{
  ADAT::Handle<FitFunction> dExpPC,Constant,symExpPC; 

  if((fit_type == "symExpPC") || (fit_type == "all"))
    symExpPC = trySymExpPCFit(t_f, t_i,minTSlice, data, m_fitComp, m_fits);
  if((fit_type == "dExpPC") || (fit_type == "all"))
    dExpPC = trydExpPCFit(t_f, t_i,minTSlice, data, m_fitComp, m_fits);
  else if((fit_type == "const") || (fit_type == "all"))
    Constant = tryConstantFit(t_f, t_i,minTSlice, data, m_fitComp, m_fits);
  else
  {
    std::cerr << __func__ << ": unknown fit type " << fit_type << "options are {all , dExpPC , const }" << std::endl;
    exit(1);
  }

  //  FitDescriptor best_dExpPC = m_fits.getBestFit(*m_fitComp);

  int rank(0); // the fits are ordered 
  FitDescriptor best = m_fits.getBestJackFit(*m_fitComp,rank);

  // zero out everything
  double tip1 = double(t_i + 1);  // this reference to double crap is getting annoying..
  m_FF = data.getYUsingNearestX(tip1);
  m_FF = ENSEM::Real(0.);
  m_A1 = m_FF;
  m_E1 = m_FF;
  m_A2 = m_FF;
  m_E2 = m_FF; 


  if(best.fitname == "FAILED")
  {
    m_chisq = 1.e10;
    m_nDoF = 1;
    m_best_fit_name = "FAILED";
    m_fit_summary = "FAILED";
    t_low = t_i; 
    t_high = t_f; 
    // no plot
  }
  else  if((fit_type == "const") || (fit_type == "all"))
  {
    if (best.ff->getFitType() == Constant->getFitType())
    {
      m_fit_type = Constant->getFitType(); 
      JackFit& bestFit = m_fits.getFit(best);
      m_FF = bestFit.getJackFitParValue("C");
      m_chisq = bestFit.getJackChisq();
      m_nDoF = bestFit.getNDoF();
      m_best_fit_name = best.fitname;

      // this gets more frustrating every time i have to look at it
      {
        std::vector<std::pair<std::string,int> > named_ints = bestFit.get_named_ints(); 
        std::vector<std::pair<std::string,int> >::const_iterator it; 
        bool found = false; 
        for(it = named_ints.begin(); it != named_ints.end(); ++it)
          if(it->first == "t_low")
          {
            found = true; 
            t_low = it->second; 
          }
        if (! found) 
          std::cerr << "Warning: something wacky is going on here" << std::endl;
        found = false;  

        for(it = named_ints.begin(); it != named_ints.end(); ++it)
          if(it->first == "t_high")
          {
            found = true; 
            t_high = it->second; 
          }
        if (! found) 
          std::cerr << "Warning: something wacky is going on here" << std::endl;
      }


      int count = 1; 
      std::map<double,FitDescriptor> list = m_fits.getFitList(*m_fitComp);
      std::stringstream ss; 
      ss << "                                           | chisq/nDoF |     Q      |  fitCrit   | " << endl;

      for( std::map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++)
      {
        JackFit& thisFit = m_fits.getFit(p->second);
        double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
        double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
        ss << setw(35) <<(p->second).fitname << "|";
        ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
        ss << setw(12) << fixed << setprecision(3) << Q <<"|";
        ss << setw(12) << scientific << setprecision(3) << p->first << "|";

        if(count == rank)
        {
          ss << "*";
        }
        else
        {
          ss << " ";
        }
        ss << " FF=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("C") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("C");
        ss << endl;
        count++;
      }

      m_fit_summary = ss.str();

      std::stringstream lab;
      lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
      lab << "; FF=" << fixed << setprecision(4) << toDouble(mean(m_FF)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(m_FF)));
      std::stringstream plot;
      plot << bestFit.makeJackFitPlotAxis(t_i - 2, t_f + 2, lab.str());

      m_axis_plot = plot.str();

      AxisPlot fancy_plot = bestFit.getJackFitPlotAxis(t_i -2 , t_f + 2, lab.str());

      m_axis_plot_component =  makeFancyPlot(fancy_plot,*this,std::string("Constant"),t_f, t_i);

    }
  }
  else if((fit_type == "dExpPC") || (fit_type == "all"))
  {
    if (best.ff->getFitType() == dExpPC->getFitType())
    {
      m_fit_type = dExpPC->getFitType(); 
      JackFit& bestFit = m_fits.getFit(best);
      m_FF = bestFit.getJackFitParValue("C");
      m_A1 = bestFit.getJackFitParValue("A1");
      m_E1 = bestFit.getJackFitParValue("E1");
      m_A2 = bestFit.getJackFitParValue("A2");
      m_E2 = bestFit.getJackFitParValue("E2");

      m_chisq = bestFit.getJackChisq();
      m_nDoF = bestFit.getNDoF();
      m_best_fit_name = best.fitname;

      int count = 1; 
      std::map<double,FitDescriptor> list = m_fits.getFitList(*m_fitComp);
      std::stringstream ss; 
      ss << "                                           | chisq/nDoF |     Q      |  fitCrit   | " << endl;

      for( std::map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++)
      {
        JackFit& thisFit = m_fits.getFit(p->second);
        double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
        double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
        ss << setw(35) <<(p->second).fitname << "|";
        ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
        ss << setw(12) << fixed << setprecision(3) << Q <<"|";
        ss << setw(12) << scientific << setprecision(3) << p->first << "|";

        if(count == rank)
        {
          ss << "*";
        }
        else
        {
          ss << " ";
        }
        ss << " FF=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("C") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("C");

        if( ((p->second).ff).operator->() == dExpPC.operator->() )
        {
          ss << ", E1=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("E1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("E1");
          ss << ", E2=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("E2") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("E2");
          ss << ", A1=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("A1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("A1");
          ss << ", A2=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("A2") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("A2");
        }
        ss << endl;count++;
      } 


      m_fit_summary = ss.str();

      std::stringstream lab;
      lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
      lab << "; FF=" << fixed << setprecision(4) << toDouble(mean(m_FF)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(m_FF)));
      std::stringstream plot;
      plot << bestFit.makeJackFitPlotAxis(t_i - 2, t_f + 2, lab.str());

      m_axis_plot = plot.str();

      AxisPlot fancy_plot = bestFit.getJackFitPlotAxis(t_i -2 , t_f + 2, lab.str());

      m_axis_plot_component =  makeFancyPlot(fancy_plot,*this,std::string("dExpPC"),t_f, t_i);

    }
  }
  else if((fit_type == "symExpPC") || (fit_type == "all"))
  {
    if (best.ff->getFitType() == symExpPC->getFitType())
    {
      m_fit_type = dExpPC->getFitType(); 
      JackFit& bestFit = m_fits.getFit(best);
      m_FF = bestFit.getJackFitParValue("C");
      m_A1 = bestFit.getJackFitParValue("A");
      m_E1 = bestFit.getJackFitParValue("E");

      m_chisq = bestFit.getJackChisq();
      m_nDoF = bestFit.getNDoF();
      m_best_fit_name = best.fitname;

      int count = 1; 
      std::map<double,FitDescriptor> list = m_fits.getFitList(*m_fitComp);
      std::stringstream ss; 
      ss << "                                           | chisq/nDoF |     Q      |  fitCrit   | " << endl;

      for( std::map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++)
      {
        JackFit& thisFit = m_fits.getFit(p->second);
        double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
        double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
        ss << setw(35) <<(p->second).fitname << "|";
        ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
        ss << setw(12) << fixed << setprecision(3) << Q <<"|";
        ss << setw(12) << scientific << setprecision(3) << p->first << "|";

        if(count == rank)
        {
          ss << "*";
        }
        else
        {
          ss << " ";
        }
        ss << " FF=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("C") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("C");

        if( ((p->second).ff).operator->() == symExpPC.operator->() )
        {
          ss << ", E=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("E") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("E");
          ss << ", A=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("A") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("A");
        }
        ss << endl;count++;
      } 


      m_fit_summary = ss.str();

      std::stringstream lab;
      lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
      lab << "; FF=" << fixed << setprecision(4) << toDouble(mean(m_FF)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(m_FF)));
      std::stringstream plot;
      plot << bestFit.makeJackFitPlotAxis(t_i - 2, t_f + 2, lab.str());

      m_axis_plot = plot.str();

      AxisPlot fancy_plot = bestFit.getJackFitPlotAxis(t_i -2 , t_f + 2, lab.str());

      m_axis_plot_component =  makeFancyPlot(fancy_plot,*this,std::string("symExpPC"),t_f, t_i);

    }
  }
  else
  {
    std::cerr << "Christian has yet again proven himself an incompetent moron.." << std::endl;
    exit(1);
  }


}





void FitThreePoint::saveFitPlot(const std::string &f) const
{
  std::ofstream out(f.c_str());
  out << m_axis_plot;
  out.close();
}


