/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : three_point_fit.cc

 * Purpose :

 * Creation Date : 28-05-2014

 * Last Modified : Thu 24 Jul 2014 05:00:15 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "three_point_fit_forms.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <algorithm> 




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



  //////////////////////////////////////////////////////////////////////
  //
  // make a plot but also show the individual components 
  std::string makeFancyPlot(const AxisPlot &pt,
      const FitThreePoint &tp,
      const std::string &mode,
      const int tf,
      const int ti)
  {

    std::cout << __func__ << " tf " << tf << " ti " << ti << std::endl;

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

    double eC,eA1,eA2,eE1,eE2;
    eC  = sqrt(ENSEM::toDouble(ENSEM::variance(C )));
    eA1 = sqrt(ENSEM::toDouble(ENSEM::variance(A1))); 
    eA2 = sqrt(ENSEM::toDouble(ENSEM::variance(A2)));
    eE1 = sqrt(ENSEM::toDouble(ENSEM::variance(E1)));
    eE2 = sqrt(ENSEM::toDouble(ENSEM::variance(E2))); 

    std::vector<double> m,mp,mm, time; 
    double min, max, mean,var, tmp;
    min = tp.dlow; 
    max = tp.dhigh; 

    double off = 0.;
    dC += off; 
    // add the form factor PLUS CONST for offset 
    mean = ENSEM::toDouble(ENSEM::mean(C)) + off; 
    var = sqrt( ENSEM::toDouble(ENSEM::variance(C)) ); 
    min = std::min(min,mean - var); 
    max = std::max(max,mean + var);
    for(int t = 0; t <= tf - ti; ++t)
    {
      time.push_back(t+ti); 
      m.push_back(mean);
      mp.push_back(mean + var);
      mm.push_back(mean - var); 
    }

    plot.addLineData(time,m,3);
    plot.addLineData(time,mm,3);
    plot.addLineData(time,mp,3);

    tmp = *std::max_element(mp.begin(),mp.end());
    max = std::max(max,tmp);
    tmp = *std::min_element(mm.begin(),mm.end());
    min = std::min(min,tmp);

    m.clear();
    mm.clear();
    mp.clear();


    // add the sink exp PLUS  FORM FACTOR SO IT LAYS CLOSE TO DATA 
    for(int t = ti; t <= tf; ++t)
    {
      dum = A1*ENSEM::exp(-E1*ENSEM::Real(double(tf-t)));
      mean = ENSEM::toDouble(ENSEM::mean(dum)) + dC;
      var = sqrt(ENSEM::toDouble(ENSEM::variance(dum)));
      m.push_back(mean);
      mm.push_back(mean - var);
      mp.push_back(mean + var);
    }

    plot.addLineData(time,m,4);
    plot.addLineData(time,mm,4);
    plot.addLineData(time,mp,4);


    //  // these things actually sum so keep track of absolute ranges
    //  std::transform(mp.begin(),mp.end(),mp.begin(),std::bind1st(std::plus<double>(),max - dC)); 

    tmp = *std::max_element(mp.begin(),mp.end());
    max = std::max(max,tmp);
    tmp = *std::min_element(mm.begin(),mm.end());
    min = std::min(min,tmp);

    m.clear();
    mm.clear();
    mp.clear();

    // add the source exp PLUS  FORM FACTOR SO IT LAYS CLOSE TO DATA 
    for(int t = ti; t <= tf; ++t)
    {
      dum = A2*ENSEM::exp(-E2*ENSEM::Real(t-ti));
      mean = ENSEM::toDouble(ENSEM::mean(dum)) +dC;
      var = sqrt(ENSEM::toDouble(ENSEM::variance(dum)));
      m.push_back(mean);
      mm.push_back(mean -var);
      mp.push_back(mean + var);
    } 

    plot.addLineData(time,m,5);
    plot.addLineData(time,mm,5);
    plot.addLineData(time,mp,5);

    // // these things actually sum so keep track of absolute ranges
    // std::transform(mp.begin(),mp.end(),mp.begin(),std::bind1st(std::plus<double>(),max -0.8*dC)); 


    tmp = *std::max_element(mp.begin(),mp.end());
    max = std::max(max,tmp);
    tmp = *std::min_element(mm.begin(),mm.end());
    min = std::min(min,tmp);

    plot.setYRange(min - 0.05 * fabs(min) ,max + 0.05 * fabs(max));
    plot.setXRange(ti - 2, tf + 2);


    double height = max - min;
    double unit = height/5.;


    std::stringstream lab;
    lab << "FF = " << std::setprecision(4) << dC << " +/- " << eC;
    plot.addLabel(tf+1,0.5*unit + min,lab.str(),3,0.8);
    lab.str("");
    lab << "A1 = " << std::setprecision(4) << dA1 << " +/- " << eA1;
    plot.addLabel(tf+1,1.5*unit + min,lab.str(),4,0.8);
    lab.str("");
    lab << "E1 = " << std::setprecision(4) << dE1 << " +/- " << eE1;
    plot.addLabel(tf+1,2.5*unit + min,lab.str(),4,0.8);
    lab.str("");
    lab << "A2 = " << std::setprecision(4) << dA2 << " +/- " << eA2;
    plot.addLabel(tf+1,3.5*unit + min,lab.str(),5,0.8);
    lab.str("");
    lab << "E2 = " << std::setprecision(4) << dE2 << " +/- " << eE2;
    plot.addLabel(tf+1,4.5*unit + min,lab.str(),5,0.8);

    return plot.getAxisPlotString(); 

  }


  // code cleanup 
  namespace 
  {


    //////////////////////////////////////////////////////////////////////
    //
    // callback black magic 
    ADAT::Handle<FitFunction> createThreePointConstant(const int tf , const int ti)
    {
      return ADAT::Handle<FitFunction>( new ThreePointConstant() ); 
    }

    ADAT::Handle<FitFunction> createThreePointDoubleExpPlusConst(const int tf, const int ti)
    {
      return ADAT::Handle<FitFunction>( new ThreePointDoubleExpPlusConst(tf,ti) ); 
    }

    ADAT::Handle<FitFunction> createThreePointSymmetricExpPlusConst(const int tf, const int ti)
    {
      return ADAT::Handle<FitFunction>( new ThreePointSymmetricExpPlusConst(tf,ti) ); 
    }

    ADAT::Handle<FitFunction> createThreePointLeftExpPlusConst(const int tf, const int ti)
    {
      return ADAT::Handle<FitFunction>( new ThreePointLeftExpPlusConst(tf,ti) ); 
    }

    ADAT::Handle<FitFunction> createThreePointRightExpPlusConst(const int tf, const int ti)
    {
      return ADAT::Handle<FitFunction>( new ThreePointRightExpPlusConst(tf,ti) ); 
    }


    //////////////////////////////////////////////////////////////////////
    //
    // poop out functions 
    ADAT::Handle<FitFunction> 
      function_factory(const std::string &s, const int tf, const int ti)
      {
        std::map<std::string, ADAT::Handle<FitFunction> (*)(const int, const int) > func_map;  
        func_map["const"] =  createThreePointConstant;
        func_map["dExpPC"] = createThreePointDoubleExpPlusConst;
        func_map["symExpPC"] = createThreePointSymmetricExpPlusConst; 
        func_map["leftExpPC"] = createThreePointLeftExpPlusConst;
        func_map["rightExpPC"] = createThreePointRightExpPlusConst;

        std::map<std::string, ADAT::Handle<FitFunction> (*)(const int, const int) >::const_iterator it;  

        it = func_map.find(s); 
        if ( it == func_map.end() )
        {
          std::cerr << __PRETTY_FUNCTION__ <<": error, unknown type " 
            << s << ", try one of the following" << std::endl; 

          for (it = func_map.begin(); it != func_map.end(); ++it)
            std::cout << it->first << std::endl; 
        }

        return (it->second)(tf,ti); 
      }


    //////////////////////////////////////////////////////////////////////
    //
    // initialize the start values for a given function type
    ADAT::Handle<FitFunction>
      init_three_point_fit_function(const std::string &s, EnsemData &data, const int t_f, const int t_i)
      {

        ADAT::Handle<FitFunction> func = function_factory(s,t_f,t_i);  

        double tmp,tmp2; 
        double C, El, Eh, Al, Ah;

        ENSEM::Real eC, eEl, eEh, eAl, eAh, et; 

        // pick C as the midpoint of the fit range
        tmp = double(t_f - t_i)/2. + t_i; 
        eC = ENSEM::mean(data.getYUsingNearestX(tmp)); 

        // log of derivative to get NEGATIVE mass, 
        // compute the slope to get sign of amplitude
        tmp = t_i + 1.;
        tmp2 = tmp + 1; 

        eEl = ENSEM::mean( ENSEM::log ( data.getYUsingNearestX(tmp) / data.getYUsingNearestX(tmp2) ) ) ; 
        et = ENSEM::Real(tmp);
        eAl = ENSEM::mean( data.getYUsingNearestX(tmp) / ENSEM::exp( eEl * et ) );  

        // finite difference slope
        if (  ENSEM::toDouble( ENSEM::mean( data.getYUsingNearestX(tmp) - data.getYUsingNearestX(tmp2) ) ) < 0 ) 
          eAl = eAl * ENSEM::Real(-1.);

        // same but sign flip 
        tmp = t_f -1; 
        tmp2 = tmp - 1; 
        eEh = ENSEM::mean( ENSEM::log ( data.getYUsingNearestX(tmp2) / data.getYUsingNearestX(tmp) ) ) ; 
        et = ENSEM::Real(t_f - tmp); 
        eAh = ENSEM::mean( data.getYUsingNearestX(tmp) / ENSEM::exp( eEh * et ) ); 


        // finite difference slope
        if (  ENSEM::toDouble( ENSEM::mean( data.getYUsingNearestX(tmp) - data.getYUsingNearestX(tmp2) ) ) < 0 ) 
          eAh = eAh * ENSEM::Real(-1.);


        // flip the sign of the mass terms here
        C = ENSEM::toDouble(eC); 
        El = -ENSEM::toDouble(eEl); 
        Eh = ENSEM::toDouble(eEh); 
        Al = ENSEM::toDouble(eAl); 
        Ah = ENSEM::toDouble(eAh);  

        if( func->getFitType() == "ThreePointConstant")
        {
          func->setDefaultParValue("C",C); 
        }
        else if (func->getFitType() == "ThreePointDoubleExpPlusConst")
        {
          func->setDefaultParValue("C",C);
          func->setDefaultParValue("E2",El);
          func->setDefaultParValue("E1",Eh); 
          func->setDefaultParValue("A2",Al);
          func->setDefaultParValue("A1",Ah);
        }
        else if (func->getFitType() == "ThreePointSymmetricExpPlusConst")
        {
          func->setDefaultParValue("C",C);
          func->setDefaultParValue("E",El);
          func->setDefaultParValue("A",Al); 
        }
        else if ( func->getFitType() == "ThreePointLeftExpPlusConst")
        {
          func->setDefaultParValue("C",C);
          func->setDefaultParValue("E",El);
          func->setDefaultParValue("A",Al); 
        }
        else if ( func->getFitType() == "ThreePointRightExpPlusConst")
        {
          func->setDefaultParValue("C",C);
          func->setDefaultParValue("E",Eh);
          func->setDefaultParValue("A",Ah); 
        }
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << ": error, unknown fit type " << func->getFitType() << std::endl;
          exit (1); 
        }  


        return func; 
      }


    //////////////////////////////////////////////////////////////////////
    //
    //  loop the fit over time separations 
    ADAT::Handle<FitFunction>
      loopFits( ADAT::Handle<FitFunction> func, 
          ADAT::Handle<FitComparator> &m_fitComp,
          const int t_i, 
          const int t_f, 
          const int minTSlice, 
          JackFitLog &m_fits, 
          EnsemData &data ) 
      { 

        if ( minTSlice < func->getNPars() ) 
        {
          std::cerr << __PRETTY_FUNCTION__ << ": need at least " << func->getNPars() << " to do a fit" << std::endl; 
          exit (1); 
        }


        // do a loop over all possible fit separations.. this is probably going to be slow
        for(int t_low = t_i; t_low < t_f; ++t_low)
          for(int t_high = t_f; t_high > t_low; --t_high)
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
            ss << func->getFitType() << ": t_low = " << t_low << " t_high = " << t_high; 
            // sanity 2 

            // effective NDoF
            //  if(m_fits.getEnsemData().getNData() >= minTSlice)
            m_fits.addFit(ss.str(),func,constructBiasParameters(m_fitComp->biasFunctionName(),t_f, t_i, t_high,t_low),named_ints);   
          }
        return func; 
      }

    void check_fitParVals( const ADAT::Handle<FitFunction> func,
        const FitParValue &fitParVals )
    {
      std::vector<std::string> allowed_pars; 

      if( func->getFitType() == "ThreePointConstant")
      {
        allowed_pars.push_back("C");
      }
      else if (func->getFitType() == "ThreePointDoubleExpPlusConst")
      {
        allowed_pars.push_back("C");
        allowed_pars.push_back("E2");
        allowed_pars.push_back("E1" );
        allowed_pars.push_back("A2");
        allowed_pars.push_back("A1");
      }
      else if (func->getFitType() == "ThreePointSymmetricExpPlusConst")
      {
        allowed_pars.push_back("C");
        allowed_pars.push_back("E");
        allowed_pars.push_back("A");
      }
      else if ( func->getFitType() == "ThreePointLeftExpPlusConst")
      {
        allowed_pars.push_back("C");
        allowed_pars.push_back("E");
        allowed_pars.push_back("A" );
      }
      else if ( func->getFitType() == "ThreePointRightExpPlusConst")
      {
        allowed_pars.push_back("C");
        allowed_pars.push_back("E");
        allowed_pars.push_back("A" );
      }
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << ": error, unknown fit type " << func->getFitType() << std::endl;
        exit (1); 
      }  

        std::vector<std::string>::const_iterator it; 
        for(int i = 0; i < fitParVals.props.size(); ++i)
        {
          it = std::find(allowed_pars.begin(),allowed_pars.end(),fitParVals.props[i].parname); 
          if( it == allowed_pars.end())
          {
            std::cout << __func__ << ": error for " << func->getFitType() 
              << " requested par " << fitParVals.props[i].parname
              << "\n is an unknown parameter, available parameters are:" << std::endl; 
            for(it = allowed_pars.begin(); it != allowed_pars.end(); ++it)
              std::cout << *it << std::endl;
            std::cout << "exiting.." << std::endl; 
            exit(1); 
          }
        }
    }

    // yikes, this is abstract 
    void update_fit_function(ADAT::Handle<FitFunction> &func, const FitParValue &vals)
    {
      for(int i = 0; i < vals.props.size(); ++i)
      {
        const FitParValue::ParElem *ptr = &vals.props[i]; 
        if(ptr->defaultParValue.use)
        {
          func->setDefaultParValue(ptr->parname,ptr->defaultParValue.value); 
        }
        if(ptr->defaultParError.use)
        {
          func->setDefaultParError(ptr->parname,ptr->defaultParError.value); 
        }
        if(ptr->paramLowerLimit.use)
        {
          func->setParamLowerLimit(ptr->parname,ptr->paramLowerLimit.value); 
        }
        if(ptr->paramUpperLimit.use)
        {
          func->setParamUpperLimit(ptr->parname,ptr->paramUpperLimit.value); 
        }
        if(ptr->fixParameter.use)
        {
          // overwrite here!!!
          func->setDefaultParValue(ptr->parname,ptr->fixParameter.value); 
          func->fixParam(ptr->parname); 
        }
      }
    }



    //////////////////////////////////////////////////////////////////////
    //
    // do a fit
    ADAT::Handle<FitFunction>
      tryThreePointFit( const std::string &s,
          const int tsnk, 
          const int tsrc, 
          const int t_f,
          const int t_i, 
          const int minTSlice,
          EnsemData &data,
          ADAT::Handle<FitComparator> &m_fitComp,
          JackFitLog &m_fits,
          const bool useVals, 
          const FitParValue fitParVals) 
      {
        ADAT::Handle<FitFunction> func = init_three_point_fit_function(s, data,  tsnk,  tsrc);

        // update values from default settings here 
        if( useVals )
        {
          check_fitParVals(func,fitParVals); 
          update_fit_function(func,fitParVals); 
        }

#if 1
        std::cout << __func__ << ": "  << func->getFitType() << " start pars\n"
          << "***************************************************************" 
          << std::endl; 
        std::vector<std::string> names = func->getParNames(); 
        std::vector<std::string>::const_iterator it; 
        for(it = names.begin(); it != names.end(); ++it)
        {
          std::cout << "  " << *it << "\n"
            <<"    val -> " << func->getDefaultParValue(*it) 
            << " +/- " << func->getDefaultParError(*it)
            <<"\n    llimit -> " << func->getParamLowerLimit(*it)
            <<"\n    ulimit -> " << func->getParamUpperLimit(*it)
            <<"\n    fixed  -> " << func->isParamFixed(*it)
            <<"\n" << std::endl;
        }
#endif



        return loopFits(func,m_fitComp,t_i,t_f,minTSlice,m_fits,data); 
      }

  }


  std::string 
    get_fit_summary(JackFitLog &m_fits, ADAT::Handle<FitComparator> &m_fitComp)
    {
      std::map<double,FitDescriptor> list = m_fits.getFitList(*m_fitComp);
      std::stringstream ss; 
      ss << "                                           | chisq/nDoF |     Q      |  fitCrit   | " << endl;

      for( std::map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++)
      {
        JackFit& thisFit = m_fits.getFit(p->second);
        ADAT::Handle<FitFunction> func = thisFit.get_fit_func(); 
        double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
        double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
        ss << setw(35) <<(p->second).fitname << "|";
        ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
        ss << setw(12) << fixed << setprecision(3) << Q <<"|";
        ss << setw(12) << scientific << setprecision(3) << p->first << "|";

        if(p == list.rbegin() )
        {
          ss << "*";
        }
        else
        {
          ss << " ";
        }

        for (int i = 0; i < func->getNPars(); ++i)
        {
          if (i != 0) 
            ss << ",";

          ss << " " << func->getParName(i) << setw(6) << fixed 
            << setprecision(3) << thisFit.getAvgFitParValue(func->getParName(i))
            << " +/- " << func->getParName(i) << setw(6) << fixed 
            << setprecision(3) << thisFit.getAvgFitParError(func->getParName(i));
        }
        ss << std::endl; 

      } 

      return ss.str();
    }


} // namespace anonomyous 








//////////////////////////////////////////////////////////////////////
//
FitThreePoint::FitThreePoint(EnsemData data,
    const int tsnk, 
    const int tsrc,
    const int t_f, 
    const int t_i,
    ADAT::Handle<FitComparator> fitComp, 
    const int minTSlice,
    const std::string &fit_type,
    const bool useVals, 
    const FitParValue fitParVals) 
: m_fits(data) , m_fitComp(fitComp)
{

  double stmp,si = t_i; 
  dlow = ENSEM::toDouble( ENSEM::mean( data.getYUsingNearestX(si) ) ); 
  for(si = t_i; si < t_f; ++si)
  {
    stmp = ENSEM::toDouble( ENSEM::mean( data.getYUsingNearestX(si) ) ); 
    dlow = std::min(dlow,stmp); 
    dhigh = std::max(dhigh,stmp); 
  }

  ADAT::Handle<FitFunction> dExpPC,Constant,symExpPC,leftExpPC,rightExpPC; 

  if( fit_type == "symExpPC" )
    symExpPC = tryThreePointFit("symExpPC",tsnk,tsrc,t_f, t_i,minTSlice, data, m_fitComp, m_fits,useVals,fitParVals);
  else if( fit_type == "leftExpPC" )
    leftExpPC = tryThreePointFit("leftExpPC",tsnk,tsrc,t_f, t_i,minTSlice, data, m_fitComp, m_fits,useVals,fitParVals);
  else if( fit_type == "rightExpPC" )
    rightExpPC = tryThreePointFit("rightExpPC",tsnk,tsrc,t_f, t_i,minTSlice, data, m_fitComp, m_fits,useVals,fitParVals);
  else if( fit_type == "dExpPC" )
    dExpPC = tryThreePointFit("dExpPC",tsnk,tsrc,t_f, t_i,minTSlice, data, m_fitComp, m_fits,useVals,fitParVals);
  else if( fit_type == "const" )
    Constant = tryThreePointFit("const",tsnk,tsrc,t_f, t_i,minTSlice, data, m_fitComp, m_fits,useVals,fitParVals);
  else if( fit_type == "all" ) 
  {
    symExpPC = tryThreePointFit("symExpPC",tsnk,tsrc,t_f, t_i,minTSlice, data, m_fitComp, m_fits,useVals,fitParVals);
    leftExpPC = tryThreePointFit("leftExpPC",tsnk,tsrc,t_f, t_i,minTSlice, data, m_fitComp, m_fits,useVals,fitParVals);
    rightExpPC = tryThreePointFit("rightExpPC",tsnk,tsrc,t_f, t_i,minTSlice, data, m_fitComp, m_fits,useVals,fitParVals);
    dExpPC = tryThreePointFit("dExpPC",tsnk,tsrc,t_f, t_i,minTSlice, data, m_fitComp, m_fits,useVals,fitParVals);
    Constant = tryThreePointFit("const",tsnk,tsrc,t_f, t_i,minTSlice, data, m_fitComp, m_fits,useVals,fitParVals);
  }
  else
  {
    std::cerr << __func__ << ": unknown fit type " << fit_type 
      << "options are {all , dExpPC , const, symExpPC, leftExpPC, rightExpPC }" << std::endl;
    exit(1);
  }

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
  else
  {
    m_fit_type = best.ff->getFitType();
    JackFit& bestFit = m_fits.getFit(best); 
    m_chisq = bestFit.getJackChisq();
    m_nDoF = bestFit.getNDoF();
    m_best_fit_name = best.fitname;

    m_fit_summary = get_fit_summary(m_fits, m_fitComp);


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


    // test for initialization then test for match
    bool matched = false; 

    if ( &*symExpPC)
    {  if (best.ff->getFitType() == symExpPC->getFitType())
      {
        matched = true; 
        m_FF = bestFit.getJackFitParValue("C");
        m_A1 = bestFit.getJackFitParValue("A");
        m_A2 = bestFit.getJackFitParValue("A");
        m_E1 = bestFit.getJackFitParValue("E");
        m_E2 = bestFit.getJackFitParValue("E");
      }
    }

    if ( &*dExpPC ) 
    { 
      if (best.ff->getFitType() == dExpPC->getFitType())
      {
        matched = true; 
        m_FF = bestFit.getJackFitParValue("C");
        m_A1 = bestFit.getJackFitParValue("A1");
        m_E1 = bestFit.getJackFitParValue("E1");
        m_A2 = bestFit.getJackFitParValue("A2");
        m_E2 = bestFit.getJackFitParValue("E2");
      }
    }

    if ( &*leftExpPC)
    { 
      if (best.ff->getFitType() == leftExpPC->getFitType())
      {
        matched = true; 
        m_FF = bestFit.getJackFitParValue("C");
        m_A2 = bestFit.getJackFitParValue("A");
        m_E2 = bestFit.getJackFitParValue("E");
      }
    }

    if ( &*rightExpPC )
    { 
      if (best.ff->getFitType() == rightExpPC->getFitType())
      {
        matched = true; 
        m_FF = bestFit.getJackFitParValue("C");
        m_A1 = bestFit.getJackFitParValue("A");
        m_E1 = bestFit.getJackFitParValue("E");
      }
    }

    if ( &*Constant ) 
    { 
      if (best.ff->getFitType() == Constant->getFitType())
      {
        matched = true; 
        m_FF = bestFit.getJackFitParValue("C");
      }
    }

    if ( !!! matched ) 
    {
      std::cerr << __PRETTY_FUNCTION__ << ":" << __FILE__ << ":" << __LINE__ 
        << ": something went wrong in this context" << std::endl;
      exit(1); 
    }

    // make the plot
    std::stringstream lab;
    lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
    lab << "; FF=" << fixed << setprecision(4) << toDouble(mean(m_FF)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(m_FF)));
    std::stringstream plot;
    plot << bestFit.makeJackFitPlotAxis(t_i - 2, t_f + 2, lab.str());

    m_axis_plot = plot.str();

    AxisPlot fancy_plot = bestFit.getJackFitPlotAxis(t_i -2 , t_f + 2, lab.str());

    m_axis_plot_component =  makeFancyPlot(fancy_plot,*this,m_best_fit_name,t_high, t_low);

  }// close if ! failed

}





void FitThreePoint::saveFitPlot(const std::string &f) const
{
  std::ofstream out(f.c_str());
  out << m_axis_plot;
  out.close();
}


