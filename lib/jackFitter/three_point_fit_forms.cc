/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : three_point_fit_forms.cc

 * Purpose :

 * Creation Date : 09-11-2012

 * Last Modified : Fri Nov  9 17:20:51 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "three_point_fit_forms.h"
#include <iostream>
#include <sstream>
#include <map>



FitThreePoint::FitThreePoint(EnsemData data, int t_f, int t_i, ADAT::Handle<FitComparator> fitComp, int minTSlice) 
  : m_fits(data) , m_fitComp(fitComp)
{

  if(minTSlice < 5) 
  {
    std::cout << __func__ << " need at least 5 data points to do a fit" << std::endl;
    exit(1);
  }

  // try the double exp plus const fit 
  ADAT::Handle<FitFunction> dExpPC (new DoubleExpPlusConst(t_f,t_i));

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
    sign = 1.;
  else
    sign = -1.;

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

      std::stringstream ss; 
      ss << "DoubleExpPlusConst: t_low = " << t_low << " t_high = " << t_high; 

      // sanity 2 
      if(m_fits.getEnsemData().getNData() >= minTSlice)
        m_fits.addFit(ss.str(),dExpPC);   
    }

  FitDescriptor best_dExpPC = m_fits.getBestFit(*m_fitComp);


  // im copying the idea of this out of zfit code.. not sure how/why/if it works
  int rank;
  FitDescriptor best = m_fits.getBestJackFit(*m_fitComp,rank);
  if(best.fitname == "FAILED")
  {
    m_chisq = 1.e10;
    m_nDoF = 1;
    m_best_fit_name = "FAILED";
    m_fit_summary = "FAILED";
    double tip1 = double(t_i + 1);  // this reference to double crap is getting annoying..
    m_FF = data.getYUsingNearestX(tip1);
    m_FF = ENSEM::Real(0.);
    m_A1 = m_FF;
    m_E1 = m_FF;
    m_A2 = m_FF;
    m_E2 = m_FF; 
    // no plot
  }
  else if (best.ff.operator->() == dExpPC.operator->())
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
    ss << "                                   | chisq/nDoF |     Q      |  fitCrit   | " << endl;

    for( std::map<double, FitDescriptor>::reverse_iterator p = list.rbegin(); p != list.rend(); p++)
    {
      JackFit& thisFit = m_fits.getFit(p->second);
      double chisq_per_ndof = thisFit.getAvgChisq() / thisFit.getNDoF();
      double Q = statQ( thisFit.getAvgChisq() , thisFit.getNDoF() );
      ss << setw(35) <<(p->second).fitname << "|";
      ss << setw(12) << fixed << setprecision(3) << chisq_per_ndof <<"|";
      ss << setw(12) << fixed << setprecision(3) << Q <<"|";
      ss << setw(12) << scientific << setprecision(3) << p->first << "|";

      if(count == rank){ ss << "*";}else{ ss << " ";}
      ss << " FF=" << setw(8) << fixed << setprecision(4) << thisFit.getAvgFitParValue("C") << " +/-" <<  setw(8) << fixed <<setprecision(4) << thisFit.getAvgFitParError("C");
      if( ((p->second).ff).operator->() == dExpPC.operator->() ){
        ss << ", E1'=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("E1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("E1");
        ss << ", E2=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("E2") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("E2");
        ss << ", A1'=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("A1") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("A1");
        ss << ", A2=" << setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParValue("A2") << " +/-" <<  setw(6) << fixed << setprecision(3) << thisFit.getAvgFitParError("A2");

      }
      ss << endl;count++;
    } 
    //    cout << ss.str();
    m_fit_summary = ss.str();

    std::stringstream lab;
    lab << "\\gx\\sp2\\ep/N\\sbdof\\eb=" << setprecision(2) << bestFit.getJackChisq() << "/" << bestFit.getNDoF(); 
    lab << "; FF=" << fixed << setprecision(4) << toDouble(mean(m_FF)) << "\\+-" <<  setprecision(4) << toDouble(sqrt(variance(m_FF)));
    m_axis_plot = bestFit.makeJackFitPlotAxis(t_i - 2, t_f + 2, lab.str());

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


