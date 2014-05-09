/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : fit_jack_threepoint.cc

 * Purpose :

 * Creation Date : 01-02-2013

 * Last Modified : Tue 06 May 2014 11:30:11 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


/*
   currently supports fitting a constant or a double exp plus const -- if you want something
   else then add it to three_point_fit_forms.h or just make your own
 */

#include "jackFitter/three_point_fit_forms.h"
#include "jackFitter/ensem_data.h"
#include "ensem/ensem.h"
#include <iostream>
#include <sstream>
#include <vector>



void doUsage(const int argc, char * argv[])
{ 
  if(argc != 11)
  {
    std::cerr << "error: usage: " << argv[0] << " <file.jack> <0-real , 1-imag> <svCutoff> <minTSlice> <fit_type> " 
      << "<fitComparator> <tlow> <thigh> <tsrc> <tsnk>" << std::endl;
    exit(2113); 
  }
}


int main(int argc , char * argv[])
{
  doUsage(argc,argv); 

  std::string dat; {std::istringstream val(argv[1]); val >> dat;}
  int isrl; {std::istringstream val(argv[2]); val >> isrl;}
  double svCutoff; {std::istringstream val(argv[3]); val >> svCutoff;}
  int minTSlice; {std::istringstream val(argv[4]); val >> minTSlice;}

  ThreePointComparatorProps_t fitProps;
  {std::istringstream val(argv[5]); val >> fitProps.fit_type;}
  {std::istringstream val(argv[6]); val >> fitProps.baseProp;}
  fitProps.biasProp = std::string("noBiasFunction"); 
  {std::istringstream val(argv[7]); val >> fitProps.tlow;}
  {std::istringstream val(argv[8]); val >> fitProps.thigh;}
  int tsrc, tsnk;
  {std::istringstream val(argv[9]); val >> tsrc;}
  {std::istringstream val(argv[10]); val >> tsnk;}


  // get data
  //
  ENSEM::EnsemVectorComplex corr2;  
  ENSEM::read(dat,corr2);
  ENSEM::EnsemVectorReal corr;

  if ( isrl == 0 ) 
    corr = ENSEM::real(corr2) ; 
  else
    corr = ENSEM::imag(corr2) ; 

  std::cout << "read file " << dat << std::endl;
  int nb = peekObs(corr,0).size();
  std::cout << "bins = " << nb << std::endl;

  
  std::vector<double> time;
  for(int t = 0; t < corr.numElem(); ++t)
    time.push_back(double(t));

  EnsemData foobar(time,corr);
  foobar.setSVCutoff(svCutoff);

  ADAT::Handle<FitComparator> fitComp = constructThreePointFitComparator(fitProps);

  FitThreePoint fit(foobar, tsnk, tsrc, fitProps.thigh,fitProps.tlow, fitComp, minTSlice, fitProps.fit_type);

  std::cout << fit.getFitSummary() << std::endl; 

  // save plots
  //
  fit.saveFitPlot( dat + std::string(".fit.ax") ) ; 

  std::string grrr = dat + std::string(".cfit.ax"); 
  std::ofstream out ( grrr.c_str() ) ; 
  out << fit.getFitPlotStringWithComponents(); 
  out.close();  

  return 0; 
}
