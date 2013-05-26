/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : fit_jack_threepoint.cc

 * Purpose :

 * Creation Date : 01-02-2013

 * Last Modified : Mon Feb  4 10:30:02 2013

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
  if(argc != 8)
  {
    std::cerr << "error: usage: " << argv[0] << " <file.jack> <svCutoff> <minTSlice> <fit_type> " 
      << "<fitComparator> <tlow> <thigh> " << std::endl;
    exit(2113); 
  }
}


int main(int argc , char * argv[])
{
  
  std::string dat; {std::istringstream val(argv[1]); val >> dat;}
  double svCutoff; {std::istringstream val(argv[2]); val >> svCutoff;}
  int minTSlice; {std::istringstream val(argv[3]); val >> minTSlice;}

  ThreePointComparatorProps_t fitProps;
  {std::istringstream val(argv[4]); val >> fitProps.fit_type;}
  {std::istringstream val(argv[5]); val >> fitProps.baseProp;}
  fitProps.biasProp = std::string("noBiasFunction"); 
  {std::istringstream val(argv[6]); val >> fitProps.tlow;}
  {std::istringstream val(argv[7]); val >> fitProps.thigh;}


  // get data 
  ENSEM::EnsemVectorReal corr;
  ENSEM::read(dat,corr);
  std::cout << "read file " << dat << std::endl;
  int nb = peekObs(corr,0).size();
  std::cout << "bins = " << nb << std::endl;

  
  std::vector<double> time;
  for(int t = 0; t < corr.numElem(); ++t)
    time.push_back(double(t));

  EnsemData foobar(time,corr);
  foobar.setSVCutoff(svCutoff);

  ADAT::Handle<FitComparator> fitComp = constructThreePointFitComparator(fitProps);

  FitThreePoint fit(foobar, fitProps.thigh, fitProps.tlow, fitComp, minTSlice, fitProps.fit_type);

  



  return 0; 
}
