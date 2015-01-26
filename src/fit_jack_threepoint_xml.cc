/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : fit_jack_threepoint.cc

 * Purpose :

 * Creation Date : 01-02-2013

 * Last Modified : Sat 06 Dec 2014 01:16:01 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "jackFitter/three_point_fit_forms.h"
#include "jackFitter/ensem_data.h"
#include "ensem/ensem.h"
#include <iostream>
#include <sstream>
#include <vector>


namespace stupid_reader
{
  template<typename T>
    void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
    {
      if(ptop.count(path) > 0)
        read(ptop,path,place);
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " 
          << f << " trying to read path, " << path
          << ", path was empty, exiting" << std::endl;
        exit(1);
      }
    }

  template<typename T>
    void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f, const T &val)
    {
      if(ptop.count(path) > 0)
        read(ptop,path,place);
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " 
          << f << " trying to read path, " << path
          << ", path was empty, reverting to default value " 
          << val << std::endl;
        place = val;
      }
    }


  struct XMLStruct
  {
    int rl; 
    int tsnk; 
    int tsrc; 
    bool RealType; 
    std::string correlator; 
    ThreePointComparatorProps_t threePointProps; 
    FitParValue fitParValues; 
  };

  void read(ADATXML::XMLReader &xml, const std::string &pth, XMLStruct &p)
  {
    ADATXML::XMLReader ptop(xml,pth); 
    doXMLRead(ptop,"rl",p.rl,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"RealType",p.RealType,__PRETTY_FUNCTION__,false); 
    doXMLRead(ptop,"tsnk",p.tsnk,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"tsrc",p.tsrc,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"correlator",p.correlator,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"threePointProps",p.threePointProps,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"fitParValues",p.fitParValues,__PRETTY_FUNCTION__); 
  }

} // stupid_reader 

void doUsage(const int argc, char * argv[])
{ 
  if(argc != 2)
  {
    std::cerr << "error: usage: " << argv[0] << " <ini.xml" << std::endl; 
    exit(2113); 
  }
}


int main(int argc , char * argv[])
{
  doUsage(argc,argv); 
  std::istringstream val(argv[1]); 
  std::string xmlini; 
  val >> xmlini; 

  stupid_reader::XMLStruct mXMLStruct; 

  try
  { 
    std::cout << "reading " << xmlini << std::endl;
    ADATXML::XMLReader xml(xmlini); 
    read(xml,"/props",mXMLStruct); 
  }
  catch(...)
  {
    std::cout << "some error occurred, debug it yourself" << std::endl;
    exit(1);
  }


  // get data
  //

  ENSEM::EnsemVectorReal corr;

  // handle complex vs real types here 
  if( !!! mXMLStruct.RealType ) 
  {
    ENSEM::EnsemVectorComplex corr2;  
    ENSEM::read(mXMLStruct.correlator,corr2);
    if ( mXMLStruct.rl == 0 ) 
      corr = ENSEM::real(corr2) ; 
    else
      corr = ENSEM::imag(corr2) ; 
  }
  else
  {
    ENSEM::read(mXMLStruct.correlator,corr); 
  }


  std::cout << "read file " << mXMLStruct.correlator << std::endl;
  int nb = peekObs(corr,0).size();
  std::cout << "bins = " << nb << std::endl;


  // make the time 
  std::vector<double> time;
  for(int t = 0; t < corr.numElem(); ++t)
  {
    time.push_back(double(t));

    double foo = ENSEM::toDouble( ENSEM::mean( ENSEM::peekObs(corr,t))); 
    double foov = ENSEM::toDouble( ENSEM::variance( ENSEM::peekObs(corr,t))); 
    std::cout << t << "  " << foo << " pm " << foov << std::endl; 
  }

  std::cout << std::endl;

  std::cout << "using SVCutOff " << mXMLStruct.threePointProps.SVCutOff << std::endl;

  // make the data
  EnsemData foobar(time,corr);
  foobar.setSVCutoff(mXMLStruct.threePointProps.SVCutOff);

  // make the comparator
  ADAT::Handle<FitComparator> fitComp = constructThreePointFitComparator(mXMLStruct.threePointProps);

  // do the fit 
  FitThreePoint fit(foobar, mXMLStruct.tsnk, mXMLStruct.tsrc, 
      mXMLStruct.threePointProps.thigh,
      mXMLStruct.threePointProps.tlow, 
      fitComp, 
      mXMLStruct.threePointProps.minTSlice, 
      mXMLStruct.threePointProps.fit_type,
      true, mXMLStruct.fitParValues);

  // what did we get?
  std::cout << fit.getFitSummary() << std::endl; 

  // save plots
  //
  // skip this guy, component plot is more useful 
  // fit.saveFitPlot( mXMLStruct.correlator + std::string(".fit.ax") ) ; 

  std::string gr = mXMLStruct.correlator + std::string(".fitComp.ax"); 
  std::ofstream fooboo (gr.c_str() ); 
  fooboo << fit.getFitPlotStringWithComponents(); 
  fooboo.close(); 

  // dont really need this 
  //   std::string grrr = mXMLStruct.correlator + std::string(".fit.summary"); 
  //   std::ofstream out ( grrr.c_str() ) ; 
  //   out << fit.getFitSummary(); 
  //   out.close();  

  // save the constant 
  ENSEM::write( mXMLStruct.correlator + std::string(".FF.jack") , fit.getFF()); 

  return 0; 
}
