/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : three_point_xml.cc

* Purpose :

* Creation Date : 28-05-2014

* Last Modified : Wed 28 May 2014 10:19:10 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "three_point_xml.h"

// xml junk

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

  template<typename T>
    void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f, const T &d)
    {
      if(ptop.count(path) > 0)
        read(ptop,path,place);
      else
      {
        place = d; 
      }
    }

} // namespace anonomyous 



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



std::string toString(const FitParValue::ParElem::elem &e)
{
  std::stringstream ss; 
  ss << ( (e.use) ? "T" :  "F" ); 
  ss << "  v " << e.value; 
  return ss.str(); 
}

std::ostream& operator<<(std::ostream &o, const FitParValue::ParElem::elem &e)
{
  return ( o << toString(e) ); 
}

void read(ADATXML::XMLReader &xml, const std::string &path, FitParValue::ParElem::elem &e)
{
  ADATXML::XMLReader ptop(xml,path);
  doXMLRead(ptop,"use",e.use,__PRETTY_FUNCTION__,false); 
  doXMLRead(ptop,"val",e.value,__PRETTY_FUNCTION__,0.); 
}



std::string toString(const FitParValue::ParElem &e)
{
  std::stringstream ss; 
  ss << e.parname << 
     ":default val " + toString(e.defaultParValue) 
    + "\ndefault err " + toString(e.defaultParError)
    + "\nl limit " + toString(e.paramLowerLimit) 
    + "\nu limit " + toString(e.paramUpperLimit); 
  return ss.str(); 
}

std::ostream& operator<<(std::ostream &o, const FitParValue::ParElem &e)
{
  return ( o << toString(e) ); 
}

void read(ADATXML::XMLReader &xml, const std::string &path, FitParValue::ParElem &e)
{
  FitParValue::ParElem::elem def; 
  def.use = false; 
  def.value = 0.;

  ADATXML::XMLReader ptop(xml,path);
  doXMLRead(ptop,"parname",e.parname,__PRETTY_FUNCTION__); 
  doXMLRead(ptop,"value",e.defaultParValue,__PRETTY_FUNCTION__,def); 
  doXMLRead(ptop,"error",e.defaultParError,__PRETTY_FUNCTION__,def); 
  doXMLRead(ptop,"llimit",e.paramLowerLimit,__PRETTY_FUNCTION__,def); 
  doXMLRead(ptop,"ulimit",e.paramUpperLimit,__PRETTY_FUNCTION__,def); 

} 


std::string toString(const FitParValue &e)
{
  std::stringstream ss; 
  for(int i = 0; i < e.props.size(); ++i)
    ss << e.props[i] << "\n";
  return ss.str(); 
}

std::ostream& operator<<(std::ostream &o, const FitParValue &e)
{
  return (o << toString(e)); 
}

void read(ADATXML::XMLReader &xml, const std::string &path, FitParValue &e)
{
  ADATXML::XMLReader ptop(xml,path); 
  doXMLRead(ptop,"props",e.props,__PRETTY_FUNCTION__); 
}

