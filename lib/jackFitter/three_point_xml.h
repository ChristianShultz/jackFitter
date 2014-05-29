#ifndef THREE_POINT_XML_H
#define THREE_POINT_XML_H 

#include "io/adat_xmlio.h"




////////////////////////////////////////////////////////////////
///////////////////     XML STUFF    ///////////////////////////
////////////////////////////////////////////////////////////////

struct ThreePointComparatorProps_t 
{
  std::string fit_type;
  std::string baseProp; 
  std::string biasProp;
  ADATXML::Array<std::string> extraProps;
  int tlow;
  int thigh; // the fit range
  int minTSlice; 
};

std::string toString(const ThreePointComparatorProps_t &);
std::ostream& operator<<(std::ostream&, const ThreePointComparatorProps_t &);
void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointComparatorProps_t &); 



////////////////////////////////////////////////////////////////
///////////////////   MINUIT STUFF    //////////////////////////
////////////////////////////////////////////////////////////////

struct FitParValue
{
  struct ParElem
  {
    struct elem
    {
      bool use; 
      double value; 
    };

    std::string parname; 
    elem defaultParValue; 
    elem defaultParError;
    elem paramLowerLimit; 
    elem paramUpperLimit; 
  };

  
  ADATXML::Array<ParElem> props;

}; 


std::string toString(const FitParValue::ParElem::elem &);
std::ostream& operator<<(std::ostream&, const FitParValue::ParElem::elem &);
void read(ADATXML::XMLReader &xml, const std::string &path, FitParValue::ParElem::elem &); 


std::string toString(const FitParValue::ParElem &);
std::ostream& operator<<(std::ostream&, const FitParValue::ParElem &);
void read(ADATXML::XMLReader &xml, const std::string &path, FitParValue::ParElem &); 


std::string toString(const FitParValue &);
std::ostream& operator<<(std::ostream&, const FitParValue &);
void read(ADATXML::XMLReader &xml, const std::string &path, FitParValue &); 


#endif /* THREE_POINT_XML_H */
