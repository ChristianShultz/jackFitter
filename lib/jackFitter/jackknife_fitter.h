#ifndef __JACKKNIFE_FITTER_H__
#define __JACKKNIFE_FITTER_H__

//jackknife wrapper around a minuit fit
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h" 
#include "Minuit2/MnUserParameterState.h" 
#include "Minuit2/MnPrint.h" 
#include "Minuit2/MnMigrad.h" 
#include "Minuit2/MnMinos.h" 
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnContours.h" 
#include "Minuit2/MnPlot.h"

#include <vector>
#include <algorithm>
#include <sstream>
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>

//#include <omp.h>

#include <itpp/itbase.h>
#include "ensem/ensem.h"
#include "adat/handle.h"

#include "ensem_data.h"
#include "linear_algebra.h"
#include "plot.h"


using namespace ENSEM;
using namespace ADAT;
using namespace std;
using namespace ROOT::Minuit2;

//*******************************************************************
//Fit Function
//*******************************************************************
class FitFunction{ //THIS IS A BASE CLASS 
  public:
    FitFunction(int nPars_);   
    virtual ~FitFunction(){};

    //return the function value with a fixed set of params - user specifies this in a concrete class
    virtual double operator()(const vector<double>& pars, double x) const = 0; //pure virtual
    //name of this type of fit
    virtual string getFitType() const = 0;


    int getNPars() const {return nPars;}
    //return the number of NON-FIXED params 
    int getNUnfixedPars() const;

    // PARAMETER NAMING
    void setParNames(vector<string> names);
    void setParName(int parNum, string name);
    vector<string> getParNames() const { return parNames; };
    string getParName(int parNum) const;
    int getParNum(string name) const;

    // DEFAULT VALUES (e.g. fit start values)
    void setDefaultParValues(vector<double> values);
    void setDefaultParValue(int parNum, double value);
    void setDefaultParValue(string name, double value);
    vector<double> getDefaultParValues() const { return defaultParValues; };
    double getDefaultParValue(int parNum) const;
    double getDefaultParValue(string name) const;

    // DEFAULT ERRORS (e.g. fit start errors)
    void setDefaultParErrors(vector<double> errors);
    void setDefaultParError(int parNum, double value);
    void setDefaultParError(string name, double value);
    vector<double> getDefaultParErrors() const { return defaultParErrors; };
    double getDefaultParError(int parNum) const;
    double getDefaultParError(string name) const;

    // PARAM FIXING
    void fixParam(int parNum);
    void fixParam(string name);
    void releaseParam(int parNum);
    void releaseParam(string name);

    bool isParamFixed(int parNum) const;
    bool isParamFixed(string name) const;

    // PARAM RANGE FIXING
    void setParamUpperLimit(int parNum, double value);
    void setParamUpperLimit(string name, double value);
    void setParamLowerLimit(int parNum, double value);
    void setParamLowerLimit(string name, double value);
    void setParamLimits(int parNum, double low, double high);
    void setParamLimits(string name, double low, double high);

    void releaseParamLimits(int parNum);
    void releaseParamLimits(string name);

    double getParamUpperLimit(int parNum) const;
    double getParamUpperLimit(string name) const;
    double getParamLowerLimit(int parNum) const;
    double getParamLowerLimit(string name) const;

    bool isParamUpperLimited(int parNum) const;
    bool isParamUpperLimited(string name) const;
    bool isParamLowerLimited(int parNum) const;
    bool isParamLowerLimited(string name) const;

  private:
    int nPars;
    vector<string> parNames;
    vector<double> defaultParValues;
    vector<double> defaultParErrors;
    vector<bool> parFixed;

    vector<double> parLowLimits;
    vector<double> parHighLimits;
    vector<bool> parUpperLimited;
    vector<bool> parLowerLimited;
};

//  CJS -- invent a failure class and create on the fly to get rid of 
//  the annoying segfault failure errors 
struct FitFunctionFailure 
: public FitFunction 
{
  FitFunctionFailure(void) : FitFunction(0) {}

  inline double operator()(const std::vector<double> &, double x) const {return -1.;}
  std::string getFitType(void) const {return std::string("FAILURE");}
};

//*******************************************************************
//chisquare class for Minuit
//*******************************************************************
class ChiSquare : public FCNBase {        
  public: 
    ChiSquare(vector<double> x_data_, vector<double> y_data_,  itpp::mat invcov_, Handle<FitFunction> ff_) 
      : x_data(x_data_), y_data(y_data_), invcov(invcov_), ff(ff_), theErrorDef(1.) {}
    ~ChiSquare() {}

    double Up() const {return theErrorDef;}
    double operator()(const vector<double>&) const; //takes the fit parameters in a vector

    void SetErrorDef(double def){theErrorDef = def;} //default 1. - see constructor

  private:
    vector<double> x_data; 
    vector<double> y_data;
    itpp::mat invcov;
    Handle<FitFunction> ff;
    double theErrorDef;
};

//************************************************************************
//perform a jackknifed set of fits
//returning an ensemble of fit params
//can use this for correlators or form-factor Qsq dependence or ...
//************************************************************************
class JackFit
{
  public:
    JackFit(EnsemData data_, Handle<FitFunction> ff_); //copies of the data  - fit function is virtual 
    ~JackFit(){};

    //run a fit on a particular (scaled) bin
    bool runBinFit(int bin, vector<double> startVals, vector<double> startErrs);

    //run a fit on the avg corr data with default start params
    bool runAvgFit();

    //run a jackknifed fit - will use results of average fit as starts
    vector<bool> runJackFit();


    vector<double> getAvgFitParValues() const {return avgParValues;};
    vector<double> getAvgFitParErrors() const {return avgParErrors;};
    double getAvgFitParValue(int n) const;
    double getAvgFitParValue(string name) const;
    double getAvgFitParError(int n) const;
    double getAvgFitParError(string name) const;

    double getAvgChisq() const {return avgChisq;};
    bool getAvgFitSuccess() const {return avgFitSuccess;};
    string getAvgFitReport() const {return avgFitReport;};


    vector<EnsemReal> getJackFitParValues() const; 
    EnsemReal getJackFitParValue(int n) const;
    EnsemReal getJackFitParValue(string name) const;

    double getJackChisq() const; //the average over the chisq on each bin

    int getNFailedFits() const;
    vector<int> getFailedFitBins() const;
    string getJackFitReport(int bin) const;
    bool isJackFitDone() const {return initJack;};
    bool getJackFitSuccess() const;

    //how many DoF for this fit?  includes reduction for reset SV & fixed params
    int getNDoF() const {return data.getNData() - data.getNResetCovSingVals() - ff->getNUnfixedPars();};

    //make a plot 
    string makeJackFitPlotAxis(EnsemFunction& weightFn, double xmin, double xmax, string label); 
    string makeJackFitPlotAxis(double xmin, double xmax, string label); 

    // CJS -- return a plot object for further editing 
    AxisPlot getJackFitPlotAxis(EnsemFunction& weightFn, double xmin, double xmax, std::string label);
    AxisPlot getJackFitPlotAxis(double xmin, double xmax, std::string label);

    // CJS -- added a vector of parameters that can be used to apply a bias function
    // in the Fit comparator routines
    void setBiasParameters(const std::vector<double> &pars);
    std::vector<double> getBiasParameters(void) const;
    void set_named_ints(const std::vector<std::pair<std::string,int> > &pars); 
    std::vector<std::pair<std::string,int> > get_named_ints(void) const;  

    // backdoor 
    ADAT::Handle<FitFunction> get_fit_func(void) {return ff;}

    // CJS -- having problems with jack fitter yet again, just crack it open 
    // and make the whole class public..
    //  private:

    EnsemData data;
    Handle<FitFunction> ff;

    //average results
    vector<double> avgParValues;
    vector<double> avgParErrors;
    double avgChisq;
    bool avgFitSuccess;
    string avgFitReport;

    //jack results
    bool initJack;
    vector<EnsemReal> scaledJackParValues;
    //vector by bin
    vector<double> jackChisqs;
    vector<bool> jackFitSuccess;
    vector<string> jackFitReports;

    // CJS -- added a vector of parameters that can be used to apply a bias function
    // in the Fit comparator routines
    bool m_have_bias_pars;
    std::vector<double> m_bias_parameters;
    std::vector<std::pair<std::string,int> > named_ints;  
};


//********************************************************************
// Method to choose the best fit from a set of attempts where
// 1. the fit function may vary, e.g. which pars are fixed, or even what function we are using
// 2. the data included may vary

//   can only use one initial EnsemData object with points turned on or off

//   FitComparator is a base class - concrete versions have an overloaded operator()(FitDescriptor fitDesc, JackFit fit) that returns a fit criterion
//   often this will just be the inverse of the chisq, but it might be funky things like 'splitN' ...



class FitDescriptor{
  public:
    FitDescriptor(Handle<FitFunction> ff_, vector<bool> activeData_, string fitname_)
      : ff(ff_), activeData(activeData_), fitname(fitname_){};

    //  FitFunction& getFF(){return *ff;}

    Handle<FitFunction> ff;  
    vector<bool> activeData;
    string fitname; 
};

//this base class has derived classes that can return a number describing the goodness of fit 
class FitComparator{
  public:
    FitComparator(){};
    virtual ~FitComparator() {}
    virtual double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const = 0 ;

    // allow for the possibility of tacking on a bias function
    virtual std::string biasFunctionName(void) const {return std::string("none");}
};

//this is an example where the inverse of the chisq / nDoF is returned (to be maximised)
class CompareFitsByChisqPerNDoF : public FitComparator{
  public:
    CompareFitsByChisqPerNDoF() : FitComparator(){};

    inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const;
};
double CompareFitsByChisqPerNDoF::operator()(const FitDescriptor& fitDesc, const JackFit& fit) const{
  return double(fit.getNDoF()) / fit.getAvgChisq();
};


//this is an example where the Q of the fit is returned (to be maximised)
class CompareFitsByQ : public FitComparator{
  public:
    CompareFitsByQ() : FitComparator(){};

    inline double operator()(const FitDescriptor& fitDesc, const JackFit& fit) const;
};
double CompareFitsByQ::operator()(const FitDescriptor& fitDesc, const JackFit& fit) const{
  return statQ( fit.getAvgChisq(), fit.getNDoF() );
};

// other possibilites, specific to whatever fit you are doing, can be added elsewhere




class JackFitLog{
  public:
    JackFitLog(EnsemData data_);

    EnsemData& getEnsemData(); //get a reference to the EnsemData to show/hide data points

    void addFit(string fitname, Handle<FitFunction> ff);  //runs an AVERAGE fit and inserts it in the map

    // allow for passing in an extra set of external bias parameters 
    void addFit(std::string fitname, ADAT::Handle<FitFunction> ff, const std::vector<double> &biasParameters);

    void addFit(std::string fitname, ADAT::Handle<FitFunction> ff, const std::vector<double> &biasParameters,
        const std::vector<std::pair<std::string,int> > &named_ints);

    void removeFit(FitDescriptor fitDesc); //removes a fit from the map

    map<double, FitDescriptor> getFitList(FitComparator& fitComp) const; //get the ordered list

    JackFit& getFit(FitDescriptor fitDesc); //return a reference to the fit - you can mess around with it

    FitDescriptor getBestFit(FitComparator& fitComp); //returns which average fit was best

    FitDescriptor getBestJackFit(FitComparator& fitComp, int& rank); //will walk down the list of fits until one is found where the jackFit succeeds

    // <CJS>  -- why do we have private members? i need to get at this for a hack 
//  private:
    // </CJS> 
    EnsemData data; //the data store
    vector<FitDescriptor> keys; 
    vector<JackFit> fits;        // a map would be better but I had problems!
};









#endif
