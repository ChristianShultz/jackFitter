#ifndef __ENSEM_ENSEM_CHISQ_FITTER_H__
#define __ENSEM_ENSEM_CHISQ_FITTER_H__

#include "jackknife_fitter.h"

//**********************************************
// take (x,y) data where both x and y are ensembles
// chisq uses a covariance that is recomputed for each fit parameter iteration
// allows for multiple ensembles to be included via a block-diagonal covariance
// fit returns parameters and error using minuit error estimation methods

// reuses the FitFunction class from jackknife_fitter.h
//***********************************************


class EnsemEnsemChiSquare : public FCNBase { //of the type required by Minuit
 public:
  EnsemEnsemChiSquare( 
		      map<string, vector< pair<EnsemReal, EnsemReal> > > data_, //string labels the ensembles 
		      Handle<FitFunction> ff_
		      ) 
    : data(data_), ff(ff_), theErrorDef(1.) 
    {
      cout << "fitting with " << ff->getFitType() << endl;
      cout << "*****************************************" << endl;
      
      for( map<string, vector< pair<EnsemReal, EnsemReal> > >::const_iterator ensemble = data.begin(); ensemble != data.end(); ensemble++){
	cout << "ensemble = " << ensemble->first << endl;
	for( vector< pair<EnsemReal, EnsemReal> >::const_iterator i = (ensemble->second).begin(); i != (ensemble->second).end(); i++){
	  cout << "x= " << toDouble(mean(i->first)) <<"+/-" << toDouble(sqrt(variance(i->first))) 
	       << ", y= "  << toDouble(mean(i->second)) <<"+/-" << toDouble(sqrt(variance(i->second)))
	       << endl;
	}
	cout << "*****************************************" << endl;
      }
      
    };
    
    double operator()(const vector<double>&) const;
    //    { cout << "Here" << endl;
    // }; //returns the chisq for a vector of fit param values

    double Up() const {return theErrorDef;};    
    void SetErrorDef(double def){theErrorDef = def;};
    
 private:
  map<string, vector< pair<EnsemReal, EnsemReal> > > data;
  Handle<FitFunction> ff;
  double theErrorDef;
};


struct min_result{
  double chisq;
  int nDoF;
  map<string, double> values;
  map<string, double> errors;
  bool minos;
  map<string, pair<double, double> > minos_errors;
  itpp::mat correlation;

  string minuit_output;
  string readout;
};


class EnsemEnsemFit{
 public:
  EnsemEnsemFit( 
		    map<string, vector< pair<EnsemReal, EnsemReal> > > data_, //string labels the ensembles 
		    Handle<FitFunction> ff_
		    )
    : data(data_), ff(ff_), chisq(data_, ff_) 
    {
      nDoF = 0;
      for( map<string, vector< pair<EnsemReal, EnsemReal> > >::const_iterator ensemble = data.begin(); ensemble != data.end(); ensemble++)
	for( vector< pair<EnsemReal, EnsemReal> >::const_iterator i = (ensemble->second).begin(); i != (ensemble->second).end(); i++){
	  nDoF++;
	}
      nDoF -= (ff->getNPars());      
    };
    
    min_result minimise(bool minos); //option to compute minos errors - slow;
    // start values etc... live inside ff

 private:
  map<string, vector< pair<EnsemReal, EnsemReal> > > data;
  Handle<FitFunction> ff;
  EnsemEnsemChiSquare chisq;
  int nDoF;
};
  
#endif
