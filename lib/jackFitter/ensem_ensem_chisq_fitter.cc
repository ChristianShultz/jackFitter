#include "ensem_ensem_chisq_fitter.h"


double EnsemEnsemChiSquare::operator()(const vector<double>& pars) const{
  //  cout << "computing chisq" << endl;
  if(pars.size() != ff->getNPars()){cerr << "chisquared can't be computed, wrong number of params" << endl; exit(1);}

  map<string, double> chisq_by_ensemble;

  for( map<string, vector< pair<EnsemReal, EnsemReal> > >::const_iterator ensemble = data.begin(); ensemble != data.end(); ensemble++){
    //  cout << "ensemble = " << ensemble->first;

    vector<EnsemReal> diff; 
    vector<double> diff_mean; 

    int count = 0;
    int Ncfgs =  ( ((ensemble->second)[0]).first ).size();
    //cout << " has " << Ncfgs << " configs" << endl; 

    for( vector< pair<EnsemReal, EnsemReal> >::const_iterator i = (ensemble->second).begin(); i != (ensemble->second).end(); i++){

      //      cout << "data point #" << count << endl;
      EnsemReal yi = i->second;
      EnsemReal xi = rescaleEnsemDown( i->first );
      EnsemReal yth = Real(0.) * yi;

      for(int cfg = 0; cfg < Ncfgs; cfg++){
	pokeEnsem(yth , Real(  (*ff)(pars, toDouble(peekEnsem(xi,cfg)) )  ) , cfg);
      }

      yth = rescaleEnsemUp(yth);
      diff.push_back(yi - yth);
      diff_mean.push_back( toDouble(mean(diff[count])) );
      count++;
    }//next data point

    int dim = count;

    //compute the covariance of ( y_i - yth(x_i;{pars}) ) for this set of pars
    itpp::mat cov(dim,dim);
    
    for(int row = 0; row < dim; row++){
      for(int col = row; col < dim; col++){
	double covar = 0.0;
	for(int cfg = 0; cfg < Ncfgs; cfg++){
	  covar += toDouble( (peekEnsem(diff[row], cfg) - Real(diff_mean[row])) * (peekEnsem(diff[col], cfg) - Real(diff_mean[col])) );
	}
	cov(row, col) = covar;
	cov(col, row) = covar;
      }
    }
    cov = cov / (Ncfgs * (Ncfgs - 1) );

    //    cout << "covariance = " << endl << cov << endl;

    itpp::mat invcov = itpp::inv(cov); 

    //compute the chisq
    double chisq = 0.0;
    for(int i = 0; i < dim; i++)
      for(int j = 0; j < dim; j++){
	chisq += diff_mean[i] * invcov(i,j) * diff_mean[j];
      }
    
    chisq_by_ensemble[ensemble->first] = chisq;

  } //next ensemble

  double chisq = 0.0;
  for(map<string, double>::iterator it = chisq_by_ensemble.begin(); it != chisq_by_ensemble.end(); it++){
    //   cout << "ensemble " << it->first << " contributes " << it->second << " to the chisq" << endl;
    chisq += (*it).second;
  }
  
  return chisq;
  
};

//******************************************************************************

min_result EnsemEnsemFit::minimise(bool minos){
  
  MnUserParameters upar;

  for(int i = 0; i < ff->getNPars(); i++)
    { 
      upar.Add(  (ff->getParName(i)).c_str(), ff->getDefaultParValue(i), ff->getDefaultParError(i));
      cout << "-> adding parameter " << (ff->getParName(i)).c_str() << " starting at " << ff->getDefaultParValue(i)  ;
      if(ff->isParamFixed(i)){upar.Fix(i); cout << " [FIXED] ";}
      else{ cout  << " +/- " << ff->getDefaultParError(i); }
      if(ff->isParamUpperLimited(i)){upar.SetUpperLimit(i, ff->getParamUpperLimit(i)); cout << " < " << ff->getParamUpperLimit(i); };
      if(ff->isParamLowerLimited(i)){upar.SetLowerLimit(i, ff->getParamLowerLimit(i)); cout << " > " << ff->getParamLowerLimit(i); };
      cout << endl;
    }
  
  
  //create the minimiser
  MnMigrad mini(chisq, upar);
  //MINIMIZE
  FunctionMinimum min = mini();

  MnMinos minosErr(chisq, min); 

  //copy results
  min_result result; result.minos = minos;

  {  stringstream tmp; tmp << min; result.minuit_output = tmp.str();  }
  stringstream log;

  if(!min.IsValid()){ log << "MIN FAILED "; }
  log << "after " << min.UserState().NFcn() << " calls," << endl;
  result.chisq = min.UserState().Fval();
  log << "chisq = " << result.chisq;
  result.nDoF = nDoF;
  log << " for " << nDoF << " degrees of freedom" << endl;

  int N = min.UserState().VariableParameters();
  for(int i=0; i < N; i++){
    string name = min.UserState().Name(i);
    result.values.insert( pair<string, double>(name, min.UserState().Value(i) ) );
    result.errors.insert( pair<string, double>(name, min.UserState().Error(i) ) ); 
    log << name << " = " << min.UserState().Value(i) << " +/- " << min.UserState().Error(i) << endl;

    //if desired, run a minos error evaluation
    if( minos && min.IsValid() ){
      cout << "****** running minos on parameter : " << name << endl;
      pair<double, double> err = minosErr(i);
      result.minos_errors.insert( pair<string, pair<double, double> >(name, err) );
      log << "       minos errors : " << err.first << " , " << err.second << endl;
    }
  }

  //"naive" parameter correlations
  itpp::mat corr(N,N);
  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++){
      if(i==j){ corr(i,j) = 1.0; }
      else{
        corr(i,j) = (min.UserState().Covariance())(i,j) / sqrt( (min.UserState().Covariance())(i,i)*(min.UserState().Covariance())(j,j)   );
      }
    }
  result.correlation = corr;
  
  log << "correlation = " << endl << corr << endl;
  result.readout = log.str();

  return result;

};

