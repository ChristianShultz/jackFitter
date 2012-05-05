#include "jackFitter/ensem_ensem_chisq_fitter.h"
#include "jackFitter/fit_forms.h"

//currently just a test suite
int main(){

  Linear* tmp = new Linear;
  Handle<FitFunction> line_fit(tmp );

  //generate some fake data
  int Ncfgs = 100;
  vector< pair<EnsemReal,EnsemReal> > vec;

  for(int i = 0; i < 10; i++)
  {
    itpp::RNG_reset(i);

    EnsemReal x,y; x.resize(Ncfgs); y.resize(Ncfgs);
    //    itpp::Normal_RNG normal(-2.0 + 4.0*double(i), 1.0);
    itpp::Normal_RNG normal(double(i), 1.0);
    for(int cfg = 0; cfg < Ncfgs; cfg++){
      pokeEnsem(x, Real( normal() ), cfg);
      pokeEnsem(y, Real(-2.0 + 4.0*double(i)), cfg);
      //      pokeEnsem(x, Real(i), cfg);
      //      pokeEnsem(y, Real( normal() ), cfg);
    }
    cout << "*** i = " << i << endl;
    cout << "x = " << toDouble(mean(x)) << " +/- " << toDouble(sqrt(variance(x))) <<endl;
    cout << "y = " << toDouble(mean(y)) << " +/- " << toDouble(sqrt(variance(y))) <<endl;


    pair<EnsemReal,EnsemReal> one = make_pair<EnsemReal,EnsemReal>(x,y) ;
    vec.push_back(one);
  }

  map<string, vector<pair<EnsemReal, EnsemReal> > > data;
  data["first"] = vec;
  cout << "made data map" << endl << endl;

  /*  EnsemEnsemChiSquare chisq(data, const_fit);
  cout << "made chisq object" << endl;

  vector<double> pars; pars.push_back(1.0);
  double c = chisq(pars);
  */

  line_fit->setDefaultParValue("slope", 0.25);
  line_fit->setDefaultParError("slope", 1.25);

  line_fit->setDefaultParValue("intercept", -0.25);
  line_fit->setDefaultParError("intercept", 2.25);
  EnsemEnsemFit fit(data, line_fit);
  
  min_result out = fit.minimise(true);
  
  cout << out.readout << endl;



  exit(0);
}
