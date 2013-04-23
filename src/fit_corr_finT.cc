#include "jackFitter/fit_correlators.h"
#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/fit_forms.h"

#include "ensem/ensem.h"
#include "adat/handle.h"

#include <ostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>



using namespace std;
using namespace ENSEM;
using namespace ADAT;

// fit Aforw * exp( -E*t) + Aback * exp( -E* (T-t) ) 
//      + CA*exp(-E1A*T)*exp( (E1A - E1B)*t )
//      + CB*exp(-E1B*T)*exp( (E1B - E1A)*t )
// using fixed E1A, E1B
// varying Aforw, Aback, CA, CB, E

int main(int argc, char *argv[]){

  if(argc != 11){ cerr << "usage: <filename> <T> <E1A> <E1B> <tmin> <tmax> <noise_cutoff> <fitCrit> <min_tslices> <covar_SV_cutoff>" << endl; exit(1); }
  string filen; {istringstream val(argv[1]); val >> filen;}

  int T; {istringstream val(argv[2]); val >> T;}
  double E1A; {istringstream val(argv[3]); val >> E1A;}
  double E1B; {istringstream val(argv[4]); val >> E1B;}

  int tmin; {istringstream val(argv[5]); val >> tmin;}
  int tmax; {istringstream val(argv[6]); val >> tmax;}
  double cutoff; {istringstream val(argv[7]); val >> cutoff;}
  string fitCrit; {istringstream val(argv[8]); val >> fitCrit;}
  int minTSlices; {istringstream val(argv[9]); val >> minTSlices;}
  double SVcutoff;  {istringstream val(argv[10]); val >> SVcutoff;}


  /*  if((fitCrit != "chisq_per_dof") && (fitCrit != "Q") && (fitCrit != "splitN") && (fitCrit != "generic") && (fitCrit != "QN"))
    { cerr << "fit criterion " << fitCrit << " is not known" << endl; 
      cerr << "known fitCrit values are : chisq_per_dof, Q, splitN, generic, QN" << endl;  
      exit(1);
    }
  */
 

  //load the data
  EnsemVectorReal corr;  {ostringstream filename;   filename << filen;   read(filename.str(), corr);}
  cout << "read file: " << filen << endl;
  int bins = peekObs(corr, 0).size();  cout << "bins = " << bins << endl;

  //make the EnsemData object
  EnsemVectorReal lambda; lambda.resize(bins); lambda.resizeObs(tmax - tmin + 1);
  vector<double> tslices;
  for(int t = tmin; t <= tmax; t++){
    pokeObs(lambda, peekObs(corr, t), t - tmin);
    tslices.push_back(t);
  }
  EnsemData corrData(tslices, lambda);

  corrData.setSVCutoff(SVcutoff);
  cout << "covariance will be inverted cutting off SV below " << corrData.getSVCutoff() << endl; 





  //make the possible fit comparators
  map<string, Handle<FitComparator> > dum;
  
  CompareFitsByChisqPerNDoF* comp1 = new CompareFitsByChisqPerNDoF; Handle<FitComparator> comp1H(comp1);
  dum.insert( make_pair("chisq_per_dof", comp1H) );
  CompareFitsByQ* comp2 = new CompareFitsByQ; Handle<FitComparator> comp2H(comp2);
  dum.insert( make_pair("Q", comp2H) );
  CompareFitsByQN* comp5 = new CompareFitsByQN; Handle<FitComparator> comp5H(comp5);
  dum.insert( make_pair("QN", comp5H) );
  


  //do the fit

  //  FitCorrelatorCoshAndConst fitCorr(corrData, T, dum[fitCrit], cutoff, minTSlices); 
  FitCorrelatorFinT fitCorr(corrData, T, E1A, E1B, dum[fitCrit], cutoff, minTSlices); 

  //output some business
  //fit summary 
  cout << endl << fitCorr.getFitSummary() << endl << endl;
  
  //the final result
  cout << "**********************************************************" << endl;
  //  cout << "   BEST FIT uses " << fitCorr.getNExp() << " coshes" << endl;
  cout << "   chisq/ndof = " << setprecision(2) << fitCorr.getChisq() << "/" << fitCorr.getNDoF() << " = " << double(fitCorr.getChisq() /  fitCorr.getNDoF() ) << endl;
  cout << "   " << fitCorr.getFitName() << endl;


  cout << "   Aforw = " << setprecision(4) << toDouble(mean(fitCorr.Aforw)) << " +/- " << toDouble(sqrt(variance(fitCorr.Aforw))) << endl;
  cout << "   Aback = " << setprecision(4) << toDouble(mean(fitCorr.Aback)) << " +/- " << toDouble(sqrt(variance(fitCorr.Aback))) << endl;
  cout << "   CA    = " << toDouble(mean(fitCorr.CA)) << " +/- " << toDouble(sqrt(variance(fitCorr.CA))) << endl;
  cout << "   CB    = " << toDouble(mean(fitCorr.CB)) << " +/- " << toDouble(sqrt(variance(fitCorr.CB))) << endl;
  cout << "   E2    = " << toDouble(mean(fitCorr.E2)) << " +/- " << toDouble(sqrt(variance(fitCorr.E2))) << endl;
  cout << "**********************************************************" << endl;

   //write the log to a file
  {
    stringstream s; s << filen << "_finT_fit.log"; 
    ofstream out; out.open(s.str().c_str());
    out << fitCorr.getFitSummary();
    out.close();
  }
  //write the plot to a file
  {  
    stringstream ss; ss << filen << "_finT_fit.ax"; 
    ofstream out; out.open(ss.str().c_str());
    out << fitCorr.getFitPlotString();
    out.close();
  }

  //write the jackknife fit files
  {  ostringstream filename;  filename << filen << "_finT_fit_Aforw.jack";  write(filename.str(), fitCorr.Aforw );  } 
  {  ostringstream filename;  filename << filen << "_finT_fit_Aback.jack";  write(filename.str(), fitCorr.Aback );  } 
  {  ostringstream filename;  filename << filen << "_finT_fit_CA.jack";  write(filename.str(), fitCorr.CA );  } 
  {  ostringstream filename;  filename << filen << "_finT_fit_CB.jack";  write(filename.str(), fitCorr.CB );  } 
  {  ostringstream filename;  filename << filen << "_finT_fit_E2.jack";  write(filename.str(), fitCorr.E2 );  } 



}
