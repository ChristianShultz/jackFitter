#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/fit_correlators.h"

#include "ensem/ensem.h"
#include "adat/handle.h"

#include <ostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
//#include <omp.h>


using namespace std;
using namespace ENSEM;
using namespace ADAT;

int main(int argc, char *argv[]){

  if(argc != 5){ cerr << "usage: <filename> <tmin> <tmax> <covar_SV_cutoff>" << endl; exit(1); }
  string filen; {istringstream val(argv[1]); val >> filen;}
  int tmin; {istringstream val(argv[2]); val >> tmin;}
  int tmax; {istringstream val(argv[3]); val >> tmax;}
  double SVcutoff;  {istringstream val(argv[4]); val >> SVcutoff;}

  //  omp_set_num_threads(1); // minuit will thread under the hood - can slow it down

  //cout << "we have access to " << omp_get_max_threads() << " threads" << endl;

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
  
  //do the fit
  FitConst fitConst(corrData, tmin, tmax); 

  //output some business
  //the final result
  cout << "**********************************************************" << endl;
  cout << "   chisq/ndof = " << setprecision(2) << fitConst.getChisq() 
       << "/" << fitConst.getNDoF() << " = " << double(fitConst.getChisq() /  fitConst.getNDoF() ) << endl;
  cout << "   c = " << setprecision(10) << toDouble(mean(fitConst.getConst())) << " +/- " << toDouble(sqrt(variance(fitConst.getConst()))) << endl;
  cout << "**********************************************************" << endl;

  //write the plot to a file
  {  
    stringstream ss; ss << filen << "_const_fit.ax"; 
    ofstream out; out.open(ss.str().c_str());
    out << fitConst.getFitPlotString();
    out.close();
  }

  //write the jackknife fit files
  EnsemReal out = fitConst.getConst();
  {  ostringstream filename;  filename << filen << "_const_fit.jack";  write(filename.str(), out );  } 



}
