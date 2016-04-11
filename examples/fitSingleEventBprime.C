#include <iostream>
#include <vector>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVector2.h"
#include "HHKinFitBprime.h"



using HHKinFit2::HHKinFitBprime;

int main(int argc, char* argv[])
{
  TLorentzVector bjet(-27.246,-100.936,95.2741,141.695);
  TLorentzVector wlep(-154.898,-286.537,30.9601,336.927);
 
  std::vector<TLorentzVector> whad_jets;
  whad_jets.push_back(TLorentzVector(-10.4016,301.142,610.043,680.788));
  whad_jets.push_back(TLorentzVector(-6.13401,123.878,400.649,419.531));

  std::vector<double> whad_jet_errors;
  whad_jet_errors.push_back(7.);
  whad_jet_errors.push_back(7.);

  HHKinFitBprime bprimefit(bjet, wlep, whad_jets,whad_jet_errors,10.,10.);
  bprimefit.set_verbosity(2);
  bprimefit.fit();

  std::cout<<"Initial B Mass "<<bprimefit.initial_Bprime().M()<<std::endl;
  std::cout<<"Final B Mass "<<bprimefit.final_Bprime().M()<<std::endl;

  return (0);
}
