#include <iostream>
#include <vector>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVector2.h"
#include "HHKinFitBprime.h"



using HHKinFit2::HHKinFitBprime;

int main(int argc, char* argv[])
{
  TLorentzVector bjet(11.8452,-105.891,-39.7275,113.83);
  TLorentzVector wlep(140.574,-732.856,-74.6357,754.237);
 
  std::vector<TLorentzVector> whad_jets;
  whad_jets.push_back(TLorentzVector(-88.0436,412.585,-93.1624,438.205));
  whad_jets.push_back(TLorentzVector(-39.0579,178.675,61.0409,193.553));

  std::vector<double> whad_jet_errors;
  whad_jet_errors.push_back(70.);
  whad_jet_errors.push_back(70.);

  HHKinFitBprime bprimefit(bjet, wlep, whad_jets,whad_jet_errors,70.,70.);

  bprimefit.fit();

  std::cout<<"Initial B Mass "<<bprimefit.initial_Bprime().M()<<std::endl;
  std::cout<<"Final B Mass "<<bprimefit.final_Bprime().M()<<std::endl;

  return (0);
}
