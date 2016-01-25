#include <iostream>
#include <vector>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVector2.h"
#include "HHKinFitBprime.h"



using HHKinFit2::HHKinFitBprime;

int main(int argc, char* argv[])
{
  TLorentzVector bjet(147.05,53.0542,112.474,195.776);
  TLorentzVector wlep(111.846,107.545,230.64,315.932);
 
  std::vector<TLorentzVector> whad_jets;
  whad_jets.push_back(TLorentzVector(-151.35,-215.554,-538.938,600.031));
  whad_jets.push_back(TLorentzVector(-172.353,-103.085,-363.394,415.385));

  std::vector<double> whad_jet_errors;
  whad_jet_errors.push_back(70.);
  whad_jet_errors.push_back(70.);

  HHKinFitBprime bprimefit(bjet, wlep, whad_jets,whad_jet_errors,70.,70.);

  bprimefit.fit();

  std::cout<<"Initial B Mass "<<bprimefit.initial_Bprime().M()<<std::endl;
  std::cout<<"Final B Mass "<<bprimefit.final_Bprime().M()<<std::endl;

  return (0);
}
