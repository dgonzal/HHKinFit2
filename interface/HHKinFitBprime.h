#ifndef HHKinFitBprime_H
#define HHKinFitBprime_H

#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVector2.h>

#ifdef HHKINFIT2
#include "HHLorentzVector.h"
#include "HHKinFit.h"
#else
#include "HHKinFit2/HHKinFit2/interface/HHLorentzVector.h"
#include "HHKinFit2/HHKinFit2/interface/HHKinFit.h"
#endif

#include <stdio.h>
#include <map>
#include <utility>
#include <vector>
#include <sstream>

namespace HHKinFit2{

  class HHKinFitBprime{
  public:
    HHKinFitBprime(TLorentzVector const& bjet,
		   TLorentzVector const& Wlep,
		   std::vector<TLorentzVector> const& Whad_jets,
		   std::vector<double> sigma_Whad_jets,
		   double sigma_Wlep,
		   double sigma_bjet);
    
    //the main action, runs over all hypotheses and performs the fit
    void fit();
    
    //Getters
    double getBestChi2();
  
  private:
    HHLorentzVector m_bjet,m_Wlep;
    std::vector<HHLorentzVector> m_Whad_jets;
    double m_bestChi2 = 99999;
    
    TLorentzVector initialHH,finalHH;
  };
}
#endif
