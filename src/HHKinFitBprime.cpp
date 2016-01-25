#include "HHKinFitBprime.h"
#include "HHKinFit.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectEConstBeta.h"
#include "HHFitObjectMET.h"
#include "HHFitObjectComposite.h"
#include "HHFitConstraintLikelihood.h"
#include "HHFitConstraint.h"
#include "HHFitConstraintEHardM.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectEConstBeta.h"
#include "HHFitConstraint4Vector.h"
#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHLimitSettingException.h"
#include "exceptions/HHCovarianceMatrixException.h"

#include "TMatrixD.h"
#include "TRandom3.h"

#include <TMath.h>
#include <cmath>
#include <cstdlib>
#include <iterator>

HHKinFit2::HHKinFitBprime::HHKinFitBprime(TLorentzVector const& bjet,
					  TLorentzVector const& Wlep,
					  std::vector<TLorentzVector> const& Whad_jets,
					  std::vector<double> sigma_Whad_jets,
					  double sigma_Wlep,
					  double sigma_bjet){

  m_bjet = HHLorentzVector(bjet);
  m_Wlep = HHLorentzVector(Wlep);
  for(auto jet : Whad_jets)
    m_Whad_jets.push_back(HHLorentzVector(jet));
  
  m_Wlep.SetMkeepE(80);
}


void HHKinFit2::HHKinFitBprime::fit(){
  //vector of fit variables
  HHFitObjectE* bjetFit = new HHFitObjectEConstBeta(m_bjet);
  HHFitObjectE* wlepFit = new HHFitObjectEConstM(m_Wlep);
  std::vector<HHFitObjectE*> whad_list;
  for(auto jet : m_Whad_jets)
    whad_list.push_back(new HHFitObjectEConstBeta(jet));
  
  //prepare composite objects
  HHFitObject* top_lep = new HHFitObjectComposite(bjetFit,wlepFit);
  HHFitObject* w_had = new HHFitObjectComposite(whad_list[0],whad_list[1]);
  //for(auto jet : whad_list)
  //  w_had->addSubobject(jet);
  HHFitObject* bprime = new HHFitObjectComposite(top_lep,w_had);

  try{
    whad_list[0]->setLowerFitLimitE(1);
    whad_list[1]->setLowerFitLimitE(1);
    whad_list[0]->setUpperFitLimitE(5000);
    whad_list[1]->setUpperFitLimitE(5000);
  }
  catch(HHLimitSettingException const& e){
     std::cout << "Exception while setting W_{had} - jet limits" << std::endl;
     std::cout << e.what() << std::endl;
  }

  try{
    wlepFit->setLowerFitLimitE(1);
    bjetFit->setLowerFitLimitE(1);
    wlepFit->setUpperFitLimitE(5000);
    bjetFit->setUpperFitLimitE(5000);
  }
  catch(HHLimitSettingException const& e){
    std::cout << "Exception while setting top limits:" << std::endl;
    std::cout << e.what() << std::endl;
  }

  //prepare constraints
  HHFitConstraint* c_toplep_mass = new HHFitConstraintEHardM(bjetFit,wlepFit, 175);
  HHFitConstraint* c_whad_mass = new HHFitConstraintEHardM(whad_list[0],whad_list[1], 80);

  HHFitConstraint* c_bjet = new HHFitConstraint4Vector(bjetFit, false, false,false, true);
  HHFitConstraint* c_jet1 = new HHFitConstraint4Vector(whad_list[0], false, false,false, true);
  HHFitConstraint* c_jet2 = new HHFitConstraint4Vector(whad_list[1], false, false,false, true);

  //fit
  HHKinFit2::HHKinFit* fitObject = new HHKinFit2::HHKinFit();
  
  fitObject->addFitObjectE(bjetFit);
  fitObject->addFitObjectE(whad_list[0]);

  fitObject->addConstraint(c_toplep_mass);
  fitObject->addConstraint(c_whad_mass);
  fitObject->addConstraint(c_bjet);
  fitObject->addConstraint(c_jet1);
  fitObject->addConstraint(c_jet2);

  fitObject->fit();

  initialHH = (TLorentzVector)bprime->getInitial4Vector();
  finalHH = (TLorentzVector)bprime->getFit4Vector();
}


