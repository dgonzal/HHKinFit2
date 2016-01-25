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
 

  bjetFit->setInitStart(100);
  bjetFit->setInitPrecision(10);
  whad_list[0]->setInitStart(100);
  whad_list[0]->setInitPrecision(10);

  fitObject->addFitObjectE(bjetFit);
  fitObject->addFitObjectE(whad_list[0]);

  fitObject->addConstraint(c_toplep_mass);
  fitObject->addConstraint(c_whad_mass);
  fitObject->addConstraint(c_bjet);
  fitObject->addConstraint(c_jet1);
  fitObject->addConstraint(c_jet2);

  try{
    fitObject->fit();
  }
  catch(HHLimitSettingException const& e){
    std::cout << e.what() << std::endl;
  }
  catch(HHKinFit2::HHEnergyRangeException const& e){
     std::cout << e.what() << std::endl;
  }
 
  initialHH = (TLorentzVector)bprime->getInitial4Vector();
  finalHH = (TLorentzVector)bprime->getFit4Vector();

  //initial vectors
  TLorentzVector wlep_ini = wlepFit->getInitial4Vector();
  TLorentzVector bjet_ini = bjetFit->getInitial4Vector();  
  TLorentzVector whad_jet1_ini = whad_list[0]->getInitial4Vector();
  TLorentzVector whad_jet2_ini = whad_list[1]->getInitial4Vector();
  TLorentzVector top_ini = top_lep->getInitial4Vector();
  TLorentzVector whad_ini = w_had->getInitial4Vector();

  //fit vectors
  TLorentzVector wlep_fit = wlepFit->getFit4Vector();
  TLorentzVector bjet_fit = bjetFit->getFit4Vector();  
  TLorentzVector whad_jet1_fit = whad_list[0]->getFit4Vector();
  TLorentzVector whad_jet2_fit = whad_list[1]->getFit4Vector();
  TLorentzVector top_fit = top_lep->getFit4Vector();
  TLorentzVector whad_fit = w_had->getFit4Vector();

  std::cout<<"Initial"<<std::endl;
  std::cout<<"W lep: Mass "<<wlep_ini.M()<<" "; wlep_ini.Print();
  std::cout<<"b jet: Mass "<<bjet_ini.M()<<" "; bjet_ini.Print();
  std::cout<<"W had jet 1: Mass "<<whad_jet1_ini.M()<<" "; whad_jet1_ini.Print();
  std::cout<<"W had jet 2: Mass "<<whad_jet2_ini.M()<<" "; whad_jet2_ini.Print();
  std::cout<<"Top: Mass "<<top_ini.M()<<" "; top_ini.Print();
  std::cout<<"Whad: Mass "<<whad_ini.M()<<" "; whad_ini.Print();

  std::cout<<"Fitted"<<std::endl;
  std::cout<<"W lep: Mass "<<wlep_fit.M()<<" "; wlep_fit.Print();
  std::cout<<"b jet: Mass "<<bjet_fit.M()<<" "; bjet_fit.Print();
  std::cout<<"W had jet 1: Mass "<<whad_jet1_fit.M()<<" "; whad_jet1_fit.Print();
  std::cout<<"W had jet 2: Mass "<<whad_jet2_fit.M()<<" "; whad_jet2_fit.Print();
  std::cout<<"Top: Mass "<<top_fit.M()<<" "; top_fit.Print();
  std::cout<<"Whad: Mass "<<whad_fit.M()<<" "; whad_fit.Print();

}


