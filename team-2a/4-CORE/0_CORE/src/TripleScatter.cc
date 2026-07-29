#include <stdexcept>
#include <sstream>
#include <map>
#include <algorithm>
#include <cctype>

//PromptGamma includes
#include "TripleScatter.h"
#include "RandomSingleton.h"

using namespace prompt_gamma_reconstruction;



TripleScatter::TripleScatter(float e1, float e2, float e3, PGVector3 pos1, PGVector3 pos2, PGVector3 pos3){

    p_raw_[0] = pos1;
    p_raw_[1] = pos2;
    p_raw_[2] = pos3;
    applyPositionUncertainty(PGVector3(0.0,0.0,0.0), 0);
    applyPositionUncertainty(PGVector3(0.0,0.0,0.0), 1);
    applyPositionUncertainty(PGVector3(0.0,0.0,0.0), 2);

    //initialize the energy deposits to the raw values which
    // do not take into account detector effects
    dE_raw_[0] = e1;
    dE_raw_[1] = e2;
    dE_[0] = e1;
    dE_[1] = e2;
    dE_[2] = e3;

    //set the energy to values WITHOUT smearing
    E_[0] = getGammaEnergy();
    E_[1] = E_[0] - dE_raw_[0];
    E_[2] = E_[0] - dE_raw_[0] - dE_raw_[1];

};

TripleScatter::~TripleScatter(){ };


PGVector3 TripleScatter::getAxis1() const{
  return (p_[0] - p_[1]).normalize();
}

float TripleScatter::getTheta1() const{
  float E0 = getGammaEnergy();
  float term1 = 1.0/E0;
  float term2 = 1.0/(E0 - dE_[0]);
  float cos_theta1 = 1 + MASS_ELECTRON_MEV*(term1 - term2);
  float theta1 = acos(cos_theta1);

  return theta1;
};

float TripleScatter::getTheta2() const{
    PGVector3 v1((p_[0] - p_[1]).normalize());
    PGVector3 v2((p_[1] - p_[2]).normalize());
    return acos( v1.dotProductNormalized(v2) );
};

void TripleScatter::applyPositionUncertainty(PGVector3 uncertainty, size_t  position_number){
  if( position_number > 2){
    stringstream ss;
    ss<<"ERROR in TripleScatter::applyPositionUncertainty: position_number "<<position_number<<" is invalid.\n";
    ss<<"position_number should be in {0,1,2} corresponding to detectors 1, 2, and 3.\n\n";
    cout<<ss.str()<<endl;
    throw runtime_error(ss.str());
  }

  p_[position_number].x = p_raw_[position_number].x + (RandomSingleton::Instance()->getRand()  - 0.5) * uncertainty.x;
  p_[position_number].y = p_raw_[position_number].y + (RandomSingleton::Instance()->getRand()  - 0.5) * uncertainty.y;
  p_[position_number].z = p_raw_[position_number].z + (RandomSingleton::Instance()->getRand()  - 0.5) * uncertainty.z;
}


void TripleScatter::applyEnergyUncertainty(string detector_type, float uncertainty_scalar, size_t energy_dep_number){
  if( energy_dep_number > 1){
    stringstream ss;
    ss<<"ERROR in TripleScatter::applyEnergyUncertainty: energy_dep_number "<<energy_dep_number<<" is invalid.\n";
    ss<<"energy_dep_number should be in {0,1} corresponding to depositions in detectors 1, and 2.\n\n";
    cout<<ss.str()<<endl;
    throw runtime_error(ss.str());
  }

  dE_[energy_dep_number] = dE_raw_[energy_dep_number] + getEnergyDepositSmear(detector_type, uncertainty_scalar, E_[energy_dep_number]);
}


float TripleScatter::getEnergyDepositSmear(size_t detector_type, float uncertainty_scalar, float incident_energy){

  double resolution = 0.0;
  const float a1 = 2.16E-3;
  const float a2 = 1.82E-6;
  const float a3 = 1.042;

  //convert incident energy to keV
  incident_energy *= 1000.0;
  switch(detector_type){
    case (0): //HPGe use formulat from Owens, 1985

      resolution = sqrt(a1*incident_energy + a2*incident_energy*incident_energy + a3);
      break;
    case(1): //CZT, use formula from Du, He et al, Evaluation of a Compton scattering camera, 2000
      resolution = 6.0 + 0.15*sqrt(incident_energy);
      break;
    default:
      stringstream ss;
      ss<<"ERROR in TripleScatter::applyEnergyUncertainty: detector_type "<<detector_type<<" is invalid.\n";
      ss<<"detector_type should be 0 for HPGe and 1 for CZT.\n\n";
      cout<<ss.str()<<endl;
      throw runtime_error(ss.str());
    };

  ///convert FWHM to standard deviation of Gaussian
  resolution /= 2.0*sqrt(2*log(2));

  ///convert to MeV from KeV
  resolution /= 1000.0;
  return uncertainty_scalar*RandomSingleton::Instance()->getRandGaus()*resolution;
}

float TripleScatter::getEnergyDepositSmear(string detector_type, float scalar, float energy){
  map<string, int> det_types_map;
  det_types_map["germanium"] = 0;
  det_types_map["hpge"] = 0;
  det_types_map["czt"] = 1;
  std::transform(detector_type.begin(), detector_type.end(), detector_type.begin(), ::tolower);
  size_t detector_type_int = det_types_map[detector_type];
  return getEnergyDepositSmear(detector_type_int, scalar, energy);
}


float TripleScatter::getGammaEnergy() const{
  float E0 = dE_[0] + 0.5*( dE_[1] + sqrt(dE_[1]*dE_[1] + 4.0*dE_[1]*MASS_ELECTRON_MEV/(1.0 - cos(getTheta2()))));
  return E0;
}


void TripleScatter::print() const{

  cout<<"\n------- triple scatter -------"<<endl;
  cout<<"gamma energy raw:     "<< E_[0] <<endl;
  cout<<"gamma energy smeared: "<<  getGammaEnergy() <<"\n";
  cout<<"energies corrected:   "<<  dE_[0] <<", "<< dE_[1]<<"\n";
  cout<<"positions raw:        "<<  p_raw_[0].print() <<", "<< p_raw_[1].print() <<", "<<p_raw_[2].print() <<"\n";
  cout<<"positions corrected:  "<<  p_[0].print() <<", "<< p_[1].print() <<", "<<p_[2].print() <<"\n";
  cout<<"theta1 (deg): "<<  getTheta1Degrees()<<"\n";
  cout<<"theta2 (deg): "<<  getTheta2Degrees()<<"\n";
  cout<<"scattering angle (ang): "<<  getConeOpeningAngle()<<"\n";
  cout<<"cone apex: "<<getConeApex().print()<<endl;
  cout<<"cone axis: "<<getAxis1().print()<<endl;
  cout<<"------------------------------\n"<<endl;

  return;
}
//////////////////////////////////////////////////////////////////////////////////////
