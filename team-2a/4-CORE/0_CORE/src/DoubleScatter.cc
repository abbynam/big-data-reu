#include <cmath>
#include <stdexcept>
#include <sstream>

//PromptGamma includes
#include "DoubleScatter.h"
#include "RandomSingleton.h"

using namespace prompt_gamma_reconstruction;

DoubleScatter::DoubleScatter(const double e1, const double e2, PGVector3 pos1, const PGVector3 pos2){

    p_raw_[0] = pos1;
    p_raw_[1] = pos2;

    applyPositionUncertainty(PGVector3(0.0,0.0,0.0), 0);
    applyPositionUncertainty(PGVector3(0.0,0.0,0.0), 1);

    //initialize the energy deposits to the raw values which
    // do not take into account detector effects
    energy_deposit_[0] = e1;
    energy_deposit_[1] = e2;
    dE_[0] = e1;
    dE_[1] = e2;

    E0_ = getGammaEnergy();
};


DoubleScatter::~DoubleScatter(){};


void DoubleScatter::applyPositionUncertainty(PGVector3 uncertainty, size_t  position_number){
    if( position_number > 2){
        stringstream ss;
        ss<<"ERROR in DoubleScatter::applyPositionUncertainty: position_number "<<position_number<<" is invalid.\n";
        ss<<"position_number should be in {0,1,2} corresponding to detectors 1, 2, and 3.\n\n";
        cout<<ss.str()<<endl;
        throw runtime_error(ss.str());
    }

    p_[position_number].x = p_raw_[position_number].x + (RandomSingleton::Instance()->getRand()  - 0.5) * uncertainty.x;
    p_[position_number].y = p_raw_[position_number].y + (RandomSingleton::Instance()->getRand()  - 0.5) * uncertainty.y;
    p_[position_number].z = p_raw_[position_number].z + (RandomSingleton::Instance()->getRand()  - 0.5) * uncertainty.z;
}


void DoubleScatter::applyEnergyUncertainty(double alpha, double beta, size_t  energy_dep_number){
    if( energy_dep_number > 1){
        stringstream ss;
        ss<<"ERROR in DoubleScatter::applyEnergyUncertainty: energy_dep_number "<<energy_dep_number<<" is invalid.\n";
        ss<<"energy_dep_number should be in {0,1} corresponding to depositions in detectors 1, and 2.\n\n";
        cout<<ss.str()<<endl;
        throw runtime_error(ss.str());
    }

    double uncertainty = sqrt(alpha + beta*E0_);
    dE_[energy_dep_number] = energy_deposit_[energy_dep_number] + uncertainty * RandomSingleton::Instance()->getRandGaus();
}



void DoubleScatter::print() const{

  cout<<"\n------- double scatter -------"<<endl;
  cout<<"gamma energy raw:     "<< E0_ <<endl;
  cout<<"gamma energy smeared: "<<  getGammaEnergy() <<"\n";
  cout<<"energies raw:         "<<  energy_deposit_[0] <<", "<< energy_deposit_[1]<<"\n";
  cout<<"energies corrected:   "<<  dE_[0] <<", "<< dE_[1]<<"\n";
  cout<<"positions raw:        "<<  p_raw_[0].print() <<", "<< p_raw_[1].print() <<"\n";
  cout<<"positions corrected:  "<<  p_[0].print() <<", "<< p_[1].print() <<"\n";
  cout<<"theta1 (deg): "<<  getTheta1Degrees()<<"\n";
  cout<<"scattering angle (ang): "<<  getConeOpeningAngle()<<"\n";
  cout<<"cone apex: "<<getConeApex().print()<<endl;
  cout<<"cone axis: "<<getAxis1().print()<<endl;
  cout<<"------------------------------\n"<<endl;

  return;
}


