#ifndef SCATTER_H_
#define SCATTER_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>

//PG includes
#include "PGVector3.h"

#define MASS_ELECTRON_MEV 0.510998910 //mass in MeV/c^2

using namespace std;
namespace prompt_gamma_reconstruction{

class  Scatter {
  protected:
    vector<float> energy_deposit_;
    vector<PGVector3> p_raw_; //vector for unsmeared detector positions
    vector<PGVector3> p_; //Same as p1_raw_ unless voxelization and uncertainty are included
    vector<float> E_; /// calculated value for the prompt gamma energy
    vector<float> dE_raw_; /// energy deposited; subclasses may include smearing
    vector<float> dE_; /// energy deposited; subclasses may include smearing

    virtual PGVector3 getAxis1() const = 0;
    virtual float getTheta1() const = 0;
    virtual float getTheta2() const = 0;

  public:
    //Scatter(const float e1, const float e2, const PGVector3 pos1, const PGVector3 pos2, const PGVector3 pos3);
    Scatter():energy_deposit_(2), p_raw_(3), p_(3), E_(3), dE_raw_(3), dE_(3){/* DO NOTHING */};
    virtual ~Scatter(){};

    inline vector<PGVector3> getScatterPositions(){return p_;};
    inline vector<PGVector3> getScatterPositionsTrue(){return p_raw_;};
    inline float getTheta1Degrees() const{return getTheta1() * 180.0/M_PI;};
    inline float getTheta2Degrees() const{return getTheta2() * 180.0/M_PI;};
    inline float getConeOpeningAngle() const {return (M_PI/2.0 >  getTheta1())? getTheta1() : M_PI - getTheta1();} ;
    PGVector3 getConeAxis() { return (M_PI/2.0 >  getTheta1()) ? getAxis1() : getAxis1() * -1.0; };

    virtual float getGammaEnergy() const = 0;
    virtual float getScatter1EnergyDeposit() const = 0;
    virtual float getScatter2EnergyDeposit() const = 0;
    virtual float getScatter3EnergyDeposit() const = 0;
    virtual float getEnergyLost() const { return 0.0; };
    virtual PGVector3 getConeApex() const = 0;
    virtual void print() const = 0;
  };
}
#endif //SCATTER_
