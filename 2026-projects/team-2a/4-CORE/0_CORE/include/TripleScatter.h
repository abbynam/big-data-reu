#ifndef TRIPLE_SCATTER_H_
#define TRIPLE_SCATTER_H_
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif // M_PI

#ifndef M_PI
#define MASS_ELECTRON_MEV 0.510998910 //mass in MeV/c^2
#endif

#include <cmath>


//PG includes
#include "PGVector3.h"
#include "Scatter.h"
#include "utilities/Random.h"


using namespace std;
namespace prompt_gamma_reconstruction{

/*! \brief Data structure for the 3 scatter Compton Camera Events
 * 
 * Applies smearing to simulated detector effects in addition
 * to storing the scattering information for a triple scatter event.
 * 
 * @author Dennis Mackin
 */
class  TripleScatter: public Scatter {
  private:
    //pg_tools::Random rand_;

  public:
    TripleScatter(float e1, float e2, float e3, PGVector3 pos1, PGVector3 pos2, PGVector3 pos3);
    virtual ~TripleScatter();

    //Apply uncertainty function are used to simulate the finite detector resolution.
    void applyPositionUncertainty(PGVector3 uncertainty, size_t  position_number);
    void applyEnergyUncertainty(string detector_type, float scalar, size_t energy_dep_number);
    
    //functions return the energy smearing value
    float getEnergyDepositSmear(string detector_type, float scalar, float incident_energy);
    float getEnergyDepositSmear(size_t detector_type, float scalar, float incident_energy);


    PGVector3 getAxis1() const;
    float getTheta1() const;
    float getTheta2() const;
    float getTheta1Degrees() const {return getTheta1() * 180.0/M_PI;};
    float getTheta2Degrees() const {return getTheta2() * 180.0/M_PI;};
    float getConeOpeningAngle() const {return (M_PI/2.0 >  getTheta1())? getTheta1() : M_PI - getTheta1();};
    inline PGVector3 getConeAxis() const{
        if(M_PI/2.0 >  getTheta1()) {
          return getAxis1();
        }else{
          return getAxis1() * -1.0;
        };
    };

    float getGammaEnergy() const;
    inline float getScatter1EnergyDeposit() const {return dE_[0];};
    inline float getScatter2EnergyDeposit() const {return dE_[1];};
    inline float getScatter3EnergyDeposit() const {return dE_[2];};
    inline float getMeasuredEnergy() const {return dE_[0] + dE_[1] + dE_[2];};
    inline float getEnergyLost() const { return getGammaEnergy() - getMeasuredEnergy(); };
    inline PGVector3 getConeApex() const {return p_[0];};
    void print() const;
  };
}
#endif //TRIPLE_SCATTER_
