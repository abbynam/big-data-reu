#ifndef DOUBLE_SCATTER_H_
#define DOUBLE_SCATTER_H_
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

//PG includes
#include "PGVector3.h"
#include "Scatter.h"

#define MASS_ELECTRON_MEV 0.510998910 //mass in MeV

using namespace std;
namespace prompt_gamma_reconstruction{

    /*! \brief Data structure for the  2 scatter Compton Camera Events
     * 
     * @author Dennis Mackin
     */
    class  DoubleScatter: public Scatter {
      private:

        double energy_deposit_[2];

        double E0_; /// calculated value for the prompt gamma energy
        double dE_[2]; /// energy deposited; subclasses may include smearing


      protected:


      public:
        DoubleScatter(const double e1, const double e2, const PGVector3 pos1, const PGVector3 pos2);
        virtual ~DoubleScatter();

        void applyPositionUncertainty(PGVector3 uncertainty, size_t  position_number);
        void applyEnergyUncertainty(double alpha, double beta, size_t  energy_dep_number);

        inline PGVector3 getAxis1() const{return (p_[0] - p_[1]).normalize();};

        inline float getTheta1() const{
            auto E = getGammaEnergy();
            auto theta = acos(1.0 + MASS_ELECTRON_MEV*(1.0/E - 1.0/(E - dE_[0])));
            return theta;
        };

        float getTheta2() const{ return -100; };
        inline float getTheta1Degrees() const {return getTheta1() * 180.0/M_PI;};

        PGVector3 getConeAxis() const{
            if(M_PI/2.0 >  getTheta1()){
                return getAxis1();
            }else{
                return getAxis1() * -1.0;
            }
        };

        float getGammaEnergy() const { return dE_[0] + dE_[1];};
        float getScatter1EnergyDeposit() const {return dE_[0];};
        float getScatter2EnergyDeposit() const {return dE_[1];};
        float getScatter3EnergyDeposit() const {return 0.0;};
        PGVector3 getConeApex() const{return p_[0];};
        void print() const;
    };
}
#endif //DOUBLE_SCATTER_
