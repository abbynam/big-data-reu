#ifndef COMPTON_SCATTER_H_
#define COMPTON_SCATTER_H_
#define _USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <ostream>
#include <sstream>

//PG INCLUDES
#include "PGVector3.h"
#include "Scatter.h"

namespace prompt_gamma_reconstruction{

/*! \brief Simple structure for storing rectangle vertices
 */
struct Rectangle { 
    std::pair<double,double> x_min_z_min ;
    std::pair<double,double> x_max_z_min ;
    std::pair<double,double> x_min_z_max ;
    std::pair<double,double> x_max_z_max ;  
};


/*! \brief Data for a Compton scattering event
 * 
 * Data for a Compton scattering event scattering event 
 * including the apex position and axis vector describing the 
 * possible source locations of the gamma
 * 
 * @author Dennis Mackin
 */
class ComptonScatter{

private:
  PGVector3 cone_apex_;
  PGVector3 cone_axis_phantom_frame_;
  PGVector3 true_origin_;
  double scatteringAngle_;
  double energy_deposition_1_;
  double energy_deposition_2_;
  double intitial_photon_energy_;

public:

  ComptonScatter(const PGVector3 &apex,
                const PGVector3 &axis,
                const PGVector3 &true_origin,
                const double &angle,
                const double &energy_deposition_1,
                const double &energy_deposition_2,
                const double &initial_energy);

  ComptonScatter(const PGVector3 &apex,
                const PGVector3 &axis,
                const double &angle,
                const double &energy_deposition_1,
                const double &energy_deposition_2,
                const double &initial_energy);

    ComptonScatter(Scatter & s):ComptonScatter(s.getConeApex(),s.getConeAxis(),
                                                                     s.getConeOpeningAngle(),
                                                                     s.getScatter1EnergyDeposit(),
                                                                     s.getScatter2EnergyDeposit(),
                                                                     s.getGammaEnergy())    {};

  inline double getDepostitedEnergy1() const{ return energy_deposition_1_; };
  inline void setDepostitedEnergy1(const double energy){ energy_deposition_1_ = energy; };
  inline double getDepostitedEnergy2() const{ return energy_deposition_2_; };
  inline void setDepostitedEnergy2(const double energy){ energy_deposition_2_ = energy; };
  inline double getInitialEnergy() const{ return intitial_photon_energy_; };
  inline void setEnergy3(const double energy){ intitial_photon_energy_ = energy; };

  inline double getScatteringAngle() const{ return scatteringAngle_; } ;
  inline void setScatteringAngle(const double angle){ scatteringAngle_ = angle; };

  inline const PGVector3 &getConeAxis() const{ return cone_axis_phantom_frame_; };
  inline const PGVector3 &getConeApex() const{ return cone_apex_; };
  inline const PGVector3 &getTrueOrigin() const{ return true_origin_; };

    inline double getAlpha( ) const{

        ///alpha is the angle between the phantom frame y axis and the cone axis. To calculate it
        /// we just need to calculate the dot product product between the y axis and the cone axis.
        /// The cone axis needs to be normalized in 2D.
        double y = cone_axis_phantom_frame_.y/sqrt(cone_axis_phantom_frame_.x*cone_axis_phantom_frame_.x
                                               + cone_axis_phantom_frame_.y*cone_axis_phantom_frame_.y
                                               + cone_axis_phantom_frame_.z*cone_axis_phantom_frame_.z);

        double alpha = acos(y);

        return alpha;
    };

};
}
#endif //COMPTON_SCATTER_H_
