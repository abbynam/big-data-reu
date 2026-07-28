
#include <stdexcept>
#include <sstream>

#include "ComptonScatter.h"

using namespace std;
namespace prompt_gamma_reconstruction{



ComptonScatter::ComptonScatter(const PGVector3 &apex,
                              const PGVector3 &axis,
                              const PGVector3 &true_origin,
                              const double &angle,
                              const double &energy_deposition_1,
                              const double &energy_deposition_2,
                              const double &initial_energy):
                              cone_apex_(apex),
                              cone_axis_phantom_frame_(axis),
                              true_origin_(true_origin),
                              scatteringAngle_(angle),
                              energy_deposition_1_(energy_deposition_1),
                              energy_deposition_2_(energy_deposition_2),
                              intitial_photon_energy_(initial_energy){
    if( energy_deposition_1_ + energy_deposition_2_ > initial_energy*1.001){
        stringstream ss;
        ss<<"ERROR: initial energy is less than the depositions.\n";
        ss<<"intial energy = "<< initial_energy<<"\n";
        ss<<"depostion 1 = "<< energy_deposition_1<<"\n";
        ss<<"depostion 2 = "<< energy_deposition_2<<"\n\n";
        cout<<ss.str()<<endl;
        throw runtime_error(ss.str());
    }

    if( angle > M_PI){
        stringstream ss;
        ss<<"ERROR: scattering angle > PI.\n";
        ss<<"Make sure angle passed to ComptonScatter is in radians.\n\n";
        cout<<ss.str()<<endl;
        throw runtime_error(ss.str());
    }
};

ComptonScatter::ComptonScatter(const PGVector3 &apex,
                              const PGVector3 &axis,
                              const double &angle,
                              const double &energy_deposition_1,
                              const double &energy_deposition_2,
                              const double &initial_energy):
                              cone_apex_(apex),
                              cone_axis_phantom_frame_(axis),
                              scatteringAngle_(angle),
                              energy_deposition_1_(energy_deposition_1),
                              energy_deposition_2_(energy_deposition_2),
                              intitial_photon_energy_(initial_energy){
    
    //Check conservation of energy
    if( energy_deposition_1_ + energy_deposition_2_ > initial_energy*1.001){
        stringstream ss;
        ss<<"ERROR: initial energy is less than the depositions.\n";
        ss<<"intial energy = "<< initial_energy<<"\n";
        ss<<"depostion 1 = "<< energy_deposition_1<<"\n";
        ss<<"depostion 2 = "<< energy_deposition_2<<"\n\n";
        cout<<ss.str()<<endl;
        throw runtime_error(ss.str());
    }

    //Scattering angle cannot be greater than PI. The likely cause
    // is angles in degrees rather than radians. 
    if( angle > M_PI){
        stringstream ss;
        ss<<"ERROR: scattering angle > PI.\n";
        ss<<"Make sure angle passed to ComptonScatter is in radians.\n\n";
        cout<<ss.str()<<endl;
        throw runtime_error(ss.str());
    }
};

}//end of namespace
