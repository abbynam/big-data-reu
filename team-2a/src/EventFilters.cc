#include "EventFilters.h"

// Standard C++ Includes
#include <cmath>
#include <valarray>
#include <vector>
#include <stdexcept>
#include <iostream>


//PG Includes
#include "PhantomVolume.h"
#include "ReconstructionParabola.h"
#include "ReconstructionEllipse.h"
#include "ConicSection.h"
#include "EventsLoader.h"
#include "DoubleScatter.h"
#include "TripleScatter.h"
#include "utilities/FileUtils.h"

using namespace std;
using namespace prompt_gamma_reconstruction;


/*** Define the event filters ***/
pair<int, string> prompt_gamma_reconstruction::find_random_in_phantom_filter(shared_ptr<ConicSection> cs, pg_tools::RunTimeParameters const *params){

    auto num_tries = params->get_int("NUM_TRIES_FOR_RANDOM");
    PGVector3 p;
    auto random_points_tried = cs->getRandomPointInPhantom(p, num_tries);
    /***** DSM 2017-09-26 turned off to allow events with scatter angle near 90 degrees.
     *          This will make finding a point slower.*/
    if(-1 == random_points_tried){
        return std::make_pair<int, string>(700, "failed to find random in phantom");
    }
    /*******/
    return std::make_pair<int, string>(0, "pass");
};



pair<int, string> prompt_gamma_reconstruction::DCA_filter(shared_ptr<ConicSection> cs, pg_tools::RunTimeParameters const *params){

    PGVector3 dca_cut_point;
    dca_cut_point.x = params->get_float("DCA_CENTER_X");
    dca_cut_point.y = params->get_float("DCA_CENTER_Y");
    dca_cut_point.z = params->get_float("DCA_CENTER_Z");
    float dca_cut = params->get_float("DCA_CUT");

    auto dca = cs->getDistanceToPoint(dca_cut_point);
    if(dca > dca_cut){
        return std::make_pair<int, string>(600, "Failed DCA cut");
    }

    return std::make_pair<int, string>(0, "pass");
};

pair<int, string> prompt_gamma_reconstruction::scatter_distance_filter(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const *params){
    auto positions = tmpScatter->getScatterPositions();

    if(positions[0].getDistanceToPoint(positions[1]) < params->get_float("SCATTER_DISTANCE")){
        return std::make_pair<int, string>(500, "scatter_distance_too_small");
    }
    return std::make_pair<int, string>(0, "pass");
};

pair<int, string> prompt_gamma_reconstruction::nan_scattering_angle_filter(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const *params){

    auto scattering_angle = tmpScatter->getConeOpeningAngle();
    if( scattering_angle != scattering_angle){ //scattering_angle in not a number
        return std::make_pair<int, string>(400, "nan_scattering_angle");
    }
    return std::make_pair<int, string>(0, "pass");
};

pair<int, string> prompt_gamma_reconstruction::energy_window_filter(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const *params){

    auto e_min = params->get_float("MIN_GAMMA_ENERGY");
    auto e_max = params->get_float("MAX_GAMMA_ENERGY");
    auto e = tmpScatter->getGammaEnergy();

    if(e < e_min || e > e_max){
        stringstream errmsg;
        errmsg << "energy outside range ["<< e_min << ", " << e_max <<"].";
        return std::make_pair<int, string>(300, errmsg.str());
    }
    return std::make_pair<int, string>(0, "pass");
};

pair<int, string> prompt_gamma_reconstruction::energy_lost_filter(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const *params){

    float max_el = params->get_float("MAX_ENERGY_LOST");
    float el = tmpScatter->getEnergyLost();
    if( el*el > max_el*max_el){ //scattering_angle in not a number, skip it
        stringstream errmsg;
        errmsg << "Energy lost > ["<< max_el <<"].";
        return std::make_pair<int, string>(100, errmsg.str());
    }
    return std::make_pair<int, string>(0, "pass");
};

pair<int, string> prompt_gamma_reconstruction::scattering_angle_filter(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const *params){

    float max_angle = params->get_float("MAX_SCATTERING_ANGLE");
    float min_angle = params->get_float("MIN_SCATTERING_ANGLE");
//    auto scattering_angle = tmpScatter->getConeOpeningAngle();
    auto scattering_angle = tmpScatter->getTheta1Degrees();
    if( scattering_angle > max_angle ) {
        stringstream errmsg;
        errmsg << "scattering angle > max [" << max_angle << "].";
//        cout<< "rejecting scatter angle " << scattering_angle << " > " << max_angle << endl;
        return std::make_pair<int, string>(100, errmsg.str());
    }

    if( scattering_angle < min_angle ) {
        stringstream errmsg;
        errmsg << "scattering angle <  min [" << min_angle << "].";
//        cout<< "rejecting scatter angle " << scattering_angle << " < " << min_angle << endl;
        return std::make_pair<int, string>(100, errmsg.str());
    }

    return std::make_pair<int, string>(0, "pass");
};