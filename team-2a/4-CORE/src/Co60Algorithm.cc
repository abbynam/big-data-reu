/* ****************************************************************************
 *  Co60Algorithm -
 *
 * \section intro_sec Overview
 *
 * Co60 algorithm is an experimental code that attempts to correct Co60 events before reconstructing them.
 *
 * The approach is to first sort the events into 1.17 and 1.33 MeV groups and then correct the scatter energies for
 * E1 and E2 based on previously derived correction factors.
 * 
 * @author Dennis Mackin
 * @date August 08, 2016
 */

// C++ Includes
#include <sstream>

// Custom Includes
#include "Co60Algorithm.h"
#include "Scatter.h"
#include "DoubleScatter.h"
#include "utilities/Random.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

string Co60Algorithm::getDataAsString() const{
    /*algo_->setConicSections(conic_sections_);
    algo_->run();
    string dose_string = algo_->getDataAsString();

    dose_string += "~\n";

    algo_->setConicSections(conic_sections_corrected_);
    algo_->run();
    dose_string += algo_->getDataAsString();*/

    string dose_string = "Turned~Off for performance reasons";

    return dose_string;
}


shared_ptr<ConicSection> Co60Algorithm::get_conic_section_(float E1, float E2, PGVector3 &p1, PGVector3 &p2, size_t event_number){

    auto tmpScatter = make_shared<DoubleScatter>(E1, E2, p1, p2);
    auto compton_scatter = shared_ptr<ComptonScatter>(new ComptonScatter(*tmpScatter));

    auto cs = shared_ptr<ConicSection>(new ConicSection(*compton_scatter, phantom_volume_ptr_, event_number));
    cs->setScatter(tmpScatter);

    return cs;
};


float matchEventToKnownEnergy(float e1, float e2){

    const float HIGH_ENERGY = 1.33;
    const float LOW_ENERGY = 1.17;

    const float COMPTON_MAX_117 = LOW_ENERGY - (0.511 * LOW_ENERGY) / (0.511 + 2.0 * LOW_ENERGY);
    const float COMPTON_MAX_133 = HIGH_ENERGY - (0.511 * HIGH_ENERGY) / (0.511 + 2.0 * HIGH_ENERGY);

    if((e1 > COMPTON_MAX_117) && (e1 <= COMPTON_MAX_133)) return HIGH_ENERGY;

    return LOW_ENERGY;
    /*if(e1 > COMPTON_MAX_117 - (COMPTON_MAX_133 - COMPTON_MAX_117)) return LOW_ENERGY;

    //Randomize the rest for now
    pg_tools::Random rand;
//    rand.SetSeed(125);
    if(rand.Rndm() > 0.5) return HIGH_ENERGY;
    return LOW_ENERGY;*/
};


//// Correct the cones based on the know gamma energies for Co60 and the known position of origin.
vector<ConicSection> Co60Algorithm::correct_conic_sections_(vector<ConicSection> &conic_sections){

    vector<ConicSection> corrected;
//    ComptonScatter const *comp_scat;
    for(size_t i = 0; i < conic_sections.size(); ++i){

        auto e1 = conic_sections[i].getScatterInfo()->getScatter1EnergyDeposit();
        auto e2 = conic_sections[i].getScatterInfo()->getScatter2EnergyDeposit();

        if(e1 < min_scatter_energy_ || e2 < min_scatter_energy_ || e1 + e2 < min_event_energy_){
            continue;
        }

        auto scat_positions = conic_sections[i].getScatterInfo()->getScatterPositions();

        assert(scat_positions.size() >= 2);

        PGVector3 v1 = scat_positions[0] - origin_;
        PGVector3 v2 = scat_positions[1] - scat_positions[0];

        float cos_theta = v1.dotProductNormalized(v2);
        float e_known = matchEventToKnownEnergy(e1, e2);
        float e1_known = e_known*e_known*(1.0 - cos_theta)/(0.511 + e_known*(1.0 - cos_theta));
        float e2_known = e_known - e1_known;

        auto cs = get_conic_section_(e1_known, e2_known, scat_positions[0], scat_positions[1], i);

        corrected.push_back(*cs);
    }

    return corrected;
};


string Co60Algorithm::get_event_record_(const ConicSection &cs, const ConicSection &cs_corrected) const{
    stringstream ss;
    ss.precision(7);

    string delimiter = ",";

    auto scatter_info = cs.getScatterInfo();
    auto comptonscatter = cs.getComptonScatter();

    auto positions = scatter_info->getScatterPositions();
    vector<float> energies = {scatter_info->getScatter1EnergyDeposit(), scatter_info->getScatter2EnergyDeposit(), scatter_info->getScatter3EnergyDeposit()};
    for( auto j = 0; j < 3; ++j){
        ss <<energies[j] << delimiter;
        ss << positions[j].x << delimiter << positions[j].y << delimiter << positions[j].z << delimiter;
    }
    ss  << scatter_info->getGammaEnergy() << delimiter;

    ss << scatter_info->getTheta1Degrees() << delimiter << scatter_info->getTheta2Degrees() << delimiter;
    ss << cs.getAlpha() * 180/M_PI << delimiter << cs.getPhi() * 180/M_PI;

    auto dca = cs.getDistanceToPoint(p_dca_);
    auto pca = cs.getPointOfClosestApproach(p_dca_);
    ss << delimiter << dca << delimiter << p_dca_.x << delimiter << p_dca_.y << delimiter << p_dca_.z
    << delimiter << pca.x << delimiter << pca.y << delimiter << pca.z;

    auto point = cs.getLikelyOrigin();
    ss  << delimiter << point.x << delimiter << point.y << delimiter << point.z << delimiter;

    //Add in the corrected values
    scatter_info = cs_corrected.getScatterInfo();
    comptonscatter = cs_corrected.getComptonScatter();
    ss  << scatter_info->getGammaEnergy() << delimiter;
    ss  << scatter_info->getScatter1EnergyDeposit() << delimiter;
    ss  << scatter_info->getScatter2EnergyDeposit() << delimiter;
    ss << scatter_info->getTheta1Degrees() << delimiter;
    ss << cs_corrected.getAlpha() * 180/M_PI << delimiter << cs_corrected.getPhi() * 180/M_PI;

    pca = cs_corrected.getPointOfClosestApproach(p_dca_);
    ss << cs_corrected.getDistanceToPoint(p_dca_) << endl;

    return ss.str();
}


string Co60Algorithm::getConicInformationAsString() const{

    stringstream ss;
    ss.precision(7);


    string header = "E1,x1,y1,z1,E2,x2,y2,z2,E3,x3,y3,z3,E,theta1,theta2,alpha,phi,dca,dca_x,dca_y,dca_z,pca_x,pca_y,pca_z,px,py,pz,"
                    "E_k,E1_c,E2_c,theta1_c,alpha_c,phi_c,dca_c";
    ss << header << endl;

    cout<<"Genertating the records . . ."<<endl;

    //Get the uncorrected event record
    for(size_t i = 0; i < conic_sections_.size(); ++i){
        ss << get_event_record_(conic_sections_[i], conic_sections_corrected_[i]);
    }

    ss << endl;
    return ss.str();
}



void Co60Algorithm::run(){
    cout<<"Co60Algorithm::run() . . ."<<endl;
    cout<<"Co60Algorithm::run() complete . . ."<<endl;
    return;
}



Image2D Co60Algorithm::getImagePlane(size_t dimension, float depth) const{
    Image2D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});

    if(0 == dimension){
        cout<<"Producing image for yz plane, z = "<< depth<< ". . ."<<endl;
    }else if(1 == dimension){
        cout<<"Producing image for xz plane, z = "<< depth<< ". . ."<<endl;
    }else if(2 == dimension){
        cout<<"Producing image for xy plane, z = "<< depth<< ". . ."<<endl;
    }else{
        stringstream err_msg;
        err_msg <<"Invalid dimension " << dimension << "for Co60Algorithm:getImagePlane." << endl;
        throw runtime_error(err_msg.str());
    }

    return I;
}

Image3D Co60Algorithm::getImageVolume(size_t dimension) const{
    Image3D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    return I;
}
