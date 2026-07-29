/* ****************************************************************************
 *  EventDataAlgorithm -
 *
 * \section intro_sec Overview
 *
 * Simple class for processing the events and calculating the scattering angles, DCA, and other data
 * intended for further processing in python or R.
 * 
 * 
 * @author Dennis Mackin
 * @date Feb. 29, 2016
 */

// C++ Includes
#include <sstream>

// Custom Includes
#include "EventDataAlgorithm.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

Image2D EventDataAlgorithm::getImagePlane(size_t dimension, float depth) const{
    Image2D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});

    return I;
}

Image3D EventDataAlgorithm::getImageVolume(size_t dimension) const{
    Image3D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    
    return I;
}

string EventDataAlgorithm::getDataAsString() const{
    return "No density data";
}

string EventDataAlgorithm::getConicInformationAsString() const{
    stringstream ss;
    ss << "E1,x1,y1,z1,E2,x2,y2,z2,E3,x3,y3,z3,E,theta1,theta2,alpha,phi,dca,dca_x,dca_y,dca_z,pca_x,pca_y,pca_z,px,py,pz" << endl;
    ss.precision(7);

    string delimiter = ",";
    for(size_t i = 0; i < conic_sections_.size(); ++i){

        auto scatter_info = conic_sections_[i].getScatterInfo();
//        auto comptonscatter = conic_sections_[i].getComptonScatter();

        auto positions = scatter_info->getScatterPositions();
        vector<float> energies = {scatter_info->getScatter1EnergyDeposit(), scatter_info->getScatter2EnergyDeposit(), scatter_info->getScatter3EnergyDeposit()};
        for( auto j = 0; j < 3; ++j){
            ss <<energies[j] << delimiter;
            ss << positions[j].x << delimiter << positions[j].y << delimiter << positions[j].z << delimiter;
        }
        ss  << scatter_info->getGammaEnergy() << delimiter;

        ss << scatter_info->getTheta1Degrees() << delimiter << scatter_info->getTheta2Degrees() << delimiter;
        ss << conic_sections_[i].getAlpha() * 180/M_PI << delimiter << conic_sections_[i].getPhi() * 180/M_PI;

        auto dca = conic_sections_[i].getDistanceToPoint(p_dca_);
        auto pca = conic_sections_[i].getPointOfClosestApproach(p_dca_);
        ss << delimiter << dca << delimiter << p_dca_.x << delimiter << p_dca_.y << delimiter << p_dca_.z
                        << delimiter << pca.x << delimiter << pca.y << delimiter << pca.z;

        auto point = conic_sections_[i].getLikelyOrigin();
        ss  << delimiter << point.x << delimiter << point.y << delimiter << point.z << "\n";
    }

    return ss.str();
}


void EventDataAlgorithm::run(){
    cout<<"EventDataAlgorithm::run() has nothing to do . . ."<<endl;
    return;
}
