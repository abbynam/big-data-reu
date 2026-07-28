/* ****************************************************************************
 *  OctaneAlgorithm -
 *
 * \section intro_sec Overview
 *
 * ImageAlgorithm defines the interface for the various reconstruction 
 * algorithms used in the Prompt Gamma imaging software.
 * 
 * 
 * @author Dennis Mackin
 * @date November 21, 2015
 */

// C++ Includes
#include <sstream>

// Custom Includes
#include "OctaneAlgorithm.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

Image2D OctaneAlgorithm::getImagePlane(size_t dimension, float depth) const{
    Image2D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    
    if(0 == dimension){
        cout<<"Producing image for yz plane, z = "<< depth<< ". . ."<<endl;
    }else if(1 == dimension){
        cout<<"Producing image for xz plane, z = "<< depth<< ". . ."<<endl;        
    }else if(2 == dimension){
        cout<<"Producing image for xy plane, z = "<< depth<< ". . ."<<endl;        
    }else{
        stringstream err_msg;
        err_msg <<"Invalid dimension " << dimension << "for OctaneAlgorithm:getImagePlane." << endl;
        throw runtime_error(err_msg.str());
    }
    
    return I;
}

Image3D OctaneAlgorithm::getImageVolume(size_t dimension) const{
    Image3D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    
    return I;
}


string OctaneAlgorithm::get_event_record_(long event_num) const{
    stringstream ss;
    ss.precision(7);

    string delimiter = ",";

    auto scatter_info = conic_sections_[event_num].getScatterInfo();
//    auto comptonscatter = conic_sections_[event_num].getComptonScatter();

    auto positions = scatter_info->getScatterPositions();
    vector<float> energies = {scatter_info->getScatter1EnergyDeposit(), scatter_info->getScatter2EnergyDeposit(), scatter_info->getScatter3EnergyDeposit()};
    for( auto j = 0; j < 3; ++j){
        ss <<energies[j] << delimiter;
        ss << positions[j].x << delimiter << positions[j].y << delimiter << positions[j].z << delimiter;
    }
    ss  << scatter_info->getGammaEnergy() << delimiter;

    ss << scatter_info->getTheta1Degrees() << delimiter << scatter_info->getTheta2Degrees() << delimiter;
    ss << conic_sections_[event_num].getAlpha() * 180/M_PI << delimiter << conic_sections_[event_num].getPhi() * 180/M_PI;

    auto dca = conic_sections_[event_num].getDistanceToPoint(p_dca_);
    if(dca != dca){
        stringstream ss;
        ss<<"ERROR: DCA is nan for point("<<p_dca_.x<<","<<p_dca_.y<<","<<p_dca_.z<<").\n";
        cout<<ss.str()<<endl;
        throw runtime_error(ss.str());
    }
    auto pca = conic_sections_[event_num].getPointOfClosestApproach(p_dca_);
    ss << delimiter << dca << delimiter << p_dca_.x << delimiter << p_dca_.y << delimiter << p_dca_.z
    << delimiter << pca.x << delimiter << pca.y << delimiter << pca.z;

    auto point = conic_sections_[event_num].getLikelyOrigin();
    ss  << delimiter << point.x << delimiter << point.y << delimiter << point.z << endl;

    return ss.str();
}


string OctaneAlgorithm::getConicInformationAsString() const{

    vector<string> records(conic_sections_.size());
    stringstream ss;
    ss << "E1,x1,y1,z1,E2,x2,y2,z2,E3,x3,y3,z3,E,theta1,theta2,alpha,phi,dca,dca_x,dca_y,dca_z,pca_x,pca_y,pca_z,px,py,pz" << endl;
    ss.precision(7);

    cout<<"Genertating the records . . ."<<endl;

    for(size_t i = 0; i < conic_sections_.size(); ++i){
        records[i] = get_event_record_(i);
    }

    cout<<"Adding records to string stream . . ."<<endl;
    for(size_t i = 0; i < records.size(); ++i){
        ss << records[i];
    }
    ss << endl;
    return ss.str();
}


//std::array<PGVector3, 8> OctaneAlgorithm::get_octant_centers(const PGVector3 &center, const double length){
std::vector<PGVector3> OctaneAlgorithm::get_octant_centers(const PGVector3 &center, const double length){
    double xvals [2] = {center.x - 0.5 * length, center.x + 0.5 * length};
    double yvals [2]  = {center.y - 0.5 * length, center.y + 0.5 * length};
    double zvals [2]  = {center.z - 0.5 * length, center.z + 0.5 * length};

    std::vector<PGVector3> points(8);
    for(size_t i = 0; i < 8; ++i) {
        points[i] = PGVector3(xvals[i & 1], yvals[(i & 2) >> 1], zvals[(i & 4) >> 2]);
    }

    return points;
};


//vector<PGVector3> OctaneAlgorithm::get_intercepts(const ConicSection &cs, const PGVector3 &center, float length, float intercept_dca){
//
//    //std::array<PGVector3, 8> centers = get_octant_centers(center, 0.5*length);
//    std::vector<PGVector3> centers = get_octant_centers(center, 0.5*length);
//    vector<PGVector3> intercepts;
//    intercepts.reserve(8*int(length*length/(intercept_dca*intercept_dca)));
//
//    for (auto bin_center: centers) {
//        auto dca = cs.getDistanceToPoint(bin_center);
//        if(dca < 0.5*length) {
//            if(0.5*length > intercept_dca) {
//                auto new_intercepts = get_intercepts(cs, bin_center, 0.5 * length, intercept_dca);
//                intercepts.insert(std::end(intercepts), std::begin(new_intercepts), std::end(new_intercepts));
//            }else {
//                intercepts.push_back(bin_center);
//            }
//        }
//    }
//    return intercepts;
//}

void OctaneAlgorithm::set_intercepts(const ConicSection &cs, const PGVector3 &center, float length, float intercept_dca){

//    std::array<PGVector3, 8> centers = get_octant_centers(center, 0.5*length);
    std::vector<PGVector3> centers = get_octant_centers(center, 0.5*length);
    for (auto bin_center: centers) {

        auto dca = cs.getDistanceToPoint(bin_center);
        if(dca <= 0.5*length) {
            if(0.5*length >= intercept_dca) {
                set_intercepts(cs, bin_center, 0.5 * length, intercept_dca);
            }else {
                PGVector3 pca = cs.getPointOfClosestApproach(bin_center);

                PGVector3 copy = pca;
                cs.transformPointToConeAxisFrameFromPhantomFrame(copy);

                if(copy.y < 0 ) continue;
                float weight = pow(pca.getDistanceToPoint(cs.getComptonScatter().getConeApex()), inverse_square_param_);
                if(phantom_volume_ptr_->is_in_volume(pca)){
                    density_estimator_ptr_->fill(pca, weight);
                }else{
                    density_estimator_ptr_->fill(bin_center, weight);
                }
            }

        }
    }
}



void OctaneAlgorithm::populate_density_matrix(const vector<ConicSection> &conics){

    size_t num_cones = conics.size();
    density_estimator_ptr_->clear();
    #pragma omp parallel for
    for(size_t i=0; i<num_cones; ++i){
        if( i % 1000 == 0) cout<<"Cone "<<i<<" . . . "<<endl;
        set_intercepts(conics[i], phantom_center_, phantom_length_, intercept_dca_);
    }
}

void OctaneAlgorithm::run(){
    cout<<"OctaneAlgorithm::run() . . ."<<endl;
    populate_density_matrix(conic_sections_);
    cout<<"OctaneAlgorithm::run() complete . . ."<<endl;
    return;
}


void OctaneAlgorithm::setConicSections(const vector<ConicSection> &conic_sections){
    conic_sections_.clear();
    for(size_t i = 0; i<conic_sections.size(); ++i){
        conic_sections_.push_back(conic_sections[i]);
    }
    cout<<"Copied "<<conic_sections.size() <<" conic sections to the OctaneAlgorithm  . . ."<<endl;
}
