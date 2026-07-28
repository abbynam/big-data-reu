/* ****************************************************************************
 *  KEMAlgorithm -
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
#include <algorithm>

// Custom Includes
#include "KEMAlgorithm.h"
#include "ASHDensity.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

Image2D KEMAlgorithm::getImagePlane(size_t dimension, float depth) const{
    Image2D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    
    if(0 == dimension){
        cout<<"Producing image for yz plane, z = "<< depth<< ". . ."<<endl;
    }else if(1 == dimension){
        cout<<"Producing image for xz plane, z  "<< depth<< ". . ."<<endl;
    }else if(2 == dimension){
        cout<<"Producing image for xy plane, z = "<< depth<< ". . ."<<endl;        
    }else{
        stringstream err_msg;
        err_msg <<"Invalid dimension " << dimension << "for KEMAlgorithm:getImagePlane." << endl;
        throw runtime_error(err_msg.str());
    }
    
    return I;
}


Image3D KEMAlgorithm::getImageVolume(size_t dimension) const{
    Image3D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0, 0.0, 0.0}, vector<size_t>{0, 0, 0});
    
    return I;
}

string KEMAlgorithm::get_event_record_(long event_num) const{
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


string KEMAlgorithm::getDataAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax,  size_t nz, float zmin, float zmax) const {

    cout << "WARNING: SBPAlgorithm only returns the SBPAlgorithm voxelation.  Ignoring the requested bin counts and ranges.\n";
    return getDataAsString();
}

string KEMAlgorithm::getDataAsString(size_t nx, size_t ny, size_t nz) const {

    cout << "WARNING: SBPAlgorithm only returns the SBPAlgorithm voxelation. Ignoring the requested bin counts. \n";
    return getDataAsString();
}

string KEMAlgorithm::getDataAsString() const {
    return getDataAsString(this->densities_);
};

string KEMAlgorithm::getDataAsString(const vector<float> hist) const {
    stringstream ss;
    ss.precision(7);

    ss << x_bins_ << " " << y_bins_ << " " << z_bins_ << endl;
    auto binedges_lambda = [&](float vmin, float vmax, size_t bins) {
        for (size_t i = 0; i <= bins; ++i) {
            ss << vmin + float(i) * (vmax - vmin) / float(bins) << ",";
        };
        ss << endl;
    };
    binedges_lambda(x_min_, x_max_, x_bins_);
    binedges_lambda(y_min_, y_max_, y_bins_);
    binedges_lambda(z_min_, z_max_, z_bins_);

    float d;

    for (size_t i = 0; i < z_bins_; ++i) {
        for (size_t j = 0; j < y_bins_; ++j) {
            for (size_t k = 0; k < x_bins_; ++k) {
                d =  hist[i*y_bins_*x_bins_ + j*x_bins_ + k];
                ss<<d<<",";
            }
        }
        ss << "\n";
    }
    return ss.str();
};


string KEMAlgorithm::getConicInformationAsString() const{

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


//void KEMAlgorithm::buildOriginConeArray(const vector<ConicSection> &cs, vector<OriginCone> &oc){
//    oc.resize(cs.size());
//
//    for(size_t i = 0; i < cs.size(); ++i){
//        oc[i].apex = cs[i].getScatterInfo()->getConeApex();
//        oc[i].axis = cs[i].getScatterInfo()->getConeAxis();
//        oc[i].cos_angle = cos(cs[i].getScatterInfo()->getConeOpeningAngle());
//
//        PGVector3 test_point(0.0, 0.0, 0.0);
//        float dca_new = oc[i].get_DCA(test_point);
//        float dca_old = cs[i].getDistanceToPoint(test_point);
//        if( (dca_new - dca_old)*(dca_new - dca_old) > 1.0E-0){
//            cout<< "ERROR: DCA does not agree: " << dca_new << ", " << cs[i].getDistanceToPoint(test_point) << endl;
//            cout << "Theta1:" << cs[i].getScatterInfo()->getTheta1Degrees() <<", angle:"<< oc[i].cos_angle*180/M_PI <<endl<<endl;
//        }
//    }
//};

OriginConesSoA KEMAlgorithm::getOriginConesSoA(const vector<ConicSection> &cs){

    vector<OriginCone> oc(cs.size());
    oc.resize(cs.size());
    OriginConesSoA origin_cones_soa;

    for(size_t i = 0; i < cs.size(); ++i){
        oc[i].apex = cs[i].getScatterInfo()->getConeApex();
        oc[i].axis = cs[i].getScatterInfo()->getConeAxis();
        oc[i].cos_angle = cos(cs[i].getScatterInfo()->getConeOpeningAngle());

        PGVector3 test_point(0.0, 0.0, 0.0);
        float dca_new = oc[i].get_DCA(test_point);
        float dca_old = cs[i].getDistanceToPoint(test_point);
        if( (dca_new - dca_old)*(dca_new - dca_old) > 1.0E-0){
            cout<< "ERROR: DCA does not agree: " << dca_new << ", " << cs[i].getDistanceToPoint(test_point) << endl;
            cout << "Theta1:" << cs[i].getScatterInfo()->getTheta1Degrees() <<", angle:"<< oc[i].cos_angle*180/M_PI <<endl<<endl;
        }
        origin_cones_soa.push_back(oc[i]);
    }
    return origin_cones_soa;
};

float KEMAlgorithm::getDensity(const PGVector3 p, const vector<OriginCone> &cones){
    float density = 0.0;
    float u = 0.0;

    //estimate density using Epanechinikov kernel
    for(size_t i=0; i < cones.size(); ++i){
        u = cones[i].get_DCA(p)*bandwidth_inv_;
//        assert(u >= 0.0f);
        density += (u < 1.0f)*0.75f*(1.0f - u*u);
    }
    return density;
}


void KEMAlgorithm::populate_system_matrix(const OriginConesSoA &conesSoA) {

    assert(system_matrix_scalar_ > 1.0);
    const size_t num_duplicates = system_matrix_scalar_;

    OriginConesSoA noise_cones = conesSoA;

    vector<PGVector3> centers = getBinCenters();
    cout<<"Sys mat has "<<system_matrix_.size()<< " voxels . . ." << endl;
    cout<<"centers mat has "<<centers.size()<< " voxels . . ." << endl;
    cout<<"noise cones "<<noise_cones.apex_x.size()<< " voxels . . ." << endl;

    vector<float> tmp_matrix = system_matrix_;
    std::fill(tmp_matrix.begin(), tmp_matrix.end(), 0.0);

    for(size_t i=0; i < num_duplicates; ++i) {

       std::random_shuffle(noise_cones.cos_angle.begin(), noise_cones.cos_angle.end());
       populate_density_matrix_cuda(tmp_matrix, &volume_grid_[0], noise_cones, bandwidth_);

        for(size_t j=0; j<system_matrix_.size(); ++j) system_matrix_[j] += tmp_matrix[j];
    }
    float max_intensity_scalar = 1.0/ *std::max_element(std::begin(system_matrix_), std::end(system_matrix_));
//     float max_intensity_scalar = 1.0/ (float(num_duplicates);
    cout <<"system matrix MAXVALUE, " << max_intensity_scalar << endl;
    for(size_t i = 0; i < system_matrix_.size(); ++i) {
        system_matrix_[i] *= max_intensity_scalar;
    }
};


void KEMAlgorithm::populate_density_matrix(const OriginConesSoA &originCones) {

    cout <<"Copying the density estimator . . ." << endl;
    vector<PGVector3> centers = getBinCenters();

    populate_density_matrix_cuda(densities_, &volume_grid_[0], originCones, bandwidth_);

    float max_intensity_scalar = 1.0/ *std::max_element(std::begin(densities_), std::end(densities_));
//    cout <<"MAXVALUE, " << max_intensity_scalar << endl;

    for(size_t i = 0; i < densities_.size(); ++i) {
        densities_[i] = std::max(0.0f, densities_[i]*max_intensity_scalar - system_matrix_[i] * float(1.1));
        //if(densities_[i] < 0.0) densities_[i] = 0.0;
    }

    max_intensity_scalar = 1.0/ *std::max_element(std::begin(densities_), std::end(densities_));

    cout <<"MAXVALUE, " << max_intensity_scalar << endl;
    for(size_t i = 0; i < densities_.size(); ++i) {
        densities_[i] *= max_intensity_scalar;
    }
};


void KEMAlgorithm::run(){
    cout<<"KEMAlgorithm::run() . . ."<<endl;
    populate_density_matrix(origin_cones_soa_);
    cout<<"KEMAlgorithm::run() complete . . ."<<endl;
    return;
};


void KEMAlgorithm::setConicSections(const vector<ConicSection> &conic_sections){
    conic_sections_.clear();
    for(size_t i = 0; i<conic_sections.size(); ++i){
        conic_sections_.push_back(conic_sections[i]);
    }
    cout<<"Copied "<<conic_sections.size() <<" conic sections to the KEMAlgorithm  . . ."<<endl;
}
