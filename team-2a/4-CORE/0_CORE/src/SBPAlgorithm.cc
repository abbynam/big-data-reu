/* ****************************************************************************
 *  VectorAlgorithm -
 *
 * \sectio-n intro_sec Overview
 *
 * ImageAlgorithm defines the interface for the various reconstruction
 * algorithms used in the Prompt Gamma imaging software.
 *
 *
 * @author Dennis Mackin
 * @date Feb 18, 2018
 */

// C++ Includes
#include <sstream>

// Custom Includes
#include "SBPAlgorithm.h"
#include "OriginCone.h"
#include "cuda_functions.h"

using namespace std;
using namespace prompt_gamma_reconstruction;


string SBPAlgorithm::getDataAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax,  size_t nz, float zmin, float zmax) const {

    cout << "WARNING: SBPAlgorithm only returns the SBPAlgorithm voxelation.  Ignoring the requested bin counts and ranges.\n";
    return getDataAsString();
}

string SBPAlgorithm::getDataAsString(size_t nx, size_t ny, size_t nz) const {

    cout << "WARNING: SBPAlgorithm only returns the SBPAlgorithm voxelation. Ignoring the requested bin counts. \n";
    return getDataAsString();
}

string SBPAlgorithm::getDataAsString() const {
//Do not allow ranges outside of the reconstruction volume

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

    PGVector3 p;
    float d;

    for (size_t i = 0; i < z_bins_; ++i) {
        for (size_t j = 0; j < y_bins_; ++j) {
            for (size_t k = 0; k < x_bins_; ++k) {
                d =  densities_[i*y_bins_*x_bins_ + j*x_bins_ + k];
                ss<<d<<",";
            }
        }
        ss << "\n";
    }
    return ss.str();
};


Image2D SBPAlgorithm::getImagePlane(size_t dimension, float depth) const{
    Image2D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    return I;
}

Image3D SBPAlgorithm::getImageVolume(size_t dimension) const{
    Image3D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    return I;
}


string SBPAlgorithm::get_event_record_(long event_num) const{
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


string SBPAlgorithm::getConicInformationAsString() const{

    vector<string> records(conic_sections_.size());
    stringstream ss;
    ss << "E1,x1,y1,z1,E2,x2,y2,z2,E3,x3,y3,z3,E,theta1,theta2,alpha,phi,dca,dca_x,dca_y,dca_z,pca_x,pca_y,pca_z,px,py,pz" << endl;
    ss.precision(7);

    cout<<"Generating the records . . ."<<endl;

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

vector<PGVector3> SBPAlgorithm::getBinCenters(){
    vector<PGVector3> centers(densities_.size());
    for(size_t i=0; i < centers.size(); ++i) centers[i] = this->getBinCenter(i);

    return centers;
}
//
//float SBPAlgorithm::getDensity(const PGVector3 p, const vector<OriginCone> &cones){
//    float density = 0.0;
//    float u = 0.0;
//
//    //estimate density using Epanechinikov kernel
//    for(size_t i=0; i < cones.size(); ++i){
//        u = cones[i].get_DCA(p)*bandwidth_inv_;
////        assert(u >= 0.0f);
//        density += (u < 1.0f)*0.75f*(1.0f - u*u);
//    }
//    return density;
//}


void SBPAlgorithm::populate_density_matrix(const vector<OriginCone> &cones){
    vector<PGVector3> centers = getBinCenters();
    OriginConesSoA cones_soa;
    cones_soa.add_cones(cones);
    populate_density_matrix_cuda(densities_, &volume_grid_[0], cones_soa, bandwidth_);

    float max_intensity_scalar = 1.0/ *std::max_element(std::begin(densities_), std::end(densities_));

    cout <<"MAXVALUE, " << max_intensity_scalar << endl;
    for(size_t i = 0; i < densities_.size(); ++i) {
        densities_[i] *= max_intensity_scalar;
    }
}

void SBPAlgorithm::run(){
    cout<<"SBPAlgorithm::run() . . ."<<endl;
    populate_density_matrix(cones_);
    cout<<"SBPAlgorithm::run() complete . . ."<<endl;
    return;
}


void SBPAlgorithm::setConicSections(const vector<ConicSection> &conic_sections){
    conic_sections_.clear();
    for(size_t i = 0; i<conic_sections.size(); ++i){
        conic_sections_.push_back(conic_sections[i]);
    }

    cones_ = OriginCone::build_origin_cone_array(conic_sections);

    cout<<"Copied "<<conic_sections.size() <<" conic sections to the SBPAlgorithm  . . ."<<endl;
}
