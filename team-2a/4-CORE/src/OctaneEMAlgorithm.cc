/* ****************************************************************************
 *  OctaneEMAlgorithm -
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
#include "OctaneEMAlgorithm.h"
#include "ASHDensity.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

Image2D OctaneEMAlgorithm::getImagePlane(size_t dimension, float depth) const{
    Image2D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    
    if(0 == dimension){
        cout<<"Producing image for yz plane, z = "<< depth<< ". . ."<<endl;
    }else if(1 == dimension){
        cout<<"Producing image for xz plane, z = "<< depth<< ". . ."<<endl;        
    }else if(2 == dimension){
        cout<<"Producing image for xy plane, z = "<< depth<< ". . ."<<endl;        
    }else{
        stringstream err_msg;
        err_msg <<"Invalid dimension " << dimension << "for OctaneEMAlgorithm:getImagePlane." << endl;
        throw runtime_error(err_msg.str());
    }
    
    return I;
}

Image3D OctaneEMAlgorithm::getImageVolume(size_t dimension) const{
    Image3D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    
    return I;
}


string OctaneEMAlgorithm::get_event_record_(long event_num) const{
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


string OctaneEMAlgorithm::getConicInformationAsString() const{

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


std::vector<PGVector3> OctaneEMAlgorithm::get_octant_centers(const PGVector3 &center, const double length){
    double xvals [2] = {center.x - 0.5 * length, center.x + 0.5 * length};
    double yvals [2]  = {center.y - 0.5 * length, center.y + 0.5 * length};
    double zvals [2]  = {center.z - 0.5 * length, center.z + 0.5 * length};

    std::vector<PGVector3> points(8);
    for(size_t i = 0; i < 8; ++i) {
        points[i] = PGVector3(xvals[i & 1], yvals[(i & 2) >> 1], zvals[(i & 4) >> 2]);
    }

    return points;
};


vector<PGVector3> OctaneEMAlgorithm::get_intercepts(const OriginCone &originCone, const PGVector3 &center, float length, float intercept_dca){

    //std::array<PGVector3, 8> centers = get_octant_centers(center, 0.5*length);
    std::vector<PGVector3> centers = get_octant_centers(center, 0.5*length);
    vector<PGVector3> intercepts;
    intercepts.reserve(8*int(length*length/(intercept_dca*intercept_dca)));
//
//    const PGVector3 scatter1 = cs.getScatterInfo()->getScatterPositions()[0];
//    const PGVector3 scatter2 = cs.getScatterInfo()->getScatterPositions()[1];
//    const float scatter_angle = cs.getScatterInfo()->getConeOpeningAngle();

    for (auto bin_center: centers) {
        auto dca = originCone.get_DCA(bin_center);
        if(dca <= 0.5*length) {
            if(0.5*length >= intercept_dca) {
                auto new_intercepts = get_intercepts(originCone, bin_center, 0.5 * length, intercept_dca);
                intercepts.insert(std::end(intercepts), std::begin(new_intercepts), std::end(new_intercepts));
            }else {
                auto pca = originCone.get_PCA(bin_center);
                if(phantom_volume_ptr_->is_in_volume(pca)) {
                    intercepts.push_back(pca);
                }
            }
        }
    }
    return intercepts;
}


void OctaneEMAlgorithm::buildOriginConeArray(const vector<ConicSection> &cs, vector<OriginCone> &oc){
    oc.resize(cs.size());

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

    }

};

void OctaneEMAlgorithm::populate_system_matrix(const vector<OriginCone> &originCones) {
/**
 * Builds a system matrix from the inverse distance of the first scatters to each voxel in the matrix
 *
 * This factor combines at least 2 effects: the relative angle subtended by the detector for each voxel,
 * proportional to the inverse of the distance squared and the cone surface area,proportional to the distance.
 * @see (none yet)
 *
 * @param conics Vector containing the detected conic sections
 * @return void.
 */
    size_t num_cones = originCones.size();
    const size_t num_duplicates = 10;

//    vector<OriginCone> noise_cones(originCones);
    vector<OriginCone> noise_cones(num_cones*num_duplicates);
    for(size_t i=0; i < noise_cones.size(); ++i) noise_cones[i] = originCones[i % num_cones];
    std::random_shuffle(noise_cones.begin(), noise_cones.end());
    for(size_t i=0; i < noise_cones.size(); ++i) {

//        float dca_before = noise_cones[i].get_DCA(PGVector3(0,0,0));
//        cout << noise_cones[i].angle <<  " " << originCones[i].angle << endl;
        noise_cones[i].cos_angle = originCones[i % num_cones].cos_angle;
//        cout << noise_cones[i].angle << " " << originCones[i].angle << endl;
//        float dca_after = noise_cones[i].get_DCA(PGVector3(0,0,0));
//        cout << dca_before << " " << dca_after << endl<<endl;
    }

    system_matrix_ptr_->clear();
    *system_matrix_ptr_ = 0.0;

    //cout <<"Copying the density estimator . . ." << endl;
    //std::shared_ptr<DensityEstimator> next_density_ptr(density_estimator_ptr_->clone());
    //std::shared_ptr<DensityEstimator> tmp_ptr;

    vector<PGVector3> intercepts;
    float weight = 0.0;
//    size_t j = 0;
    size_t k = 0;

    #pragma omp parallel for private(intercepts, weight, k) num_threads(number_threads_)
    for (size_t i = 0; i < noise_cones.size(); ++i) {

        intercepts = get_intercepts(noise_cones[i], phantom_center_, phantom_length_, intercept_dca_);
        weight = 1.0/static_cast<float>(num_cones * intercepts.size());
        if (i % 1000 == 0) cout << "Cone " << i << " has " << intercepts.size() << " intercepts . . . " << endl;

        for (k = 0; k < intercepts.size(); ++k) {
            system_matrix_ptr_->fill(intercepts[k], 10000.0*weight);
        }
    }
};


void OctaneEMAlgorithm::populate_density_matrix(const vector<OriginCone> &originCones) {

    auto num_cones = originCones.size();
    density_estimator_ptr_->clear();
    *density_estimator_ptr_ = 1.0;

    cout <<"Copying the density estimator . . ." << endl;
    std::shared_ptr<DensityEstimator> next_density_ptr(density_estimator_ptr_->clone());
    std::shared_ptr<DensityEstimator> tmp_ptr;

    float sum = {0.0};
    float weight = 1.0;
    vector<PGVector3> intercepts;
    vector<float> densities;
//    float system_matrix_scaler = 0.0;
    for (size_t m = 0; m < number_iterations_; ++m) {
        cout <<"Starting iteration " << m << "(" << number_threads_ <<" threads) . . ." << endl;

        #pragma omp parallel for private(sum, intercepts, densities, weight) num_threads(number_threads_)
        for (size_t i = 0; i < num_cones; ++i) {
            if (i % 1000 == 0) cout << "Cone " << i << " . . . " << endl;
            intercepts = get_intercepts(originCones[i], phantom_center_, phantom_length_, intercept_dca_);
            densities = density_estimator_ptr_->getDensities(intercepts);

            for (auto d: densities) sum += d;

            for (size_t j = 0; j < densities.size(); ++j) {

//                system_matrix_scaler = (m == 0 ? system_matrix_ptr_->getDensity(intercepts[j]): 1.0); //only use the system matrix in first iteration

                weight = densities[j]/(system_matrix_ptr_->getDensity(intercepts[j]) * sum);

//                weight = densities[j]/sum;
//                system_matrix_scaler = 1.0;
//                cout << system_matrix_scaler << endl;
                if(!(weight > 0.0)){
                    stringstream s;
                    s << "Invalid weight: " << weight << endl;
                    throw runtime_error(s.str());
                }
                next_density_ptr->fill(intercepts[j], weight);
            }
            sum = 0.0;
        }

        tmp_ptr = density_estimator_ptr_;
        density_estimator_ptr_ = next_density_ptr;
        next_density_ptr = tmp_ptr;
        *next_density_ptr = 0.0;
    }
}

void OctaneEMAlgorithm::run(){
    cout<<"OctaneEMAlgorithm::run() . . ."<<endl;
    populate_density_matrix(origin_cones_);
    cout<<"OctaneEMAlgorithm::run() complete . . ."<<endl;
    return;
}


void OctaneEMAlgorithm::setConicSections(const vector<ConicSection> &conic_sections){
    conic_sections_.clear();
    for(size_t i = 0; i<conic_sections.size(); ++i){
        conic_sections_.push_back(conic_sections[i]);
    }
    cout<<"Copied "<<conic_sections.size() <<" conic sections to the OctaneEMAlgorithm  . . ."<<endl;
}
