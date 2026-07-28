/* ****************************************************************************
 *  DCALIneAlgorithm -
 *
 * \section intro_sec Overview
 *
 * DCA line algorithm produces an image by plotting the position on each cone that is closest
 * to a line defined by two points specified in the paramaters file.
 * 
 * @author Dennis Mackin
 * @date August 08, 2016
 */

// C++ Includes
#include <sstream>

// Custom Includes
#include "DCALineAlgorithm.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

Image2D DCALineAlgorithm::getImagePlane(size_t dimension, float depth) const{
    Image2D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    
    if(0 == dimension){
        cout<<"Producing image for yz plane, z = "<< depth<< ". . ."<<endl;
    }else if(1 == dimension){
        cout<<"Producing image for xz plane, z = "<< depth<< ". . ."<<endl;        
    }else if(2 == dimension){
        cout<<"Producing image for xy plane, z = "<< depth<< ". . ."<<endl;        
    }else{
        stringstream err_msg;
        err_msg <<"Invalid dimension " << dimension << "for DCALineAlgorithm:getImagePlane." << endl;
        throw runtime_error(err_msg.str());
    }
    
    return I;
}

Image3D DCALineAlgorithm::getImageVolume(size_t dimension) const{
    Image3D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    
    return I;
}

string DCALineAlgorithm::getDataAsString() const{
    auto dose_string = density_estimator_ptr_->get3DDose();
    return dose_string;
}


string DCALineAlgorithm::get_event_record_(long event_num) const{
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
    auto pca = conic_sections_[event_num].getPointOfClosestApproach(p_dca_);
    ss << delimiter << dca << delimiter << p_dca_.x << delimiter << p_dca_.y << delimiter << p_dca_.z
    << delimiter << pca.x << delimiter << pca.y << delimiter << pca.z;

    auto point = conic_sections_[event_num].getLikelyOrigin();
    ss  << delimiter << point.x << delimiter << point.y << delimiter << point.z << endl;

    return ss.str();
}

string DCALineAlgorithm::getConicInformationAsString() const{

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


void DCALineAlgorithm::populate_density_matrix(const vector<ConicSection> &conics){

    size_t num_cones = conics.size();

    //Create extra space for cones that intercept line and have two PCA's
    vector<PGVector3> positions(num_cones*2);
    vector<float> weights(num_cones*2);
    vector<float> DCAs(num_cones*2,0);

    #pragma omp parallel for
    for(size_t i=0; i<num_cones; ++i){
        if( i % 1000 == 0) cout<<"Cone "<<i<<" . . . "<<endl;
        auto results = conics[i].getDistanceToLine(p1_, p2_);
        positions[i] = results.second;
        weights[i] = 1.0; //@todo Find where cone intercepts line twice and add both positions with weight 1/2.
        DCAs[i] = results.first;
    }
    for(size_t i=0; i<num_cones; ++i){
//        cout<<positions[i].print()<<endl;
        try{
            if(phantom_volume_ptr_->is_in_volume(positions[i])){
                this->density_estimator_ptr_->fill(positions[i], 1.0);
            }
        } catch(std::exception e){
            cout<< e.what()<< endl;
        }
    }
}

void DCALineAlgorithm::run(){
    cout<<"DCALineAlgorithm::run() . . ."<<endl;
    populate_density_matrix(conic_sections_);
    cout<<"DCALineAlgorithm::run() complete . . ."<<endl;
    return;
}


void DCALineAlgorithm::setConicSections(const vector<ConicSection> &conic_sections){
    conic_sections_ = conic_sections;
    cout<<"Copied "<<conic_sections.size() <<" conic sections to the DCALineAlgorithm  . . ."<<endl;
}
