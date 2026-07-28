/* ****************************************************************************
 *  PointAlgorithm PointAlgorithm -
 *
 * \section intro_sec Overview
 *
 * PointAlgorithm is the improved version of SOEAlgorthm. It makes no
 * attempt to be loyal to SOE algorithm details outlined in the Sitek et al or Mackin et al papers.
 *
 *
 * @author Dennis Mackin
 * @date December 30, 2015
 */

// C++ Includes
#include <sstream>

// Custom Includes
#include "SOEAlgorithm.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

Image2D SOEAlgorithm::getImagePlane(size_t dimension, float depth) const{
    Image2D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});

    if(0 == dimension){
        cout<<"Producing image for yz plane, z = "<< depth<< ". . ."<<endl;
    }else if(1 == dimension){
        cout<<"Producing image for xz plane, z = "<< depth<< ". . ."<<endl;
    }else if(2 == dimension){
        cout<<"Producing image for xy plane, z = "<< depth<< ". . ."<<endl;
    }else{
        stringstream err_msg;
        err_msg <<"Invalid dimension " << dimension << "for SOEAlgorithm:getImagePlane." << endl;
        throw runtime_error(err_msg.str());
    }

    return I;
}

Image3D SOEAlgorithm::getImageVolume(size_t dimension) const{
    Image3D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});

    return I;
}

string SOEAlgorithm::getDataAsString() const{

    auto dose_string = density_estimator_ptr_->get3DDose();
    return dose_string;
}

string SOEAlgorithm::getDataAsString(size_t nx, size_t ny, size_t nz) const{

    auto dose_string = density_estimator_ptr_->get3DDose(nx, ny, nz);
    return dose_string;
}

void SOEAlgorithm::setConicSections(const vector<ConicSection> &conic_sections){
    conic_sections_ = conic_sections;
//    set_inverse_square_params_();
}
//
//void SOEAlgorithm::set_inverse_square_params_() {
//    for(auto iter = conic_sections_.begin(); iter != conic_sections_.end(); ++iter){
//        iter->setInverseSquareParam(this->inverse_square_parameter_);
//    }
//};

void SOEAlgorithm::run(){
    // Variable Declaration
    size_t changeCount = 0;

    PGVector3 oldEventPos, randomPosition;
    double oldDensity, newDensity;

    // Loop through iterations
    auto num_cones = conic_sections_.size();
    printf("--- Number of MCMC iterations: %zu ---\n", number_of_iterations_);
    printf("--- Number of events: %lu ---\n", num_cones);


    size_t i = 0;
    for (; i < number_of_iterations_; i++) {
        changeCount = 0;

        long long index = i;
#pragma omp parallel for reduction(+:changeCount), \
                private(index, oldEventPos, randomPosition, oldDensity, newDensity) num_threads(number_of_threads_)
        for(size_t jEventNum=0; jEventNum< num_cones; ++jEventNum){

            oldEventPos = conic_sections_[jEventNum].getLikelyOrigin();
            oldDensity = density_estimator_ptr_->getDensity(oldEventPos);

            long number_tries;
            number_tries = conic_sections_[jEventNum].getRandomPointInPhantom(randomPosition, number_tries_for_random_);

            if( -1 == number_tries) continue; //did not find random point so continue to next point

            newDensity = density_estimator_ptr_->getDensity(randomPosition);

            index = static_cast<long long>(num_cones)*i + jEventNum;
            float rand = RandomSingleton::Instance()->getRandIndex(index);

            //Acquire the mutex and make the change if the new density
            // is greater than the old density times a random raised to
            // power temperature.
            if(newDensity >= rand*oldDensity){
                conic_sections_[jEventNum].setLikelyOrigin(randomPosition);
                #pragma omp critical
                {
                    density_estimator_ptr_->updateMatrix(oldEventPos,randomPosition, 1.0);
                }
                ++changeCount;
            }
        }//end for

        //Give status update to the terminal
        if ( ( i < 100 &&  (i+1) % 10 == 0) ||
             ( i < 1000 &&  (i+1) % 100 == 0) ||
             ( i < 10000 &&  (i+1) % 1000 == 0) ||
             ( i < 25000 &&  (i+1) % 5000 == 0) ||
             ( i < 100000 &&  (i+1) % 10000 == 0)  ){
            printf("Iteration: %zu, time %ld, Number of Position Changes: %zu, ratio: %.3f",
                   i+1, get_running_time(), changeCount, static_cast<double>(changeCount)/static_cast<double>(num_cones) );
            cout<<endl;//print it out now!

        }
    }
    printf("--- Total Iterations: %zu, time %ld, Number of Position Changes: %zu, %lu, ratio: %.3f\n",
           number_of_iterations_, get_running_time(), changeCount, num_cones, (double)changeCount/(double)(num_cones));
}


void SOEAlgorithm::populate_density_matrix_(const vector<ConicSection> &conicSections, DensityEstimator &density_estimator){
    density_estimator.clear();
    density_estimator = 1.0;

    //vector<ConicSection>::iterator conic_iter = conicSections.begin();
    auto conic_iter = conicSections.begin();
    for(/* */; conic_iter != conicSections.end(); ++conic_iter){
        auto origin = conic_iter->getLikelyOrigin();
        auto weight = conic_iter->getWeight();
        density_estimator.fill(origin, weight);
    }
}

string SOEAlgorithm::getConicInformationAsString() const{
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
        ss  << delimiter << point.x << delimiter << point.y << delimiter << point.z << endl;
        //ss << "\n";
    }

    return ss.str();
}