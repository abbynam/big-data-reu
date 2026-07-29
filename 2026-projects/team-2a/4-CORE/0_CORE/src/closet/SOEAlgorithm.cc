/* ****************************************************************************
 *  SOEAlgorithm SOEAlgorithm -
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
#include "SOEAlgorithm.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

Image2D SOEAlgorithm::getImagePlane(int dimension, float depth) const{
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

Image3D SOEAlgorithm::getImageVolume(int dimension) const{
    Image3D I(vector<float>{-3.0,-3.0,-3.0}, vector<float>{0.0,0.0,0.0}, vector<size_t>{0,0,0});
    
    return I;
}

string SOEAlgorithm::getDataAsString() const{
    
    auto dose_string = density_estimator_ptr_->get3DDose();
    return dose_string;
}

void SOEAlgorithm::setConicSections(const vector<ConicSection> &conic_sections){
    conic_sections_ = conic_sections;    
}  
        
void SOEAlgorithm::run(){
// Variable Declaration
    int changeCount = 0;

    PGVector3 oldEventPos, randomPosition;
    double oldDensity, newDensity;

    // Loop through iterations
    auto num_cones = conic_sections_.size();
    printf("--- Number of MCMC iterations: %d ---\n", number_of_iterations_);
    printf("--- Number of events: %lu ---\n", num_cones);

//    int NUM_THREADS = run_time_parameters_.get_int("NUMBER_OF_THREADS");
    
    double old_prior = 1.0; //Gaussian likelihood weight of current representative point
    double new_prior = 1.0; //Gaussian likelihood weight of random test point
    double event_weight = 0.0; //Better (logistic regression) events weight more
    double C1 = -0.5/(gaussian_width_*gaussian_width_);

    int i = 0;
    for (; i < number_of_iterations_; i++) {
        changeCount = 0;

        long long index = i;
        #pragma omp parallel for reduction(+:changeCount), \
                private(index, oldEventPos, randomPosition, oldDensity, newDensity, event_weight, old_prior, new_prior ) num_threads(number_of_threads_)
        for(size_t jEventNum=0; jEventNum< num_cones; ++jEventNum){

            oldEventPos = conic_sections_[jEventNum].getLikelyOrigin();
            old_prior = exp(C1*( (offset_x_ - oldEventPos.x)*(offset_x_-oldEventPos.x)+(offset_y_ - oldEventPos.y)*(offset_y_ - oldEventPos.y)));

            event_weight = conic_sections_[jEventNum].getWeight();

            oldDensity = (density_estimator_ptr_->getDensity(oldEventPos) - event_weight)*old_prior;

            int number_tries;
            number_tries = conic_sections_[jEventNum].getRandomPointInPhantom(randomPosition, number_tries_for_random_);

            if( -1 == number_tries) continue; //did not find random point so continue to next point

            new_prior = exp(C1*((offset_x_ - randomPosition.x)*(offset_x_ - randomPosition.x)+(offset_y_ - randomPosition.y)*(offset_y_ - randomPosition.y)));

            newDensity = density_estimator_ptr_->getDensity(randomPosition)*new_prior;

            index = static_cast<long long>(num_cones)*i + jEventNum;
            float rand = RandomSingleton::Instance()->getRandIndex(index);

            //Acquire the mutex and make the change if the new density
            // is greater than the old density times a random raised to 
            // power temperature.
            if( newDensity >= pow(rand, temperature_)*( oldDensity )){
                conic_sections_[jEventNum].setLikelyOrigin(randomPosition);
                #pragma omp critical
                {
                    density_estimator_ptr_->updateMatrix(oldEventPos,randomPosition,event_weight);
                }
                ++changeCount;
            }
        }//end for

        //Give status update to the terminal
        if ( ( i < 1000 &&  (i+1) % 100 == 0) ||
             ( i < 10000 &&  (i+1) % 1000 == 0) ||
             ( i < 25000 &&  (i+1) % 5000 == 0) ||
             ( i < 100000 &&  (i+1) % 10000 == 0)  ){
            printf("Iteration: %d, time %ld, Number of Position Changes: %d, ratio: %.3f",
                    i+1, get_running_time(), changeCount, static_cast<double>(changeCount)/static_cast<double>(num_cones) );
            cout<<endl;//print it out now!

            //Save data for analysis
//            save_density_info_(width_num*iterations + i+1, true);
        }
    }
//    save_density_info_(width_num*iterations+iterations, true);
    printf("--- Total Iterations: %d, time %ld, Number of Position Changes: %d, %lu, ratio: %.3f\n",
                number_of_iterations_, get_running_time(), changeCount, num_cones, (double)changeCount/(double)(num_cones));
}


void SOEAlgorithm::populate_density_matrix_(const vector<ConicSection> &conicSections, DensityEstimator &density_estimator){
  density_estimator.clear();
  density_estimator = 0.0;
  
  //vector<ConicSection>::iterator conic_iter = conicSections.begin();
  auto conic_iter = conicSections.begin();
  for(/* */; conic_iter != conicSections.end(); ++conic_iter){
      auto origin = conic_iter->getLikelyOrigin();
      auto weight = conic_iter->getWeight();
	  density_estimator.fill(origin, weight);
  }
}