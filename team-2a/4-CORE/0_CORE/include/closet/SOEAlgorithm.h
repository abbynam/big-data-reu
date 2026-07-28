#ifndef _SOE_ALGORITHM
#define _SOE_ALGORITHM

//C++ standard includes
#include <vector>
#include <ctime>
#include <algorithm>
#include <memory>
#include <cstdlib>

// Custom includes
#include "ConicSection.h"
#include "RunTimeParameters.h"
#include "DensityEstimator.h"
#include "ImageAlgorithm.h"

using namespace std;
namespace prompt_gamma_reconstruction{
    
    class SOEAlgorithm: public ImageAlgorithm{
    public:
        SOEAlgorithm(const vector<ConicSection> &conic_sections, shared_ptr<const PhantomVolume> phantom_volume): phantom_volume_(phantom_volume),
                gaussian_width_(10000.0), offset_x_(0.0), offset_y_(0.0), number_tries_for_random_(1000), number_of_threads_(1), temperature_(1.0)
        { 
            setConicSections(conic_sections);
        };
        ~SOEAlgorithm() { };
        
        
        //Get 2D image data 
        Image2D getImagePlane(int dimension, float depth) const;
        
        //Get 3D image data
        Image3D getImageVolume(int dimension) const;
        
        //Flexible function that can return the data as a string in any format
        string getDataAsString() const;
        
        //SETTERS
        void setNumberOfIterations(const int iterations){number_of_iterations_ = iterations;};
        void setConicSections(const vector<ConicSection> &conic_sections);
        void setDensityEstimator(shared_ptr<DensityEstimator> de){
            density_estimator_ptr_ = de;
            populate_density_matrix_(conic_sections_, *density_estimator_ptr_);
        }
        
        void setVolueConstraints(float gaussian_width, float offset_x, float offset_y){
            gaussian_width_ = gaussian_width;
            offset_x_ = offset_x;
            offset_y_ = offset_y;
        }
        
        void setTuningParameters(int threads, double temperature){
            number_of_threads_ = threads;
            temperature_ = temperature;
        }
        
        void setNumberOfTriesForRandom(int tries){ number_tries_for_random_ = tries; };
        void run();
        
    private:
        void populate_density_matrix_(const vector<ConicSection> &conicSections, DensityEstimator &density_estimator);

        //Image reconstruction steps 
        void calculate_image_();
        int number_of_iterations_;

        //PROPERTIES
        vector<ConicSection> conic_sections_; ///vector to store pointers to the parabolas and ellipses
        shared_ptr<DensityEstimator> density_estimator_ptr_; 
        shared_ptr<const PhantomVolume> phantom_volume_;
        std::time_t start_time_;
        double gaussian_width_;
        double offset_x_;
        double offset_y_;
        int number_tries_for_random_;
        int number_of_threads_;
        double temperature_;
    };
};
#endif // _SOE_ALGORITHM
