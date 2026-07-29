#ifndef _RECONSTRUCTOR_TEMPLATE
#define _RECONSTRUCTOR_TEMPLATE
//Standard includes
#include <vector>
#include <ctime>
#include <algorithm>
#include <memory>

//Custom includes
#include "ConicSection.h"
#include "RunTimeParameters.h"
#include "DensityEstimator.h"
#include "ImageAlgorithm.h"
#include "AbstractResultsFormatter.h"

using namespace std;
namespace prompt_gamma_reconstruction{
    
    
class ReconstructorTemplate{
    public:
        ReconstructorTemplate(){
            /* DO NOTHING*/
        };
        ~ReconstructorTemplate(){};
        void run();
        
        //SETTERS
        
        void setImageAlgorithm(shared_ptr<ImageAlgorithm> &IA){
            image_algorithm_ptr_ = IA;
        };
        
        void setResultsFormatter(shared_ptr<AbstractResultsFormatter> &results_formatter_ptr){
            results_formatter_ptr_ = results_formatter_ptr;
        };        
        
    private:
        std::time_t get_running_time(){ return std::time(nullptr) - start_time_;};

        //Image reconstruction steps 
        void generateImages();
        void saveResults();
        

        //PROPERTIES
        //pg_tools::RunTimeParameters params_; //Parameters from config file
        vector<ConicSection> conic_sections_; ///vector to store pointers to the parabolas and ellipses
        shared_ptr<DensityEstimator> density_estimator_ptr_;  
        shared_ptr<ImageAlgorithm> image_algorithm_ptr_;
        shared_ptr<AbstractResultsFormatter> results_formatter_ptr_;  
        std::time_t start_time_;
    };
};
#endif // _RECONSTRUCTOR_TEMPLATE