#include "ReconstructorTemplate.h"
#include "AbstractResultsFormatter.h"

namespace prompt_gamma_reconstruction{
    
    void  ReconstructorTemplate::run(){
        start_time_ = std::time(nullptr);
        results_formatter_ptr_->saveOutput();
    }; //Kicks off the density estimation
    
    void ReconstructorTemplate::saveResults(){};

};