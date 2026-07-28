#ifndef PHANTOM_VOLUME_BUILDER_H_
#define PHANTOM_VOLUME_BUILDER_H_
#define _USE_MATH_DEFINES

#include "RandomSingleton.h"
#include "PhantomVolume.h"
#include "RunTimeParameters.h"

namespace prompt_gamma_reconstruction{
    
 /*! \brief Builds phanntom volume 
  * Builds phanntom volume from parameters defined in the run time parameters file. 
  * 
  * @author Dennis Mackin
  * @date Nov. 29, 2015
  * @date Jan. 7, 2015 -- Changed to use min and max rather than length
  */  
class PhantomVolumeBuilder{
    public:
        static shared_ptr<const PhantomVolume> build(const pg_tools::RunTimeParameters &params){
            
            auto phantom_ptr = shared_ptr<PhantomVolume>(new PhantomVolume());

            phantom_ptr->x_bins_ = params.get_int("X_BINS");
            if(phantom_ptr->x_bins_ % 2 == 0) phantom_ptr->x_bins_++; //use odd number of bins to avoid assymmetry around 0

            phantom_ptr->x_max = params.get_double("X_MAX");
            phantom_ptr->x_min = params.get_double("X_MIN");

            phantom_ptr->y_bins_ = params.get_int("Y_BINS");
            if(phantom_ptr->y_bins_ % 2 == 0) phantom_ptr->y_bins_++; //use odd number of bins to avoid asymmetry around 0

            phantom_ptr->y_max = params.get_double("Y_MAX");
            phantom_ptr->y_min = params.get_double("Y_MIN");

            phantom_ptr->z_bins_ = params.get_int("Z_BINS");
            if(phantom_ptr->z_bins_ % 2 == 0) phantom_ptr->z_bins_++; //use odd number of bins to avoid assymmetry around 0

            phantom_ptr->z_max = params.get_double("Z_MAX");
            phantom_ptr->z_min = params.get_double("Z_MIN");

            return phantom_ptr;            
        };

    static shared_ptr<const PhantomVolume> buildOctanePhantom(float x, float y, float z, float length, size_t bins){

        auto phantom_ptr = shared_ptr<PhantomVolume>(new PhantomVolume());

        phantom_ptr->x_bins_ = bins;
//        if(phantom_ptr->x_bins_ % 2 == 0) phantom_ptr->x_bins_++; //use odd number of bins to avoid assymmetry around 0

        phantom_ptr->x_max = x + 0.5*length;
        phantom_ptr->x_min = x - 0.5*length;

        phantom_ptr->y_bins_ = bins;
//        if(phantom_ptr->y_bins_ % 2 == 0) phantom_ptr->y_bins_++; //use odd number of bins to avoid asymmetry around 0

        phantom_ptr->y_max = y + 0.5*length;
        phantom_ptr->y_min = y - 0.5*length;

        phantom_ptr->z_bins_ = bins;
//        if(phantom_ptr->z_bins_ % 2 == 0) phantom_ptr->z_bins_++; //use odd number of bins to avoid assymmetry around 0

        phantom_ptr->z_max = z + 0.5*length;
        phantom_ptr->z_min = z - 0.5*length;

        return phantom_ptr;
    };

};

}//end of namespace
#endif //PHANTOM_VOLUME_BUILDER_
