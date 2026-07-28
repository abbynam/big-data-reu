#ifndef _SOE_ALGORITHM_BUILDER
#define _SOE_ALGORITHM_BUILDER

//C++ standard includes
#include <vector>
#include <ctime>
#include <algorithm>
#include <memory>
#include <cstdlib>

// Custom includes
#include "RunTimeParameters.h"
#include "AbstractImageAlgorithmBuilder.h"
#include "SOEAlgorithm.h"
#include "ASHDensity.h"
#include "DensityMatrix.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{
    
    class SOEAlgorithmBuilder: public AbstractImageAlgorithmBuilder{
    public:
        SOEAlgorithmBuilder(){ };
        virtual ~SOEAlgorithmBuilder(){};
        
        shared_ptr<ImageAlgorithm> build(const vector<ConicSection> &conics, const RunTimeParameters &params) const{
            
            auto phantom_volume_ptr = PhantomVolumeBuilder::build(params);
            
            auto IA = make_shared<SOEAlgorithm>(conics, phantom_volume_ptr);
            IA->setNumberOfIterations(params.get_int("ITERATIONS"));
            
            auto de = get_density_estimator(params, phantom_volume_ptr);
            IA->setDensityEstimator(de);
            
            IA->setVolueConstraints(params.get_double("GAUSS_WIDTH"), params.get_double("OFFSET_X"), params.get_double("OFFSET_Y"));
            IA->setTuningParameters(params.get_int("NUMBER_OF_THREADS"), params.get_double("TEMPERATURE"));
            IA->run();
            
            cout<<"Built SOE Algorithm object . . ."<<endl;
            return IA;
        };        
        
    private: 
        shared_ptr<DensityEstimator> get_density_estimator(const RunTimeParameters &params, shared_ptr<const PhantomVolume> pv_ptr) const{
            
            int density_estimator_type = params.get_int("DENSITY_ESTIMATOR_TYPE");

            if(2 == density_estimator_type){//using averaged shifted histograms
                cout<<"Using Averaged Shifted Histograms for density estimation . . ."<<endl;
                int number_of_shifts = params.get_int("NUMBER_OF_SHIFTS");
                return shared_ptr<DensityEstimator>(
                        new ASHDensity(number_of_shifts,
                                       pv_ptr->x_min, pv_ptr->x_max, params.get_int("X_BINS"),
                                       pv_ptr->y_min, pv_ptr->y_max, params.get_int("Y_BINS"),
                                       pv_ptr->z_min, pv_ptr->z_max, params.get_int("Z_BINS")) );
            }else{
                cout<<"Using a standard 3D histogram for density estimation . . ."<<endl;
                return  shared_ptr<DensityEstimator>(
                    new DensityMatrix(
                            pv_ptr->x_min, pv_ptr->x_max, params.get_int("X_BINS"),
                            pv_ptr->y_min, pv_ptr->y_max, params.get_int("Y_BINS"),
                            pv_ptr->z_min, pv_ptr->z_max, params.get_int("Z_BINS")) );
            }
            
            return nullptr;
        };
    };
}
#endif // _SOE_ALGORITHM_BUILDER