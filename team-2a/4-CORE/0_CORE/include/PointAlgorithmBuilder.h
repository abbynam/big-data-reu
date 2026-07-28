#ifndef _POINT_ALGORITHM_BUILDER
#define _POINT_ALGORITHM_BUILDER

//C++ standard includes
#include <vector>
#include <ctime>
#include <algorithm>
#include <memory>
#include <cstdlib>

// Custom includes
#include "RunTimeParameters.h"
#include "AbstractImageAlgorithmBuilder.h"
#include "PointAlgorithm.h"
#include "ASHDensity.h"
#include "DensityMatrix.h"
#include "DensityEstimatorFactory.h"
#include "PhantomVolumeBuilder.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{
    
    class PointAlgorithmBuilder: public AbstractImageAlgorithmBuilder{
    public:
        PointAlgorithmBuilder(){ };
        virtual ~PointAlgorithmBuilder(){};
        
        shared_ptr<ImageAlgorithm> build(const vector<ConicSection> &conics, const RunTimeParameters &params) const{
            
            auto phantom_volume_ptr = PhantomVolumeBuilder::build(params);

            auto is_param = params.get_float("INVERSE_SQUARE_PARAM");
            auto IA = make_shared<PointAlgorithm>(conics, phantom_volume_ptr, is_param);
            IA->setNumberOfIterations(params.get_int("POINT_ITERATIONS"));
            IA->setEventMultiplier(params.get_int("EVENT_MULTIPLIER"));
            auto dca_x = params.get_float("DCA_CENTER_X");
            auto dca_y = params.get_float("DCA_CENTER_Y");
            auto dca_z = params.get_float("DCA_CENTER_Z");
            IA->setDCACenter(PGVector3(dca_x, dca_y, dca_z));


            auto de = DensityEstimatorFactory::createDensityEstimator(params, phantom_volume_ptr);
            IA->setDensityEstimator(de);

            IA->setTuningParameters(params.get_int("NUMBER_OF_THREADS"), params.get_double("TEMPERATURE"));
            IA->run();
            
            cout<<"Built SOE Algorithm object . . ."<<endl;
            return IA;
        };        
        

    };
}
#endif // _POINT_ALGORITHM_BUILDER