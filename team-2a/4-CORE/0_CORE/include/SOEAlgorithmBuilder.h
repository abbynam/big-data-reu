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
#include "DensityEstimatorFactory.h"

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
            IA->setNumberOfIterations(params.get_int("POINT_ITERATIONS"));

            auto de = DensityEstimatorFactory::createDensityEstimator(params, phantom_volume_ptr);
            IA->setDensityEstimator(de);

            IA->setTuningParameters(params.get_int("NUMBER_OF_THREADS"));
            IA->run();
            
            cout<<"Built SOE Algorithm object . . ."<<endl;
            return IA;
        };        

    };
}
#endif // _SOE_ALGORITHM_BUILDER