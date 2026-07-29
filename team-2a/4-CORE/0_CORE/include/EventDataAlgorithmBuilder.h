#ifndef _EVENT_DATA_ALGORITHM_BUILDER
#define _EVENT_DATA_ALGORITHM_BUILDER

//C++ standard includes
#include <vector>
#include <ctime>
#include <algorithm>
#include <memory>
#include <cstdlib>

// Custom includes
#include "RunTimeParameters.h"
#include "ImageAlgorithm.h"
#include "EventDataAlgorithm.h"
#include "PhantomVolume.h"
#include "PhantomVolumeBuilder.h"
#include "DensityEstimatorFactory.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{
        class EventDataAlgorithmBuilder: public AbstractImageAlgorithmBuilder{
    public:
        EventDataAlgorithmBuilder(){ cout<<"Ready to build EventData algorithm object . . ."<<endl; };
        ~EventDataAlgorithmBuilder(){};

        shared_ptr<ImageAlgorithm> build(const vector<ConicSection> &conics, const RunTimeParameters &params) const {

            cout<<"Building EventData algorithm object . . ."<<endl;
            auto phantom_volume_ptr = PhantomVolumeBuilder::build(params);

            auto EDA = make_shared<EventDataAlgorithm>(conics, phantom_volume_ptr);
            auto dca_x = params.get_float("DCA_CENTER_X");
            auto dca_y = params.get_float("DCA_CENTER_Y");
            auto dca_z = params.get_float("DCA_CENTER_Z");
            EDA->setDCACenter(PGVector3(dca_x, dca_y, dca_z));

            cout<<"Preparing the density estimator . . ."<<endl;
            auto de = DensityEstimatorFactory::createDensityEstimator(params, phantom_volume_ptr);
            EDA->setDensityEstimator(de);

            cout<<"Ready to run . . ."<<endl;
            EDA->run();
            
            return EDA;
        };        
    };
};
#endif // _EVENT_DATA_ALGORITHM_BUILDER