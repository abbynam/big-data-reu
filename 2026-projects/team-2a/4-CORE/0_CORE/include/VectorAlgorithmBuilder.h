#ifndef _VECTOR_ALGORITHM_BUILDER
#define _VECTOR_ALGORITHM_BUILDER

//C++ standard includes
#include <vector>
#include <ctime>
#include <algorithm>
#include <memory>
#include <cstdlib>

// Custom includes
#include "AbstractImageAlgorithmBuilder.h"
#include "RunTimeParameters.h"
#include "ImageAlgorithm.h"
#include "VectorAlgorithm.h"
#include "PhantomVolume.h"
#include "PhantomVolumeBuilder.h"
#include "ASHDensity.h"
#include "DensityMatrix.h"
#include "DensityEstimatorFactory.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{
        class VectorAlgorithmBuilder: public AbstractImageAlgorithmBuilder{
    public:
        VectorAlgorithmBuilder(){cout<<"Ready to build Vector algorithm object . . ."<<endl; };
        ~VectorAlgorithmBuilder(){};


        shared_ptr<ImageAlgorithm> build(const vector<ConicSection> &conics, const RunTimeParameters &params) const {

            cout<<"Building Vector algorithm object . . ."<<endl;
            auto phantom_volume_ptr = PhantomVolumeBuilder::buildOctanePhantom(params.get_float("PHANTOM_CENTER_X"),
                                                                               params.get_float("PHANTOM_CENTER_Y"),
                                                                               params.get_float("PHANTOM_CENTER_Z"),
                                                                               params.get_int("PHANTOM_LENGTH"),
                                                                               params.get_int("PHANTOM_BINS"));
//            auto phantom_volume_ptr = PhantomVolumeBuilder::build(params);

            auto OA = make_shared<VectorAlgorithm>(conics, phantom_volume_ptr);
            OA->setInterceptDCA(params.get_float("INTERCEPT_DCA"));
            OA->setPhantomLength(params.get_float("PHANTOM_LENGTH"));
            OA->setPhantomCenter(params.get_float("PHANTOM_CENTER_X"), params.get_float("PHANTOM_CENTER_Y"), params.get_float("PHANTOM_CENTER_Z"));
            OA->setSAD(params.get_float("SOURCE_AXIS_DISTANCE"));
            OA->setConeLengthCorrectionFactor(params.get_float("CONE_LENGTH_CORRECTION"));

            auto dca_x = params.get_float("DCA_CENTER_X");
            auto dca_y = params.get_float("DCA_CENTER_Y");
            auto dca_z = params.get_float("DCA_CENTER_Z");
            OA->setDCACenter(PGVector3(dca_x, dca_y, dca_z));


            cout<<"Preparing the density estimator . . ."<<endl;
            auto de = DensityEstimatorFactory::createDensityEstimator(params, phantom_volume_ptr);
            OA->setDensityEstimator(de);

            cout<<"Ready to run . . ."<<endl;
            OA->run();

            return OA;
        };        
    };
};
#endif // _VECTOR_ALGORITHM_BUILDER