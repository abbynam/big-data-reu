#ifndef _OCTANE_EM_ALGORITHM_BUILDER
#define _OCTANE_EM_ALGORITHM_BUILDER

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
#include "OctaneEMAlgorithm.h"
#include "PhantomVolume.h"
#include "PhantomVolumeBuilder.h"
#include "ASHDensity.h"
#include "DensityMatrix.h"
#include "DensityEstimatorFactory.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{
        class OctaneEMAlgorithmBuilder: public AbstractImageAlgorithmBuilder{
    public:
        OctaneEMAlgorithmBuilder(){cout<<"Ready to build Octane algorithm object . . ."<<endl; };
        ~OctaneEMAlgorithmBuilder(){};

        shared_ptr<ImageAlgorithm> build(const vector<ConicSection> &conics, const RunTimeParameters &params) const {

            cout<<"Building Octane algorithm object . . ."<<endl;
            float bin_width = params.get_float("BIN_WIDTH");
            float min_length = params.get_float("PHANTOM_LENGTH");
            assert(min_length >= 2.0);
            float length = bin_width;
            while(length < min_length) length *= 2.0;
            cout<<"Phantom length: " << length <<" . . . "<<endl;
            size_t bins = length/bin_width;
            cout<<"Phantom length: " << length <<" bins: " << bins << ". . . "<<endl;
            auto phantom_volume_ptr = PhantomVolumeBuilder::buildOctanePhantom(params.get_float("PHANTOM_CENTER_X"),
                                                                               params.get_float("PHANTOM_CENTER_Y"),
                                                                               params.get_float("PHANTOM_CENTER_Z"),
                                                                               length,
                                                                               bins);
//            auto phantom_volume_ptr = PhantomVolumeBuilder::build(params);

            auto OA = make_shared<OctaneEMAlgorithm>(conics, phantom_volume_ptr);
            OA->setInterceptDCA(bin_width);
//            OA->setSystemMatrix(system_matrix_ptr);
            //OA->setInterceptDCA(params.get_float("INTERCEPT_DCA"));
            OA->setPhantomLength(length);
            OA->setPhantomCenter(params.get_float("PHANTOM_CENTER_X"), params.get_float("PHANTOM_CENTER_Y"), params.get_float("PHANTOM_CENTER_Z"));
            //OA->setSAD(params.get_float("SOURCE_AXIS_DISTANCE"));
            OA->setInverseSquareParameter(params.get_float("INVERSE_SQUARE_PARAM"));
            OA->setNumberIterations(params.get_float("OCTANE_ITERATIONS"));
            OA->setNumberThreads(params.get_float("NUMBER_OF_THREADS"));

            auto dca_x = params.get_float("DCA_CENTER_X");
            auto dca_y = params.get_float("DCA_CENTER_Y");
            auto dca_z = params.get_float("DCA_CENTER_Z");
            OA->setDCACenter(PGVector3(dca_x, dca_y, dca_z));


            cout<<"Preparing the density estimator . . ."<<endl;
            auto de = DensityEstimatorFactory::createDensityEstimator(params, phantom_volume_ptr);
            OA->setDensityEstimator(de);
            OA->setSystemMatrixEstimator(de);

            //auto system_matrix_ptr = DensityEstimatorFactory::createDensityMatrix(phantom_volume_ptr, 80, 80, 80);
//            auto system_matrix_ptr = DensityEstimatorFactory::createDensityEstimator(params, phantom_volume_ptr);
//            OA->setSystemMatrix(system_matrix_ptr);

            cout<<"Ready to run . . ."<<endl;
            OA->run();

            return OA;
        };        
    };
};
#endif // _OCTANE_EM_ALGORITHM_BUILDER