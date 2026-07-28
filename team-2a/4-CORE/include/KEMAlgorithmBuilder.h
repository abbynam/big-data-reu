#ifndef _KEM_ALGORITHM_BUILDER
#define _KEM_ALGORITHM_BUILDER

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
#include "KEMAlgorithm.h"
#include "PhantomVolume.h"
#include "PhantomVolumeBuilder.h"
#include "ASHDensity.h"
#include "DensityMatrix.h"
#include "DensityEstimatorFactory.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{
        class KEMAlgorithmBuilder: public AbstractImageAlgorithmBuilder{
    public:
        KEMAlgorithmBuilder(){cout<<"Ready to build Octane algorithm object . . ."<<endl; };
        ~KEMAlgorithmBuilder(){};

        shared_ptr<ImageAlgorithm> build(const vector<ConicSection> &conics, const RunTimeParameters &params) const {
            cout<<"Building KEM algorithm object . . ."<<endl;
            auto algo = make_shared<KEMAlgorithm>(conics);

            algo->setBandwidth(params.get_float("KERNEL_BANDWIDTH"));
            //algo->setNumberIterations(params.get_int("OCTANE_ITERATIONS"));
            algo->setNumberThreads(params.get_float("NUMBER_OF_THREADS"));
            algo->setSystemMatrixScalar(params.get_float("SYSTEM_MATRIX_SCALAR"));

            auto dca_x = params.get_float("DCA_CENTER_X");
            auto dca_y = params.get_float("DCA_CENTER_Y");
            auto dca_z = params.get_float("DCA_CENTER_Z");
            algo->setDCACenter(PGVector3(dca_x, dca_y, dca_z));
            algo->buildPhantomVolume(params);

            cout<<"Ready to run . . ."<<endl;
            algo->run();

            return algo;
        };        
    };
};
#endif // _KEM_ALGORITHM_BUILDER