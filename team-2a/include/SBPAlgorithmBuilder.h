#ifndef _SBP_ALGORITHM_BUILDER
#define _SBP_ALGORITHM_BUILDER

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
#include "SBPAlgorithm.h"
#include "PhantomVolume.h"
#include "PhantomVolumeBuilder.h"
#include "ASHDensity.h"
#include "DensityMatrix.h"
#include "DensityEstimatorFactory.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{
        class SBPAlgorithmBuilder: public AbstractImageAlgorithmBuilder{
    public:
        SBPAlgorithmBuilder(){cout<<"Ready to build SBP algorithm object . . ."<<endl; };
        ~SBPAlgorithmBuilder(){};


        shared_ptr<ImageAlgorithm> build(const vector<ConicSection> &conics, const RunTimeParameters &params) const {

            cout<<"Building SBP algorithm object . . ."<<endl;
            auto algo = make_shared<SBPAlgorithm>(conics);

            algo->setBandwidth(params.get_float("KERNEL_BANDWIDTH"));

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
#endif // _SBP_ALGORITHM_BUILDER