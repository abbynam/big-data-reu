#ifndef _DCA_LINE_ALGORITHM_BUILDER
#define _DCA_LINE_ALGORITHM_BUILDER

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
#include "DCALineAlgorithm.h"
#include "PhantomVolume.h"
#include "PhantomVolumeBuilder.h"
#include "ASHDensity.h"
#include "DensityMatrix.h"
#include "DensityEstimatorFactory.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{
        class DCALineAlgorithmBuilder: public AbstractImageAlgorithmBuilder{
    public:
        DCALineAlgorithmBuilder(){cout<<"Ready to build Octane algorithm object . . ."<<endl; };
        ~DCALineAlgorithmBuilder(){};

        shared_ptr<ImageAlgorithm> build(const vector<ConicSection> &conics, const RunTimeParameters &params) const {

            cout<<"Building DCALineAlgorithm object . . ."<<endl;

            auto phantom_volume_ptr = PhantomVolumeBuilder::build(params);

            auto algo = make_shared<DCALineAlgorithm>(conics, phantom_volume_ptr);
            auto p = params.get_csv_values("BEAM_LINE_POINT1");
            PGVector3 p1(p[0], p[1], p[2]);
            assert(p.size() == 3);
            p = params.get_csv_values("BEAM_LINE_POINT2");
            PGVector3 p2(p[0], p[1], p[2]);
            assert(p.size() == 3);
            algo->setPoints(p1, p2);
            algo->setNumberOfThreads(params.get_int("NUMBER_OF_THREADS"));


            cout<<"Preparing the density estimator . . ."<<endl;
            auto de = DensityEstimatorFactory::createDensityEstimator(params, phantom_volume_ptr);
            algo->setDensityEstimator(de);

            algo->run();
            cout<<"Built DCALineAlgorithm object . . ."<<endl;

            return algo;
        };        
    };
};
#endif // _DCA_LINE_ALGORITHM_BUILDER