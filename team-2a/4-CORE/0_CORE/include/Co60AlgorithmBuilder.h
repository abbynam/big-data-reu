#ifndef _CO60_ALGORITHM_BUILDER
#define _CO60_ALGORITHM_BUILDER

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
#include "Co60Algorithm.h"
#include "PhantomVolume.h"
#include "PhantomVolumeBuilder.h"
#include "ASHDensity.h"
#include "DensityMatrix.h"
#include "DensityEstimatorFactory.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{
        class Co60AlgorithmBuilder: public AbstractImageAlgorithmBuilder{
    public:
        Co60AlgorithmBuilder(){cout<<"Ready to build Octane algorithm object . . ."<<endl; };
        ~Co60AlgorithmBuilder(){};

        shared_ptr<ImageAlgorithm> build(const vector<ConicSection> &conics, const RunTimeParameters &params) const {

            cout<<"Building Co60Algorithm object . . ."<<endl;

            auto phantom_volume_ptr = PhantomVolumeBuilder::build(params);

            auto algo = make_shared<Co60Algorithm>(conics, phantom_volume_ptr);

            auto p = params.get_csv_values("ORIGIN");
            assert(p.size() == 3);
            PGVector3 origin(p[0], p[1], p[2]);
            algo->setOrigin(origin);

            algo->setMinScatterEnergy(params.get_float("MIN_ENERGY_SCATTER"));
            algo->setMinEventEnergy(params.get_float("MIN_ENERGY_EVENT"));


            algo->setNumberOfThreads(params.get_int("NUMBER_OF_THREADS"));

            OctaneAlgorithmBuilder b;
            auto ia = b.build(conics, params);
            algo->setImageAlgorithm(ia);

            return algo;
        };        
    };
};
#endif // _CO60_ALGORITHM_BUILDER