#ifndef _ABSTRACT_IMAGE_ALGORITHM_BUILDER
#define _ABSTRACT_IMAGE_ALGORITHM_BUILDER

//C++ standard includes
#include <vector>
#include <ctime>
#include <algorithm>
#include <memory>
#include <cstdlib>

// Custom includes
#include "RunTimeParameters.h"
#include "ImageAlgorithm.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{
        class AbstractImageAlgorithmBuilder{
    public:
        AbstractImageAlgorithmBuilder(){ };
        virtual ~AbstractImageAlgorithmBuilder(){};
        virtual shared_ptr<ImageAlgorithm> build(const vector<ConicSection> &conics, const RunTimeParameters &params) const = 0;        
    };
};
#endif // _ABSTRACT_IMAGE_ALGORITHM_BUILDER