#ifndef _CUDA_FUNCTIONS
#define _CUDA_FUNCTIONS

#include<vector>
#include "PGVector3.h"
#include "DensityMatrix.h"
#include "OriginCone.h"
#include "OriginConesSoA.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

void populate_density_matrix_cuda(vector<float> &density_matrix, const float [], const OriginConesSoA &cone_soa, const float bandwidth);

#endif // _CUDA_FUNCTIONS