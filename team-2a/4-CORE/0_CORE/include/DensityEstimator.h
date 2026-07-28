#ifndef DENSITY_ESTIMATOR_H_
#define DENSITY_ESTIMATOR_H_
#define _USE_MATH_DEFINES

//standard C++ includes
#include <cmath>
#include <cstdio>
#include <valarray>
#include <vector>
#include <stdexcept>

#include <memory>


//PromptGamma includes
#include "ComptonScatter.h"
#include "PGVector3.h"

using namespace std;
namespace prompt_gamma_reconstruction{

/*! \brief Abstract base class for density estimators 
 *  which can used for image reconstruction.
 *  e.g. histograms, kernel methods, ASH, etc.
 * 
 * @author Dennis Mackin
 */    
class DensityEstimator{

  public:

    virtual ~DensityEstimator() {};

    virtual float getDensity(const PGVector3 &pos) const = 0;
    virtual vector<float> getDensities(const vector<PGVector3> &positions) const = 0;

    virtual void fill(const PGVector3 &pos, float weight) = 0;

    virtual void updateMatrix(const PGVector3 &oldPos, const PGVector3 &newPos) = 0;
    virtual void updateMatrix(const PGVector3 &oldPos, const PGVector3 &newPos, const float &weight) = 0;
    virtual std::shared_ptr<DensityEstimator> clone() const = 0;
    virtual DensityEstimator &operator=(float rhs)=0;

    virtual bool is_in_volume(const PGVector3 &point) const = 0;

    virtual void clear() = 0;
    virtual void print() const = 0;
    virtual void dump_bins() const {print();};
    virtual PGVector3 getBinCenter(size_t bin) const = 0;

    virtual string get3DDose() const { return string("get3DDose() has not been implemented for this density estimator.");};
    virtual string get3DDose(size_t nx, size_t ny, size_t nz) const { return string("get3DDose() has not been implemented for this density estimator.");};
    virtual string get3DDose(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax,  size_t nz, float zmin, float zmax) const
            { return string("get3DDose() has not been implemented for this density estimator.");}
};

}
#endif //DENSITY_ESTIMATOR_H_
