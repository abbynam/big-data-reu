#ifndef _OCTANE_ALGORITHM
#define _OCTANE_ALGORITHM

//C++ standard includes
#include <vector>
#include <ctime>
#include <algorithm>
#include <memory>
#include <cstdlib>

// Custom includes
#include "ConicSection.h"
#include "RunTimeParameters.h"
#include "DensityEstimator.h"
#include "ImageAlgorithm.h"

using namespace std;
namespace prompt_gamma_reconstruction{
    
    class OctaneAlgorithm: public ImageAlgorithm{
    public:
        OctaneAlgorithm(const vector<ConicSection> &conics, shared_ptr<const PhantomVolume> phantom_volume):
                phantom_volume_ptr_(phantom_volume), phantom_center_(0.0,0.0,0.0), inverse_square_param_(1.0){
                setConicSections(conics);
        };
        ~OctaneAlgorithm() { };
        
        //Get 2D image data 
        Image2D getImagePlane(size_t dimension, float depth) const;
        
        //Get 3D image data
        Image3D getImageVolume(size_t dimension) const;

        string getConicInformationAsString() const;
        string getSystemMatrixAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const {

            //return system_matrix_ptr_->get3DDose(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
            return "getSystemMatrixAsString is not implemented in OctaneAlgorithm.\n";
        };

        void setConicSections(const vector<ConicSection> &conic_sections);
        void setInterceptDCA(float d){intercept_dca_ = d; }
        void setPhantomCenter(float x, float y, float z){ phantom_center_ = PGVector3(x,y,z);};
        void setPhantomLength(float l){ phantom_length_ = l;};
        void setInverseSquareParameter(float rhs){ inverse_square_param_ = rhs;};
        string get_event_record_(long event_num) const;
        
        void run();        
        
    private:
        void populate_density_matrix(const vector<ConicSection> &conics);
        std::vector<PGVector3>  get_octant_centers(const PGVector3 &center, const double length);
        vector<PGVector3> get_intercepts(const ConicSection &cs, const PGVector3 &center, float length, float intercept_dca);
        void set_intercepts(const ConicSection &cs, const PGVector3 &center, float length, float intercept_dca);

        size_t number_iterations_;

        //PROPERTIES
        vector<ConicSection> conic_sections_;

        shared_ptr<const PhantomVolume> phantom_volume_ptr_;
        std::time_t start_time_;
        float intercept_dca_;
        float phantom_length_;
        PGVector3 phantom_center_;
        float inverse_square_param_;
    };
};
#endif // _OCTANE_ALGORITHM
