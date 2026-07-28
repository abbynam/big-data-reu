#ifndef _OCTANE_EM_ALGORITHM
#define _OCTANE_EM_ALGORITHM

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
#include "DensityMatrix.h"
#include "ImageAlgorithm.h"
#include "OriginCone.h"

using namespace std;
namespace prompt_gamma_reconstruction{
    
    class OctaneEMAlgorithm: public ImageAlgorithm{
    public:
        OctaneEMAlgorithm(const vector<ConicSection> &conics, shared_ptr<const PhantomVolume> phantom_volume):
                phantom_volume_ptr_(phantom_volume), phantom_center_(0.0,0.0,0.0), SAD_(100.0), inverse_square_param_(1.0), number_threads_(1){

            setConicSections(conics);
            buildOriginConeArray(conic_sections_, origin_cones_);

        };
        ~OctaneEMAlgorithm() { };
        
        //Get 2D image data 
        Image2D getImagePlane(size_t dimension, float depth) const;
        
        //Get 3D image data
        Image3D getImageVolume(size_t dimension) const;

        string getConicInformationAsString() const;
        void setConicSections(const vector<ConicSection> &conic_sections);

        void buildOriginConeArray(const vector<ConicSection> &cs, vector<OriginCone> &oc);
        void setInterceptDCA(float d){intercept_dca_ = d; }
        void setPhantomCenter(float x, float y, float z){ phantom_center_ = PGVector3(x,y,z);};
        void setPhantomLength(float l){ phantom_length_ = l;};
        void setSAD(float sad){ SAD_ = 1.0/sad;};
        void setInverseSquareParameter(float rhs){ inverse_square_param_ = rhs;};
        void setNumberThreads(const size_t n){ number_threads_ = n;};
        void setNumberIterations(const size_t n){ number_iterations_ = n;};
        string get_event_record_(long event_num) const;

        void setSystemMatrixEstimator(const shared_ptr<DensityEstimator> density_ptr){
            system_matrix_ptr_ = density_ptr->clone();
            populate_system_matrix(origin_cones_);
        };

        string getSystemMatrixAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const {
            system_matrix_ptr_->print();
            return system_matrix_ptr_->get3DDose(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
        };

        void run();        
        
        
    private:
        void populate_density_matrix(const vector<OriginCone> &originCones);
        void populate_system_matrix(const vector<OriginCone> &originCones);
        std::vector<PGVector3>  get_octant_centers(const PGVector3 &center, const double length);
        vector<PGVector3> get_intercepts(const OriginCone &originCone, const PGVector3 &center, float length, float intercept_dca);

        //PROPERTIES
        vector<ConicSection> conic_sections_;
        vector<OriginCone> origin_cones_;

        shared_ptr<const PhantomVolume> phantom_volume_ptr_;
        float SAD_;
        shared_ptr<DensityEstimator> system_matrix_ptr_;
        float inverse_square_param_;
        size_t number_iterations_;
        size_t number_threads_;
        std::time_t start_time_;
        float intercept_dca_;
        float phantom_length_;
        PGVector3 phantom_center_;


    };
};
#endif // _OCTANE_EM_ALGORITHM
