#ifndef _KEM_ALGORITHM
#define _KEM_ALGORITHM

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
#include "OriginConesSoA.h"
#include "cuda_functions.h"

using namespace std;
namespace prompt_gamma_reconstruction{
    
    class KEMAlgorithm: public ImageAlgorithm{
    public:
        KEMAlgorithm(const vector<ConicSection> &conics): number_threads_(1), bandwidth_(1.0){

            setConicSections(conics);

            origin_cones_soa_ = getOriginConesSoA(conic_sections_);
        };
        ~KEMAlgorithm() { };
        
        //Get 2D image data 
        Image2D getImagePlane(size_t dimension, float depth) const;
        
        //Get 3D image data
        Image3D getImageVolume(size_t dimension) const;

        string getConicInformationAsString() const;
        void setConicSections(const vector<ConicSection> &conic_sections);

        float getDensity(const PGVector3 p, const vector<OriginCone> &cones);

        inline PGVector3 getBinCenter(size_t bin){
            if(bin >= x_bins_*y_bins_*z_bins_) cout << "bin: " << bin << endl;
            assert(bin < x_bins_*y_bins_*z_bins_);
            PGVector3 p(0.0, 0.0, 0.0);

            size_t zbin = bin/(x_bins_*y_bins_);
            size_t ybin = (bin - zbin*(x_bins_*y_bins_))/x_bins_;
            size_t xbin = bin - zbin*(x_bins_*y_bins_) - ybin*x_bins_;

            float x_step = (x_max_ - x_min_)/float(x_bins_);
            float y_step = (y_max_ - y_min_)/float(y_bins_);
            float z_step = (z_max_ - z_min_)/float(z_bins_);

            p.x = x_min_ + (xbin + 0.5)*x_step;
            p.y = y_min_ + (ybin + 0.5)*y_step;
            p.z = z_min_ + (zbin + 0.5)*z_step;

            return p;
        };

        vector<PGVector3> getBinCenters(){
            vector<PGVector3> centers(densities_.size());
            for(size_t i=0; i < centers.size(); ++i) centers[i] = this->getBinCenter(i);

            return centers;
        };

        void buildPhantomVolume(const pg_tools::RunTimeParameters &params){

            x_bins_ = params.get_int("X_BINS");
            x_max_ = params.get_double("X_MAX");
            x_min_ = params.get_double("X_MIN");

            y_bins_ = params.get_int("Y_BINS");
            y_max_ = params.get_double("Y_MAX");
            y_min_ = params.get_double("Y_MIN");

            z_bins_ = params.get_int("Z_BINS");
            z_max_ = params.get_double("Z_MAX");
            z_min_ = params.get_double("Z_MIN");

            for(float f: {x_min_, x_max_, (float)x_bins_, y_min_, y_max_, (float)y_bins_, z_min_, z_max_, (float)z_bins_})
                volume_grid_.push_back(f);

            densities_ = vector<float>(x_bins_*y_bins_*z_bins_, 1);
            system_matrix_ = vector<float>(densities_);
            cout<<"Sys mat now has "<<system_matrix_.size()<< " voxels . . ." << endl;
            populate_system_matrix(origin_cones_soa_);
        }

//        void buildOriginConeArray(const vector<ConicSection> &cs, vector<OriginCone> &oc);
        OriginConesSoA getOriginConesSoA(const vector<ConicSection> &cs);
        void setInterceptDCA(float d){intercept_dca_ = d; }
//        void setPhantomCenter(float x, float y, float z){ phantom_center_ = PGVector3(x,y,z);};
//        void setPhantomLength(float l){ phantom_length_ = l;};
//        void setSAD(float sad){ SAD_ = 1.0/sad;};
        void setBandwidth(float b){bandwidth_ = b;}
//        void setInverseSquareParameter(float rhs){ inverse_square_param_ = rhs;};
        void setNumberThreads(const size_t n){ number_threads_ = n;};
        void setNumberIterations(const size_t n){ number_iterations_ = n;};
        void setSystemMatrixScalar(const float s){ system_matrix_scalar_ = s;};
        string get_event_record_(long event_num) const;

        string getDataAsString() const;
        string getDataAsString(const vector<float> hist) const;
        string getDataAsString(size_t nx, size_t ny, size_t nz) const;
        string getDataAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const;

//        void setSystemMatrixEstimator(const shared_ptr<DensityEstimator> density_ptr){
//            system_matrix_ptr_ = density_ptr->clone();
//            populate_system_matrix(origin_cones_);
//        };

        string getSystemMatrixAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const {

            return getDataAsString(system_matrix_);
            return "getSystemMatrixAsString is not implemented in SBPAlgorithm.\n";
        };

        void run();        
        
        
    private:
        void populate_density_matrix(const OriginConesSoA &originCones);
        void populate_system_matrix(const OriginConesSoA &conesSoA);
//        std::vector<PGVector3>  get_octant_centers(const PGVector3 &center, const double length);
//        vector<PGVector3> get_intercepts(const OriginCone &originCone, const PGVector3 &center, float length, float intercept_dca);

        //PROPERTIES
        vector<ConicSection> conic_sections_;
        OriginConesSoA origin_cones_soa_;
        vector<float> densities_;
        vector<float> system_matrix_;
        float system_matrix_scalar_;

        shared_ptr<const PhantomVolume> phantom_volume_ptr_;
        float SAD_;
        shared_ptr<DensityEstimator> system_matrix_ptr_;
        float inverse_square_param_;
        size_t number_iterations_ = 1;
        size_t number_threads_;
        std::time_t start_time_;
        float intercept_dca_;
        float phantom_length_;
        float bandwidth_inv_;
        float bandwidth_;
        PGVector3 phantom_center_;

        float x_min_;
        float x_max_;
        size_t x_bins_;
        float y_min_;
        float y_max_;
        size_t y_bins_;
        float z_min_;
        float z_max_;
        size_t z_bins_;
        vector<float> volume_grid_;
    };
};
#endif // _KEM_ALGORITHM
