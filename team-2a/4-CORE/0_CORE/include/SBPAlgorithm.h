#ifndef _SBP_ALGORITHM
#define _SBP_ALGORITHM

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
#include "OriginCone.h"



using namespace std;
namespace prompt_gamma_reconstruction{

    class SBPAlgorithm: public ImageAlgorithm{
    public:
        SBPAlgorithm(const vector<ConicSection> &conics):bandwidth_(1.0){
                setConicSections(conics);
        };
        ~SBPAlgorithm() { };

        string getDataAsString() const;
        virtual string getDataAsString(size_t nx, size_t ny, size_t nz) const;
        virtual string getDataAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const;

        //Get 2D image data
        Image2D getImagePlane(size_t dimension, float depth) const;

        //Get 3D image data
        Image3D getImageVolume(size_t dimension) const;

        string getConicInformationAsString() const;
        string getSystemMatrixAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const {

            //return system_matrix_ptr_->get3DDose(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
            return "getSystemMatrixAsString is not implemented in SBPAlgorithm.\n";
        };
        void setConicSections(const vector<ConicSection> &conic_sections);
        vector<PGVector3> getBinCenters();
        void setBandwidth(float b){bandwidth_ = b; bandwidth_inv_ = 1.0/bandwidth_;}
//        float getDensity(const PGVector3 p, const vector<OriginCone> &cones);

        string get_event_record_(long event_num) const;

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

            densities_ = vector<float>(x_bins_*y_bins_*z_bins_, 12.0f);

            for(float f: {x_min_, x_max_, (float)x_bins_, y_min_, y_max_, (float)y_bins_, z_min_, z_max_, (float)z_bins_})
                volume_grid_.push_back(f);
        }

        void run();


    private:
        void populate_density_matrix(const vector<OriginCone> &cones);
//        std::vector<PGVector3>  get_octant_centers(const PGVector3 &center, const double length);
//        vector<PGVector3> get_intercepts(const ConicSection &cs, const PGVector3 &center, float length, float intercept_dca);
//        void set_intercepts(const ConicSection &cs, const PGVector3 &center, float length, float intercept_dca);

        size_t number_iterations_;

        //PROPERTIES
        vector<ConicSection> conic_sections_;
        vector<OriginCone> cones_;
        vector<float> densities_;

        float bandwidth_;
        float bandwidth_inv_;
        std::time_t start_time_;

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

//    __global__ void helloFromGPU(void){
//
//        for(size_t i=0; i<100; ++i){
//            printf("Hello World %d from GPU (%d, %d, %f)!\n", i, blockIdx.x, threadIdx.x, cos(M_PI * (float)threadIdx.x/(1.0f+(float)blockIdx.x)));
//        }
//    }
};
#endif // _SBP_ALGORITHM
