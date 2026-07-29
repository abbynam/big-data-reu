#ifndef _IMAGE_ALGORITHM
#define _IMAGE_ALGORITHM

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

using namespace std;
namespace prompt_gamma_reconstruction{
    
    template<typename T>
    struct Image2D_T{
        Image2D_T(const vector<T> &start, const vector<T> &step, const vector<size_t> &n) :
            start_pos(start), step_size(start), number_pixels(n){ };
        
        const size_t number_of_dimensions = 2;
        vector<T> start_pos; //center of first pixel
        vector<T> step_size; //distance between pixel centers
        vector<size_t> number_pixels{number_of_dimensions, 0u};    
        
        vector<std::tuple<size_t, size_t, T>> image_data;  
    }; 
    typedef Image2D_T<float> Image2D;
    
    template<typename T>
    struct Image3D_T{
        Image3D_T(const vector<T> &start, const vector<T> &step, const vector<size_t> &n) :
            start_pos(start), step_size(start), number_bins(n) {};
        
        const size_t number_of_dimensions = 3;
        vector<T> start_pos{number_of_dimensions, 0}; //center of first pixel
        vector<T> step_size{number_of_dimensions, 0};
        vector<size_t> number_bins{number_of_dimensions, 0u};    
        
        vector<std::tuple<size_t, size_t, size_t, T>> image_volume;  
    };        
    typedef Image3D_T<float> Image3D;
    
    class ImageAlgorithm{
    public:
        ImageAlgorithm(): p_dca_(0.0,0.0,0.0), start_time_(std::time(nullptr)) { };
        virtual ~ImageAlgorithm() { };
        
        virtual void run() = 0;
        
        //Get 2D image data 
        virtual Image2D getImagePlane(size_t dimension, float depth) const = 0;
        
        //Get 3D image data
        virtual Image3D getImageVolume(size_t dimension) const = 0;
        
        //Flexible function that can return the data as a string in any format
        virtual string getDataAsString() const {return density_estimator_ptr_->get3DDose(); };
        virtual string getDataAsString(size_t nx, size_t ny, size_t nz) const { return density_estimator_ptr_->get3DDose(nx, ny, nz);};
        virtual string getDataAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const {
            cout<<"ImageAlgorithm::getDataAsString: getting dose from density estimator . . ." << endl;
            density_estimator_ptr_->print();
            return density_estimator_ptr_->get3DDose(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
        };
//        virtual string getSystemMatrixAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const {
//            cout<<"ImageAlgorithm::getDataAsString: getting dose from denstiy estimator . . ." << endl;
//            return "Not yet implemented.";
//        };

        virtual string getSystemMatrixAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const = 0;
        //virtual string getConicInformationAsString() const { return "Not yet implemented."; };
        virtual string getConicInformationAsString() const = 0;
        virtual string getConicInformationAsString(size_t nx, size_t ny, size_t nz) const { return getConicInformationAsString(); };
        
        
        //SETTERS
        virtual void setConicSections(const vector<ConicSection> &conic_sections) = 0;   
        
        void setRegionOfReconstruction(const PhantomVolume &phantom, vector<unsigned int> bins){
            phantom_volume_ = phantom;
            bins_ = bins;
        }
        
        void setDensityEstimator(const shared_ptr<DensityEstimator> density_ptr){
            density_estimator_ptr_ = density_ptr;
        }

        void setDCACenter(PGVector3 p){ p_dca_ = p;};
        
    protected:
        std::time_t get_running_time(){ return std::time(nullptr) - start_time_;};
        shared_ptr<DensityEstimator> density_estimator_ptr_;
        PGVector3 p_dca_;

    private:

        void calculate_image_();

        //PROPERTIES
        PhantomVolume phantom_volume_;
        vector<unsigned int> bins_;

        std::time_t start_time_;
    };
};
#endif // _IMAGE_ALGORITHM