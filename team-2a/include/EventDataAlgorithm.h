#ifndef _EVENT_DATA
#define _EVENT_DATA

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
    
    class EventDataAlgorithm: public ImageAlgorithm{
    public:
        EventDataAlgorithm(const vector<ConicSection> &conics, shared_ptr<const PhantomVolume> phantom_volume):
                phantom_volume_ptr_(phantom_volume){
                setConicSections(conics);
        };
        ~EventDataAlgorithm() { };
        
        //Get 2D image data 
        Image2D getImagePlane(size_t dimension, float depth) const;
        
        //Get 3D image data
        Image3D getImageVolume(size_t dimension) const;
        
        //Flexible function that can return the data as a string in any format
        string getDataAsString() const;
        string getSystemMatrixAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const {

            //return system_matrix_ptr_->get3DDose(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
            return "getSystemMatrixAsString is not implemented in EventDataAlgorithm.\n";
        };
        string getConicInformationAsString() const;
        void setConicSections(const vector<ConicSection> &conic_sections) {conic_sections_ = conic_sections;};

//        void setDCACenter(PGVector3 p){ p_dca_ = p;}; ///Posize_t used for DCA calculations
        
        void run();

    private:
        //PROPERTIES
        vector<ConicSection> conic_sections_;
        shared_ptr<const PhantomVolume> phantom_volume_ptr_;
//        PGVector3 p_dca_;
        std::time_t start_time_;
    };
};
#endif // _EVENT_DATA
