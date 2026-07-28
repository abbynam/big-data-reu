#ifndef _DCA_LINE_ALGORITHM
#define _DCA_LINE_ALGORITHM

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
    
    class DCALineAlgorithm: public ImageAlgorithm{
    public:
        DCALineAlgorithm(const vector<ConicSection> &conics, shared_ptr<const PhantomVolume> phantom_volume):
                phantom_volume_ptr_(phantom_volume){
                setConicSections(conics);
        };
        ~DCALineAlgorithm() {cout<<"Destroying DCALineAlgorithm . . ."<<endl;};
        
        //Get 2D image data 
        Image2D getImagePlane(size_t dimension, float depth) const;
        
        //Get 3D image data
        Image3D getImageVolume(size_t dimension) const;
        
        
        //Flexible function that can return the data as a string in any format
        string getDataAsString() const;
        string getSystemMatrixAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const {

            //return system_matrix_ptr_->get3DDose(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
            return "getSystemMatrixAsString is not implemented in DCALineAlgorithm.\n";
        };
        string getConicInformationAsString() const;
        void setConicSections(const vector<ConicSection> &conic_sections);
        void setPoints(const PGVector3 &p1, const PGVector3 &p2){ p1_ = p1; p2_ = p2;};
        void setPoint1(const PGVector3 &p){ p1_ = p;};
        void setPoint2(const PGVector3 &p){ p2_ = p;};
        void setNumberOfThreads(const size_t &n){number_of_threads_ = n;};
        string get_event_record_(long event_num) const;
        
        void run();        
        
        
    private:
        void populate_density_matrix(const vector<ConicSection> &conics);

        //PROPERTIES
        vector<ConicSection> conic_sections_;
        //shared_ptr<DensityEstimator> density_estimator_ptr_;
        shared_ptr<const PhantomVolume> phantom_volume_ptr_;
        std::time_t start_time_;
        PGVector3 p1_; //p1 and p2 determine the line used for the recontstruction
        PGVector3 p2_;
        size_t number_of_threads_;
    };
};
#endif // _DCA_LINE_ALGORITHM
