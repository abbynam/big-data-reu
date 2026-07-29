#ifndef _CO60_ALGORITHM
#define _CO60_ALGORITHM

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
    
    class Co60Algorithm: public ImageAlgorithm{
    public:
        Co60Algorithm(const vector<ConicSection> &conics, shared_ptr<const PhantomVolume> phantom_volume):
                phantom_volume_ptr_(phantom_volume){

                cout<<"Creating Co60Algorithm . . ."<<endl;
                setConicSections(conics);
        };
        ~Co60Algorithm() {cout<<"Destroying Co60Algorithm . . ."<<endl;};
        
        //Get 2D image data 
        Image2D getImagePlane(size_t dimension, float depth) const;
        
        //Get 3D image data
        Image3D getImageVolume(size_t dimension) const;
        
        
        //Flexible function that can return the data as a string in any format
        string getDataAsString() const;
        string getConicInformationAsString() const;
        string getSystemMatrixAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const {

            //return system_matrix_ptr_->get3DDose(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
            return "getSystemMatrixAsString is not implemented in Co60Algorithm.\n";
        };

        void setConicSections(const vector<ConicSection> &conic_sections){
            conic_sections_ = conic_sections;
            conic_sections_corrected_ = correct_conic_sections_(conic_sections_);
        };
        void setOrigin(const PGVector3 &p){ origin_ = p;};

        void setNumberOfThreads(const size_t &n){number_of_threads_ = n;};
        void setMinScatterEnergy(const float e){min_scatter_energy_ = e;};
        void setMinEventEnergy(const float e){min_event_energy_ = e;};

        void setImageAlgorithm(shared_ptr<ImageAlgorithm> algo){ algo_ = algo;}
        void run();

    private:
        vector<ConicSection> correct_conic_sections_(vector<ConicSection> &conic_sections);
        shared_ptr<ConicSection> get_conic_section_(float E1, float E2, PGVector3 &p1, PGVector3 &p2, size_t event_number);

        string get_event_record_(const ConicSection &cs, const ConicSection &cs_corrected) const;

        //PROPERTIES
        vector<ConicSection> conic_sections_;
        vector<ConicSection> conic_sections_corrected_;
        shared_ptr<const PhantomVolume> phantom_volume_ptr_;
        shared_ptr<ImageAlgorithm> algo_;
        std::time_t start_time_;
        PGVector3 origin_;
        float min_scatter_energy_;
        float min_event_energy_;
        size_t number_of_threads_;
    };
};
#endif // _CO60_ALGORITHM
