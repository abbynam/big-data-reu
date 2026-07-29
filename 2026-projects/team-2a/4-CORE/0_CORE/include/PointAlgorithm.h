//
// Created by dsmackin on 12/30/15.
//

#ifndef TRUST11_POINTALGORITHM_H
#define TRUST11_POINTALGORITHM_H


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

    class PointAlgorithm: public ImageAlgorithm{
    public:
        PointAlgorithm(const vector<ConicSection> &conic_sections, shared_ptr<const PhantomVolume> phantom_volume, float inverse_square_parameter): phantom_volume_(phantom_volume),
                                                                                                                  gaussian_width_(10000.0),
                                                                                                                    offset_x_(0.0), offset_y_(0.0),
                                                                                                                    number_tries_for_random_(1000),
                                                                                                                    number_of_threads_(1),
                                                                                                                     temperature_(1.0),
                                                                                                                    inverse_square_parameter_(inverse_square_parameter)
        {
            setConicSections(conic_sections);
            set_inverse_square_params_();
        };
        ~PointAlgorithm() { };


        //Get 2D image data
        Image2D getImagePlane(size_t dimension, float depth) const;

        //Get 3D image data
        Image3D getImageVolume(size_t dimension) const;

        //Flexible function that can return the data as a string in any format
//        string getDataAsString() const;
//        string getDataAsString(size_t nx, size_t ny, size_t nz) const;
//        string getDataAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const;
        string getConicInformationAsString() const;
        string getSystemMatrixAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const {

            //return system_matrix_ptr_->get3DDose(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
            return "getSystemMatrixAsString is not implemented in PointAlgorithm.\n";
        };
        //SETTERS
        void setNumberOfIterations(const size_t iterations){number_of_iterations_ = iterations;};

        void setConicSections(const vector<ConicSection> &conic_sections);

        void setDensityEstimator(shared_ptr<DensityEstimator> de){
            density_estimator_ptr_ = de;
            populate_density_matrix_(conic_sections_, *density_estimator_ptr_);
        }

        void setEventMultiplier(size_t multiplier){
            cout<<"Creating "<<(multiplier - 1)<<" duplicates of the events."<<endl;
            event_multiplier_ = multiplier;
            auto original_size = conic_sections_.size();
            cout<<"Old size: "<< conic_sections_.size()<<endl;

            conic_sections_.reserve(event_multiplier_ * original_size);

            cout<<"New size (after reserve): "<< conic_sections_.size()<<endl;
            for(size_t i=1; i<event_multiplier_; ++i){
                std::copy_n(conic_sections_.begin(), original_size, std::back_inserter(conic_sections_));
            }

            cout<<"New size (after fill): "<< conic_sections_.size()<<endl;
        }


        void setTuningParameters(size_t threads, float temperature){
            number_of_threads_ = threads;
            temperature_ = temperature;
        }

        void set_inverse_square_params_() {
            for(auto iter = conic_sections_.begin(); iter != conic_sections_.end(); ++iter){
                iter->setInverseSquareParam(this->inverse_square_parameter_);
            }
        };

        void setNumberOfTriesForRandom(size_t tries){ number_tries_for_random_ = tries; };
        void run();

    private:
        void populate_density_matrix_(const vector<ConicSection> &conicSections, DensityEstimator &density_estimator);

        //Image reconstruction steps
        void calculate_image_();
        size_t number_of_iterations_;

        //PROPERTIES
        vector<ConicSection> conic_sections_; ///vector to store pointers to the parabolas and ellipses
//        shared_ptr<DensityEstimator> density_estimator_ptr_;
        shared_ptr<const PhantomVolume> phantom_volume_;
        std::time_t start_time_;
        float gaussian_width_;
        float offset_x_;
        float offset_y_;
        size_t number_tries_for_random_;
        size_t number_of_threads_;
        float temperature_;
        size_t event_multiplier_;
        float inverse_square_parameter_;
    };
};

#endif //TRUST11_POINTALGORITHM_H
