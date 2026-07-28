//
// Created by dsmackin on 12/30/15.
//

#ifndef TRUST11_SOEALGORITHM_H
#define TRUST11_SOEALGORITHM_H


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

    class SOEAlgorithm: public ImageAlgorithm{
    public:
        SOEAlgorithm(const vector<ConicSection> &conic_sections, shared_ptr<const PhantomVolume> phantom_volume): phantom_volume_(phantom_volume),
                                                                                                                  number_tries_for_random_(1000),
                                                                                                                  number_of_threads_(1),
                                                                                                                  inverse_square_parameter_(1.0)
        {
            setConicSections(conic_sections);
        };
        ~SOEAlgorithm() { };


        //Get 2D image data
        Image2D getImagePlane(size_t dimension, float depth) const;

        //Get 3D image data
        Image3D getImageVolume(size_t dimension) const;

        //Flexible function that can return the data as a string in any format
        string getDataAsString() const;
        string getSystemMatrixAsString(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax, size_t nz, float zmin, float zmax) const {

            //return system_matrix_ptr_->get3DDose(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
            return "getSystemMatrixAsString is not implemented in VectorAlgorithm.\n";
        };
        string getDataAsString(size_t nx, size_t ny, size_t nz) const;
        string getConicInformationAsString() const;

        //SETTERS
        void setNumberOfIterations(const size_t iterations){number_of_iterations_ = iterations;};

        void setConicSections(const vector<ConicSection> &conic_sections);

        void setDensityEstimator(shared_ptr<DensityEstimator> de){
            density_estimator_ptr_ = de;
            populate_density_matrix_(conic_sections_, *density_estimator_ptr_);
        }

        void setTuningParameters(size_t threads){
            number_of_threads_ = threads;
        }

        void setNumberOfTriesForRandom(size_t tries){ number_tries_for_random_ = tries; };
        void run();

    private:
        void populate_density_matrix_(const vector<ConicSection> &conicSections, DensityEstimator &density_estimator);
        void set_inverse_square_params_();

        //Image reconstruction steps
        void calculate_image_();
        size_t number_of_iterations_;

        //PROPERTIES
        vector<ConicSection> conic_sections_; ///vector to store pointers to the parabolas and ellipses
//        shared_ptr<DensityEstimator> density_estimator_ptr_;
        shared_ptr<const PhantomVolume> phantom_volume_;
        std::time_t start_time_;

        size_t number_tries_for_random_;
        size_t number_of_threads_;
        float inverse_square_parameter_;


    };
};

#endif //TRUST11_SOEALGORITHM_H
