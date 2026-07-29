#ifndef CONIC_ENSEMBLE_H_
#define CONIC_ENSEMBLE_H_
//
// ConicEnsemble.cc
// 


// standard libraies
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <string>

//ROOT libraries
#include "TTree.h"
#include "TTimeStamp.h"

// private libraries
#include "ComptonScatter.h"
#include "ConicSection.h"
#include "KernelConic.h"


using namespace std;


namespace prompt_gamma_reconstruction{
    
    // Class declarations
    class ConicEnsemble {

    private:

        //data structures
        vector<KernelConic> kernel_conic_sections_; ///vector to store pointers to the parabolas and ellipses

        //Density estimation parameters
        int number_cones_;
        double h_; //smoothing parameter
        double sigma_; //kernel width
        double coefficient_; //coefficient for each conic
        double exponential_denominator_; //store 1/(h*sigma) so we don't have to calculate for every cone.
        shared_ptr<const PhantomVolume> phantom_volume_;

        //methods
        double getDensityForPoint_(const double &distance) const;
        int loadConicSections_(TTree *tree);
        void readInConicSections_(const string &root_input_file_path);
        string setupOutputFolder_(const string &outputFolderPath);
        
        /**Calculate and store the Gaussian denominator to convert division 
         * into multiplication for performance improvement.
         */
        inline void setExponentialDenominator_(){exponential_denominator_ = 1.0/(h_*sigma_);};
        
    public:  
        double getDensity(const PGVector3 &point, const int number_conics) const;
        void setSmoothingParameterForEnsemble(int number_cones){setSmoothingParameterForEnsemble(0.0,number_cones);} ;
        void setSmoothingParameterForEnsemble(double h, int number_cones);
        
        inline void setCoefficient(){ coefficient_ = 1.0/(sigma_ *h_ *sqrt(2*M_PI)); setExponentialDenominator_(); };
        void setWidthForEnsemble(double width);
        size_t size() const{return kernel_conic_sections_.size();};

        ConicEnsemble(const string &root_input_file_path, double smoothing, double width);
        ~ConicEnsemble();
    };
}//end of namespace
#endif //CONIC_ENSEMBLE_H_
