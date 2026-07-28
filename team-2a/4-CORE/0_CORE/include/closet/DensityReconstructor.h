#ifndef DENSITY_RECONSTRUCTOR_H_
#define DENSITY_RECONSTRUCTOR_H_

// standard libraies
#include <iostream>

#include <cstdio>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <string>
#include <valarray>
#include <memory>

//boost  libraries


// Root libraries
// #include <TGClient.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TTimeStamp.h>
#include <TH3.h>
#include <TImage.h>
#include <TH2D.h>


// private libraries
#include "ComptonScatter.h"
#include "ConicSection.h"
#include "ConicEnsemble.h"
#include "ReconstructionEllipse.h"
#include "ReconstructionParabola.h"
#include "StringUtils.h"
#include "DensityMatrix.h"
#include "utilities/FileUtils.h"
#include "RunTimeParameters.h"


using namespace std;
using namespace pg_tools;

namespace prompt_gamma_reconstruction{
// Class declarations
class DensityReconstructor {

private:

    //data structures
    shared_ptr<ConicEnsemble> ensemble_; ///ensemble is the collection of triple compton scatters from which the image is to be reconstructed
    double smoothing_parameter_;//kernel density smoothing parameter
    double width_; //kernel width used instead of calculated standard deviation

    //TOOLS
    shared_ptr<TTimeStamp> start_time_;

    //PARAMETERS retrieved at run time from parameters file
    string parameters_file_path_;
    RunTimeParameters run_time_parameters_;
    int sharpening_iterations_; //number of iterations for the steps where the density is smoothed out
    int smoothing_iterations_;  // number of tries to find the densest (aka sharpest) point
    int smooth_sharp_loops_; //number of times to repeat the smoothing then sharpening steps
    int number_kernels_; //Number of events to use in density estimation
    int num_cones_;
    string events_file_path_;


    //I/O
    string output_root_file_path_;
    string output_folder_path_;
    TFile *output_root_file_;

    //Utilities
    void open_output_file_(const string &name);
    void produce_output_(const string &output_name);
    string setup_output_folder_(const string &outputFolderPath);
    void produce_1D_plot_(TH1D &hist);
    void produce_2D_plot_(TH2D &hist);
    TCanvas *setup_canvas_();


public:
    DensityReconstructor(const string &parameters_file_path);
    ~DensityReconstructor();
    void StartExperiment();
    inline double getRunningTime() const{ return TTimeStamp().AsDouble() - start_time_->AsDouble();} ;
};

}//end of namespace
#endif //DENSITY_RECONSTRUCTOR_H_
