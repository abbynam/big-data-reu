#ifndef STOCHASTIC_RECONSTRUCTOR_H_
#define STOCHASTIC_RECONSTRUCTOR_H_

// standard libraies
#include <vector>
#include <string>
#include <memory>

// PromptGamma libraries
#include "ComptonScatter.h"
#include "ConicSection.h"
#include "ReconstructionParabola.h"
#include "DensityEstimator.h"
#include "DensityMatrix.h"
#include "RunTimeParameters.h"
#include "ReconstructionEllipse.h"

//ROOT
#include "TTimeStamp.h"
#include "TFile.h"

using namespace std;
using namespace pg_tools;

namespace prompt_gamma_reconstruction{

/*! \brief Implements Stochastic origin ensembles algorithm for Compton Camera Reconstruction
 *  
 * Performs Compton Camera reconstruction by implementing the following steps.
 * (2) through (4) are managed by StartCalc():
 * <ol>
 * <li> load_events_() -- Reads the Compton Camera scatter events and converts
 *          them into Compton scattering cones </li>
 * <li> setup_density_matrix_() builds a DensityEstimator object</li>
 * <li> populate_density_matrix_() loads the Compton scatters into the density matrix</li>
 * <li> calculate_density_matrix_() performs the iterations of the SOE algorithm</li>
 * </ol>
 * @author Dennis Mackin
 */    
class StochasticReconstructor{

private:
    
    string parameters_file_path_;
    RunTimeParameters run_time_parameters_; //Parameters from config file
    
    //Data storage
    //vector<shared_ptr<ConicSection> > conic_sections_; ///vector to store pointers to the parabolas and ellipses
    vector<ConicSection> conic_sections_; ///vector to store pointers to the parabolas and ellipses
    vector<ComptonScatter> bad_compton_scatters_; ///vector to store the compton scatters which fail reconstruction
    vector <ConicSection *> _badEvents; ///Store pointer to events for which an initial spot in the fiducial region could not be found
    TTimeStamp *start_time_;
    
    ///Density estimator base class
    shared_ptr<DensityEstimator> density_estimator_ptr_;
    
    ///pointer to ROOT 3D histogram
    TH3F *density_hist_; ///pointer to ROOT 3D histogram

    //properties
    string event_file_path_; ///Pathh to Compton Camera output file
    int max_number_cones_; //limit number cones used in reconstruction

    int iterations_; //number of iterations for the steps where the density is smoothed out
//    int number_kernels_; //Number of events to use in density estimation
    int number_tries_for_random_; //Maximum attempts to find random on cone & in phantom
    int num_scatters_; //Total number of scatters in in the input file
    int num_cones_; //Number of gamma scatter cones used in reconstruction

    PhantomVolume phantom_volume_;  /// Volume of the reconstruct,

    //Number of bins for the phantom
    int x_bins_;
    int y_bins_;
    int z_bins_;
    
    // dimensions of the phantom
    double x_length_;
    double y_length_;
    double z_length_;

    // I/O files
    string output_root_file_path_;
    string output_folder_path_;
    TFile *output_root_file_;

    //utilities
    char str_buffer_[501];//used for sprint string building

    //Program step methods
    void load_events_();
    void setup_density_matrix_(float bins_scalar);
    //void populate_density_matrix_(const vector<shared_ptr<ConicSection> > &conicSections, DensityEstimator &density_estimator);
	void populate_density_matrix_(const vector<ConicSection> &conicSections, DensityEstimator &density_estimator);
    void  calculate_density(int num_randoms, int width_number, float temperature);

    //Output control methods
    string setup_output_folder_(const string &outputFolderPath);    
    void open_output_file_(const int numIterations);
    void save_density_info_(const int iterations, Bool_t saveEventData);
    
public:
    StochasticReconstructor(const string &parameters_file_path);
    ~StochasticReconstructor();
    void start_calc(); //Kicks off the density estimation
    
    /*! \brief Returns the time since object creation */
    inline double get_running_time(){ return TTimeStamp().AsDouble() - start_time_->AsDouble();};
};

}//end of namespace
#endif //STOCHASTIC_RECONSTRUCTOR_H_

