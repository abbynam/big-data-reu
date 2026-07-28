#ifndef HIST_BASED_RECONSTRUCTOR_H_
#define HIST_BASED_RECONSTRUCTOR_H_
//
// HistBasedReconstructor.cpp
//
// Constructs 3D histogram of random points from scattering cones.
// Since billions of points are generated, the matrix should be able to
// handle submillimeter bin sizes.


// standard libraies
#include <iosfwd>
#include <vector>
#include <string>
#include <valarray>
#include <memory>


// Root libraries
// #include <TGClient.h>
#include <TFile.h>
#include <TTree.h>
#include <TTimeStamp.h>
#include <TH3F.h>
#include <TImage.h>
#include <TH2F.h>


// PromptGamma libraries
#include "ComptonScatter.h"
#include "ConicSection.h"
#include "ReconstructionEllipse.h"
#include "ReconstructionParabola.h"
#include "StringUtils.h"
#include "DensityMatrix.h"

#include "utilities/FileUtils.h"
#include "RunTimeParameters.h"



using namespace std;
using namespace pg_tools;

namespace prompt_gamma_reconstruction{

    
/*! \brief Experimental reconstruction algorithm based on histograms
 * 
 * Rather than using a single representative point for each detected gamma
 * ray. A histogram is filled using a large number of weighted random points
 * from each cone.
 * 
 * When a large number of events have been detected this method has a lot
 * of smearing. It does not seem to be effective at precisely determining 
 * the proton beam range, and it has not been used beyond initial testing.
 * 
 * @author Dennis Mackin
 */   
class HistBasedReconstructor {

private:

    //data structures
    vector<shared_ptr<ConicSection> > conic_sections_; ///vector to store pointers to the parabolas and ellipses
    vector<ComptonScatter> bad_compton_scatters_; ///vector to store the compton scatters which fail reconstruction
    vector <ConicSection *> _badEvents; ///Store pointer to events fo which an initial spot in the fiducial region could not be found
    TTimeStamp start_time_;

    string parameters_file_path_;
    RunTimeParameters run_time_parameters_;
    RandomSingleton *ptr_random_;
    RandomSqrtSingleton *ptr_random_sqrt_;
    shared_ptr<DensityMatrix> density_matrix1_ptr_;
    shared_ptr<DensityMatrix> density_matrix2_ptr_;

    string event_file_path_; /// path to file with gamma dections
    int max_number_cones_ ; /// Ignore cones exceeding this limit
    int density_points_per_cone_; /// number of random points in the density matrix per cone
    int test_points_per_cone_; /// number of random points used to findest densest per cone
    int number_cones_; /// smaller of max_number_cones_ and the actually number of cones in the event file
    int move_type_; /// determines what criteria is used to choose between points on a cone during an iteration.
    int density_estimator_type_; /// Holds density estimator which can be a 3D histogram, an average shifted histogram,



    shared_ptr<PhantomVolume> phantom_volume_;

    Int_t bins_x_;
    Int_t bins_y_;
    Int_t bins_z_;
    double x_length_;
    double y_length_;
    double z_length_;
    double sum_of_weights;


    //   I/O
    string output_root_file_path_;
    string output_folder_path_;
    TFile *output_root_file_;


    //utilities
    char buffer_[501];

    // Recon Functions
    void choose_densest_points_();
    void choose_random_points_();

    //void  fillHistogram(int number_randoms);
    void  saveHistogram(int number_randoms);

    void open_output_file_( int numIterations);
    void save_density_info_( int test_points);
    void read_event_file_();
    int load_conic_sections_(TTree *tree);
    void load_events_();
    string setup_output_folder_(const string &outputFolderPath);
    void setup_density_matrix_();
    void populate_density_matrix_(const vector<shared_ptr<ConicSection> > &conicSections, DensityEstimator &density_estimator);

    void random_walk_algo(const vector<shared_ptr<ConicSection> > &conicSections,
                          size_t first_cone_index, size_t last_cone_index,
                          const DensityMatrix &old_de, DensityMatrix &new_de,
                          float step_size);

    void repopulate_density_matrix_parallel(
                const vector<shared_ptr<ConicSection> > &conicSections,
                const DensityMatrix &old_de,
                DensityMatrix &new_de,
                int number_threads);

    void repopulate_density_matrix_(const vector<shared_ptr<ConicSection> > &conicSections,
                const DensityMatrix &old_de,
                DensityMatrix &new_de);

    inline double getRunningTime(){ return TTimeStamp().AsDouble() - start_time_.AsDouble();};

public:
    HistBasedReconstructor(const string &parameters_file_path);

    ~HistBasedReconstructor();
    void StartExperiment();
};

}//end of namespace
#endif //HIST_BASED_RECONSTRUCTOR_H_
