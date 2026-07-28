#ifndef NEAREST_NEIGHBOR_RECONSTRUCTOR_H_
#define NEAREST_NEIGHBOR_RECONSTRUCTOR_H_
//
// NearestNeighborReconstructor.cpp
//
// Constructs 3D histogram of random points from scattering cones.
// Since billions of points are generated, the matrix should be able to
// handle submillimeter bin sizes.


// standard libraies
#include <iosfwd>
#include <cstdio>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <string>
#include <valarray>
#include <cassert>
#include <memory>


// Root libraries
// #include <TGClient.h>
#include <TCanvas.h>
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
#include "DCAKernelDensityEstimator.h"
#include "utilities/FileUtils.h"
#include "RunTimeParameters.h"


using namespace std;
using namespace pg_tools;

namespace prompt_gamma_reconstruction{

/*! \brief Experimental reconstruction algorithm for point sources.
 * 
 * Based on K nearest neighbors algorithms. Random points are selected 
 * in the phantom to be reconstructed. Then, for each the DCA to the 
 * random point is calculated. The cone values are averaged and then 
 * the process is repeated.
 * 
 * This method is promising for point sources and may have applications
 * for beams.
 * 
 * @author Dennis Mackin
 */ 
class NearestNeighborReconstructor {

private:

    //data structures
    vector<shared_ptr<ConicSection> > conic_sections_; ///vector to store pointers to the parabolas and ellipses
    vector<ComptonScatter> bad_compton_scatters_; ///vector to store the compton scatters which fail reconstruction
    //vector <ConicSection *> _badEvents; ///Store pointer to events fo which an initial spot in the fiducial region could not be found
    TTimeStamp start_time_;

    string parameters_file_path_;
    RunTimeParameters run_time_parameters_;
    RandomSingleton *ptr_random_;
    RandomSqrtSingleton *ptr_random_sqrt_;
    shared_ptr<DCAKernelDensityEstimator> density_estimator_ptr_;


    string event_file_path_; /// path to file with gamma dections
    size_t max_number_cones_ ; /// Ignore cones exceeding this limit
    size_t number_cones_;
    int density_points_per_cone_; /// number of random points in the density matrix per cone
    int test_points_per_cone_; /// number of random points used to findest densest per cone
    int move_type_; /// determines what criteria is used to choose between points on a cone during an iteration.
    int density_estimator_type_; /// Holds density estimator which can be a 3D histogram, an average shifted histogram,

    shared_ptr <PhantomVolume> phantom_volume_ptr_;

    Int_t bins_x_;
    Int_t bins_y_;
    Int_t bins_z_;
    double x_length_;
    double y_length_;
    double z_length_;


    //   I/O
    string output_root_file_path_;
    string output_folder_path_;
    TFile *output_root_file_;


    //utilities
    char buffer_[501];


    void  saveHistogram(int number_randoms);

    void open_output_file_( int numIterations);
    void save_density_info_( int test_points);
    void read_event_file_();
    int load_conic_sections_(TTree *tree);
    void load_events_();
    string setup_output_folder_(const string &outputFolderPath);
    void setup_density_matrix_();

    PGVector3 set_pca_values(vector<shared_ptr<ConicSection> > conic_sections, PGVector3 &point);
    PGVector3 get_average_point(vector<shared_ptr<ConicSection> > conic_sections);
    void populate_density_matrix_(const vector<shared_ptr<ConicSection> > &conicSections, DensityEstimator &density_estimator);
    void repopulate_density_matrix_(const vector<shared_ptr<ConicSection> > &conicSections,
                                                          const DensityEstimator &old_de,
                                                          DensityEstimator &new_de)  ;
    inline double getRunningTime(){ return TTimeStamp().AsDouble() - start_time_.AsDouble();};

public:
    NearestNeighborReconstructor(const string &parameters_file_path);

    ~NearestNeighborReconstructor();
    void StartExperiment();
};

}//end of namespace
#endif //NEAREST_NEIGHBOR_RECONSTRUCTOR_H_
