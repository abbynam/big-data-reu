#ifndef DCA_KERNEL_RECONSTRUCTOR_H_
#define DCA_KERNEL_RECONSTRUCTOR_H_
//
// DCAKernelReconstructor.cpp
//
// Estimates the d


// standard libraies
#include <vector>
#include <string>
#include <valarray>
#include <memory>

// Root libraries
#include <TFile.h>
#include <TTree.h>
#include <TTimeStamp.h>
#include <TH3F.h>
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
// Class declarations
class DCAKernelReconstructor {

private:

  vector<ConicSection> conic_sections_; ///vector to store pointers to the parabolas and ellipses

  string parameters_file_path_;
  RunTimeParameters run_params_;
  shared_ptr<DCAKernelDensityEstimator> density_estimator_ptr_;
  
  string event_file_path_; /// path to file with gamma dections
  size_t max_number_cones_ ; /// Ignore cones exceeding this limit
  int density_points_per_cone_; /// number of random points in the density matrix per cone
  int test_points_per_cone_; /// number of random points used to findest densest per cone
  size_t number_cones_; /// smaller of max_number_cones_ and the actually number of cones in the event file


  shared_ptr<PhantomVolume> phantom_volume_ptr_;

  int bins_x_;
  int bins_y_;
  int bins_z_;
  float x_length_;
  float y_length_;
  float z_length_;


  //   I/O
  string output_root_file_path_;
  string output_folder_path_;
  TFile *output_root_file_;
  void open_output_file_( int numIterations);
  void save_results_( int test_points);
  void read_event_file_();
  int load_conic_sections_(TTree *tree);
  void load_events_();
  string setup_output_folder_();  



   //utilities
   char buffer_[501];
   TTimeStamp start_time_;
   inline double getRunningTime(){ return TTimeStamp().AsDouble() - start_time_.AsDouble();};
  
   //Gamma emission density related properties and methods
   void setup_density_matrix_();
   void populate_density_matrix_(const vector<shared_ptr<ConicSection> > &conicSections, DensityEstimator &density_estimator);
   void repopulate_density_matrix_(const vector<shared_ptr<ConicSection> > &conicSections,
							const DensityEstimator &old_de,
							DensityEstimator &new_de)  ;


public:
  DCAKernelReconstructor(const string &parameters_file_path);
  ~DCAKernelReconstructor();
  void StartExperiment();
};

}//end of namespace
#endif //DCA_KERNEL_RECONSTRUCTOR_H_
