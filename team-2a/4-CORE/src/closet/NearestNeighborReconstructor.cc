#include <sstream>
#include "omp.h"

//ROOT Includes
#include "TH1F.h"
#include "TH3F.h"

//PromptGamma includes
#include "NearestNeighborReconstructor.h"
#include "DCAKernelDensityEstimator.h"
#include "ConicEnsemble.h"
#include "RandomSingleton.h"
#include "RandomSqrtSingleton.h"
#include "EventsLoader.h"
#include "PGSimulatedEventsLoader.h"
#include "SimulatedCCEventsLoader.h"
#include "CSVEventsLoader.h"
#include "DetectorEffectsEventsLoader.h"
#include "utilities/utilities.h"
#include "EventsLoaderFactory.h"

using namespace prompt_gamma_reconstruction;
using namespace pg_tools;


NearestNeighborReconstructor::NearestNeighborReconstructor(const string &parameters_file_path):
    start_time_(TTimeStamp()), 
            parameters_file_path_(parameters_file_path),
            run_time_parameters_(parameters_file_path)
    {

    event_file_path_ = run_time_parameters_["EVENT_FILE_PATH"];
    max_number_cones_ = std::atoi(run_time_parameters_["MAX_NUM_CONES"].c_str());
    density_points_per_cone_ = std::atoi(run_time_parameters_["DENSITY_POINTS_PER_CONE"].c_str());


    setup_density_matrix_();

    output_folder_path_ =  setup_output_folder_(run_time_parameters_["OUTPUT_FOLDER_PATH"]);
    output_root_file_path_ = output_folder_path_ ;
    output_root_file_path_ +=  run_time_parameters_["FOLDER_DELIMITER"];
    output_root_file_path_ += run_time_parameters_["RUN_NAME"] + ".root"; //set in _openOutputFile

    load_events_();
    number_cones_ = this->conic_sections_.size();

    printf("processing %d cones . . .\n",  static_cast<int>(number_cones_));

    int totalVoxels = bins_x_ * bins_y_ * bins_z_;
    double eventsPerVox = (double)number_cones_ / (double)totalVoxels;
    printf("--- Number of Voxels: %d ---\n", totalVoxels);
    printf("--- Events per Voxel: %.2f ---\n", eventsPerVox);
    printf("--- Run Time: %f sec ---\n", getRunningTime());

    cout<<"completed construction of NearestNeighborReconstructor . . ."<<endl;
}


NearestNeighborReconstructor::~NearestNeighborReconstructor() {
    if(output_root_file_){
      output_root_file_->Write();
      output_root_file_->Close();
      delete output_root_file_;
    }

    cout<<"completed destruction of NearestNeighborReconstructor . . ."<<endl;
}


void NearestNeighborReconstructor::setup_density_matrix_(){

    x_length_ = strtod(run_time_parameters_["X_LENGTH"].c_str(), 0);
    bins_x_ = StringUtils::strtoi(run_time_parameters_["X_BINS"].c_str());

    phantom_volume_ptr_->x_max = x_length_ / 2.0;
    phantom_volume_ptr_->x_min = -1.0 * phantom_volume_ptr_->x_max;

    y_length_ = strtod(run_time_parameters_["Y_LENGTH"].c_str(), 0);
    bins_y_ = StringUtils::strtoi(run_time_parameters_["Y_BINS"].c_str());

    phantom_volume_ptr_->y_max = y_length_ / 2.0;
    phantom_volume_ptr_->y_min = -1.0 * phantom_volume_ptr_->y_max;

    z_length_ = strtod(run_time_parameters_["Z_LENGTH"].c_str(), 0);
    bins_z_ = StringUtils::strtoi(run_time_parameters_["Z_BINS"].c_str());

    phantom_volume_ptr_->z_max = z_length_ / 2.0;
    phantom_volume_ptr_->z_min = -1.0 * phantom_volume_ptr_->z_max;

}

void NearestNeighborReconstructor::load_events_(){

  int data_file_format = atoi(run_time_parameters_["DATA_FILE_FORMAT"].c_str());

  string gamma_tree_name = run_time_parameters_["GAMMA_TREE_NAME"];
  double percent_cone_phantom_overlap = strtod(run_time_parameters_["CONE_PHANTOM_OVERLAP_PERCENTAGE"].c_str(), 0);
  int number_tries_per_random_point = static_cast<int>(100 * (1.0/percent_cone_phantom_overlap));
  
  auto eventLoaderPtr = EventsLoaderFactory::create(data_file_format, event_file_path_, &run_time_parameters_,
                                                    phantom_volume_ptr_);
  eventLoaderPtr->LoadEvents(conic_sections_, number_tries_per_random_point);
  
  number_cones_ = static_cast<int>(conic_sections_.size());
  if(max_number_cones_ < number_cones_){
    conic_sections_.erase(conic_sections_.begin() + max_number_cones_, conic_sections_.end());
    cout<<"WARNING: Only using "<< max_number_cones_ <<" of " << conic_sections_.size() << " possible cones . . ." << endl;
    number_cones_ = max_number_cones_;
  }

}

PGVector3 NearestNeighborReconstructor::set_pca_values(vector<shared_ptr<ConicSection> > conic_sections, PGVector3 &point){
  size_t number_cones = conic_sections.size();
  PGVector3 pca;
  PGVector3 origin(0.0,0.0,0.0);
  PGVector3 average_point(0.0,0.0,0.0);  
  float dca;
  float cone_weight = 0.0;
  float total_weight = 0.0;  
  float epsilon =  1.0/static_cast<float>(number_cones);

  for(size_t iCone = 0; iCone < number_cones; ++iCone){
    dca = conic_sections[iCone]->getDistanceToPoint(point);
    pca = conic_sections[iCone]->getPointOfClosestApproach(point);
    conic_sections[iCone]->setLikelyOrigin(pca);
    
    cone_weight =  1.0/(dca + epsilon);
    average_point +=  PGVector3(conic_sections[iCone]->getLikelyOrigin()) * cone_weight;
    total_weight += cone_weight;
  }
  average_point *= 1.0/(total_weight);
  return average_point;  
}

PGVector3 NearestNeighborReconstructor::get_average_point(vector<shared_ptr<ConicSection> > conic_sections){
  size_t number_cones = conic_sections.size();
  PGVector3 average_point(0.0,0.0,0.0);

  for(size_t iCone = 0; iCone < number_cones; ++iCone){
    average_point += conic_sections[iCone]->getLikelyOrigin();
  }
  average_point *= 1.0/static_cast<float>(number_cones);
  return average_point;
}

// Draws function graphics in randomly choosen interval
void NearestNeighborReconstructor::StartExperiment() {

    auto convergence_tolerance = run_time_parameters_.get_double("CONVERGENCE_TOLERANCE");
  PGVector3 central_point(-90.0,-80.0,90.0);
  central_point = get_average_point(conic_sections_);
  int iterations = atoi(run_time_parameters_["ITERATIONS"].c_str());
  for(int i=0; i<iterations; ++i){
    central_point = set_pca_values(conic_sections_,central_point);
    central_point.y = 0.0;
    
    auto new_central_point = get_average_point(conic_sections_);
    if(central_point.getDistanceToPoint(new_central_point) < convergence_tolerance){
        break;
    }else{
        cout<<"Central point: "<<central_point.print()<< ", "<< new_central_point.print()<<", "<<central_point.getDistanceToPoint(new_central_point)<<endl;
    }
    central_point = new_central_point;
    cout<<"Central point: "<<central_point.print()<<endl;
  }
  save_density_info_(0);
}

void NearestNeighborReconstructor::open_output_file_( int number_iterations) {

  //////// Open file to store data
  vector<string> parts;
  vector<string> fileNameRootParts;
  string fileNameRoot;
  StringUtils::split( event_file_path_, "/", parts );
  if(parts.size() > 0){
    StringUtils::split( parts.back(), ".root", fileNameRootParts );
  }
  if(fileNameRootParts.size() > 0){
    fileNameRoot = fileNameRootParts.front();
  }
  sprintf(buffer_,"%s%s%s_%s_%d.root",
	  output_folder_path_.c_str(),
	  run_time_parameters_["FOLDER_DELIMITER"].c_str(),
	  run_time_parameters_["RUN_NAME"].c_str(),
	  fileNameRoot.c_str(),
	  number_iterations);

  output_root_file_ = new TFile(buffer_,"RECREATE");
  if(! output_root_file_->IsOpen()){
    string error_message("ERROR: failed to open ");
    error_message += buffer_;
    error_message += ".\nABORTING\n\n";
    throw runtime_error(error_message);
  };
}

void NearestNeighborReconstructor::save_density_info_( int numIterationsPerformed) {
  open_output_file_( numIterationsPerformed);
  assert(output_root_file_->IsOpen());
  output_root_file_->cd();

  /// Declare ntuples to store cone and final event data
  TTree output_tree("promptGammaOutput", "Gamma Reconstruction Information");

  float final_point[3], origin_true[3];
  float p1[3], p2[3], p3[3], p1_true[3], p2_true[3], p3_true[3];

  float incident_energy[3], incident_energy_true[3], initial_energy_true;
  float energy_deposit[3], energy_deposit_true[3];

  float  apex[3], axis[3];
  float angle[2], angle_true[2];
  float distance_closest_approach=0.0;
  PGVector3 pca;
  float point_closest_approach[3];


  output_tree.Branch("final_point", &final_point, "final_x/F:final_y/F:final_z/F");
  if( 0 == numIterationsPerformed){
    ///POSITION

    output_tree.Branch("origin", &origin_true, "origin_x/F:origin_y/F:origin_z/F");
    output_tree.Branch("position1", &p1, "p1_x/F:p1_y/F:p1_z/F");
    output_tree.Branch("position2", &p2, "p2_x/F:p2_y/F:p2_z/F");
    output_tree.Branch("position3", &p3, "p3_x/F:p3_y/F:p3_z/F");
    output_tree.Branch("position1true", &p1_true, "p1_x_true/F:p1_y_true/F:p1_z_true/F");
    output_tree.Branch("position2true", &p2_true, "p2_x_true/F:p2_y_true/F:p2_z_true/F");
    output_tree.Branch("position3true", &p3_true, "p3_x_true/F:p3_y_true/F:p3_z_true/F");

    ///ENERGY
    output_tree.Branch("initial_energy_true", &initial_energy_true, "initial_energy/F");
    output_tree.Branch("incident_energy", &incident_energy, "E0/F:E1/F:E2/F");
    output_tree.Branch("incident_energy_true", &incident_energy_true, "E0_true/F:E1_true/F:E2_true/F");

    output_tree.Branch("energy_deposit", &energy_deposit, "deposit_1/F:deposit_2/F");
    output_tree.Branch("energy_deposit_true", &energy_deposit_true, "deposit_1_true/F:deposit_2_true/F");

    ///SCATTERING ANGLES
    output_tree.Branch("angle", &angle, "angle_1/F:angle_2/F");
    output_tree.Branch("angle_true", &angle_true, "angle_1_true/F:angle_2_true/F");

    ///CALCULATED VALUES
    output_tree.Branch("distance_closest_approach", &distance_closest_approach, "distance_closest_approach/F");
    output_tree.Branch("point_closest_approach", &point_closest_approach, "pca_x/F:pca_y/F:pca_z/F");
    output_tree.Branch("apex", &apex, "apex_x/F:apex_y/F:apex_z/F");
    output_tree.Branch("axis", &axis, "axis_x/F:axis_y/F:axis_z/F");
  };


  // Create matrices to store final density data
  sprintf(buffer_,"density_iter_%05d", numIterationsPerformed);
  TH3F* density_hist_ = new TH3F(buffer_, "3D histogram of density matrix;z (mm);z (mm);y (mm)",
                                 bins_z_, phantom_volume_ptr_->z_min, phantom_volume_ptr_->z_max,
                                 bins_x_, phantom_volume_ptr_->z_min, phantom_volume_ptr_->z_max,
                                 bins_y_, phantom_volume_ptr_->y_min, phantom_volume_ptr_->y_max) ;
  sprintf(buffer_,"Event Origin Density (%d iterations);x (mm);y (mm);z (mm)", numIterationsPerformed);
  density_hist_->SetTitle(buffer_);


  PGVector3 origin(0.0,0.0,0.0);
  // Loop through cones to store data (apexPos, scatVect, scatAng, finalEventPos) into ntuples
  vector<shared_ptr<ConicSection> >::iterator iter = conic_sections_.begin();
  for( ; iter != conic_sections_.end(); ++iter){
    //populate variables mapped to tree branches
    ComptonScatter tmpScatter = (*iter)->getComptonScatter();
    MC_Truth mc_truth = (*iter)->getMCTruth();
    shared_ptr<Scatter> scatter_info = (*iter)->getScatterInfo(); 

    origin_true[0] = mc_truth.origin_true[0];
    origin_true[1] = mc_truth.origin_true[1];
    origin_true[2] = mc_truth.origin_true[2];

    final_point[0] = (float)(*iter)->getLikelyOrigin().x;
    final_point[1] = (float)(*iter)->getLikelyOrigin().y;
    final_point[2] = (float)(*iter)->getLikelyOrigin().z;

    ///Set the position information
    p1[0] = scatter_info->getScatterPositions()[0].x;
    p1[1] = scatter_info->getScatterPositions()[0].y;
    p1[2] = scatter_info->getScatterPositions()[0].z;

    p2[0] = scatter_info->getScatterPositions()[1].x;
    p2[1] = scatter_info->getScatterPositions()[1].y;
    p2[2] = scatter_info->getScatterPositions()[1].z;

    p3[0] = scatter_info->getScatterPositions()[2].x;
    p3[1] = scatter_info->getScatterPositions()[2].y;
    p3[2] = scatter_info->getScatterPositions()[2].z;

    p1_true[0] = scatter_info->getScatterPositionsTrue()[0].x;
    p1_true[1] = scatter_info->getScatterPositionsTrue()[0].y;
    p1_true[2] = scatter_info->getScatterPositionsTrue()[0].z;

    p2_true[0] = scatter_info->getScatterPositionsTrue()[1].x;
    p2_true[1] = scatter_info->getScatterPositionsTrue()[1].y;
    p2_true[2] = scatter_info->getScatterPositionsTrue()[1].z;

    p3_true[0] = scatter_info->getScatterPositionsTrue()[2].x;
    p3_true[1] = scatter_info->getScatterPositionsTrue()[2].y;
    p3_true[2] = scatter_info->getScatterPositionsTrue()[2].z;

    ///Set the energy information
    initial_energy_true = mc_truth.initial_energy;

    copy_vector_to_array(mc_truth.incident_energy, incident_energy_true);
    incident_energy[0] = scatter_info->getGammaEnergy();
    incident_energy[1] = scatter_info->getGammaEnergy() - scatter_info->getScatter1EnergyDeposit();
    incident_energy[2] = scatter_info->getGammaEnergy() - scatter_info->getScatter1EnergyDeposit() - scatter_info->getScatter2EnergyDeposit();
    copy_vector_to_array(mc_truth.incident_energy, incident_energy_true);

    energy_deposit[0] = scatter_info->getScatter1EnergyDeposit();
    energy_deposit[1] = scatter_info->getScatter2EnergyDeposit();
    energy_deposit[2] = 0.0;
    copy_vector_to_array(mc_truth.energy_deposition, energy_deposit_true);

    ///Set Scattering Information
    angle[0] = scatter_info->getTheta1Degrees();
    angle[1] = scatter_info->getTheta2Degrees();
    copy_vector_to_array(mc_truth.scattering_angle, angle_true);

    ///Set Calculated values
    axis[0] = (float)tmpScatter.getConeAxis().x;
    axis[1] = (float)tmpScatter.getConeAxis().y;
    axis[2] = (float)tmpScatter.getConeAxis().z;

    apex[0] = (float)tmpScatter.getConeApex().x;
    apex[1] = (float)tmpScatter.getConeApex().y;
    apex[2] = (float)tmpScatter.getConeApex().z;

    distance_closest_approach = (*iter)->getDistanceToPoint( PGVector3(origin_true[0], origin_true[1], origin_true[2]) );
    pca = (*iter)->getPointOfClosestApproach( PGVector3(origin_true[0], origin_true[1], origin_true[2]) );
    pca.set_xyz_array(point_closest_approach);

    output_tree.Fill();
  }

  printf("--- Storing data to root file: %s (%f)---\n", output_root_file_path_.c_str(), getRunningTime());
  output_tree.Write();
  // Write data to file


  printf("--- Storing data to root file: %s ---\n", output_root_file_path_.c_str());
  output_root_file_->Close();
  printf("-------- Saved State Histograms for %d iterations (%f)--------------\n", numIterationsPerformed, getRunningTime());
}
//////////////////////////////////////////////////////////////////////////////////////
string NearestNeighborReconstructor::setup_output_folder_(const string &outputFolderPath){
//////////////////////////////////////////////////////////////////////////////////////
///
/// Creates a new directory to store inputs and outputs for a looPDE run. This is useful
/// for keeping track of runs and their inputs and outputs.
/// - Creates the directory to store run information
/// - Copies parameters file to run folder
/// - Copies signal training data file to run folder
/// - Copies signal test data file to run folder
/// - Copies background training data file to run folder
/// - Copies background test data file to run folder
///
/// If the ouput folder already exists, then the exisiting files will be overwritten
/// except for the root output file which will be updated.
///
/// @param outpuFolderPath New run output folder will be created in this path
/// @returns the path to the new output file
//////////////////////////////////////////////////////////////////////////////////////

  stringstream ss;
  ss<<run_time_parameters_["OUTPUT_FOLDER_PATH"]<<"_DCA";
//   ss<<"_"<<test_points_per_cone_;
//   ss<<"_"<<bins_x_<<"x"<<bins_y_<<"x"<<bins_z_;
//   ss<<"_move" << move_type_;
//   ss<<"_densityEstimator" << density_estimator_type_;

  string output_folder(ss.str());
  char command[501];
  //check to see if folder exists
  sprintf(command, "%s %s",
	  run_time_parameters_["MAKE_DIR_COMMAND"].c_str(),
	  output_folder.c_str());
  cout<<"Creating directory "<<output_folder<<"."<<endl;
  printf("%s\n", command);

  int return_code = system(command);
  if(return_code != 0){
    printf("failed to make %s. \nDoes it exist? If not, there is a problem . . .\n",
      output_folder.c_str());
  }else{
    cout<<"Created directory "<< output_folder<<". . . \n";
  }

  //copy the parametersFile and the data files to the new folder
  FileUtils::fileCopy(parameters_file_path_, output_folder, run_time_parameters_["COPY_COMMAND"]);

  return output_folder;
}//end of SetUpOutputFolder
//////////////////////////////////////////////////////////////////////////////////////

