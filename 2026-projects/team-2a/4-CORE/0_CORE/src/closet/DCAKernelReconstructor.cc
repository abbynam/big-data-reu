// //
// DCAKernelReconstructor
//
// v1: Created based on StochasticReconstructor
//  ---- added by Dennis Mackin 2010-09-20
// v4: designed to compile to executable
// v5: replaced GUI interface with a config file

// Automatically compiles code / correct format: root reconGUIv3.C+
//  - must be in same directory (would like to change this)
#include <sstream>
//#include "omp.h"

//ROOT Includes
#include "TH1F.h"
#include "TH3F.h"
#include "TLine.h"

//PromptGamma includes
#include "DCAKernelReconstructor.h"
#include "DCAKernelDensityEstimator.h"
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


DCAKernelReconstructor::DCAKernelReconstructor(const string &parameters_file_path):
        start_time_(TTimeStamp()),
        parameters_file_path_(parameters_file_path),
        run_params_(parameters_file_path){

    max_number_cones_ = static_cast<size_t>(run_params_.get_int("MAX_NUM_CONES"));

    output_folder_path_ =  setup_output_folder_();
    output_root_file_path_ = output_folder_path_ + run_params_["FOLDER_DELIMITER"] + run_params_["RUN_NAME"] + ".root"; 

    setup_density_matrix_();
    load_events_();

    printf("processing %zu cones . . .\n",  this->conic_sections_.size());

    int totalVoxels = bins_x_ * bins_y_ * bins_z_;
    double eventsPerVox = (double)this->conic_sections_.size() / (double)totalVoxels;

    printf("--- Number of Voxels: %d ---\n", totalVoxels);
    printf("--- Events per Voxel: %.2f ---\n", eventsPerVox);
    printf("--- Run Time: %f sec ---\n", getRunningTime());

    cout<<"completed construction of DCAKernelReconstructor . . ."<<endl;
}


DCAKernelReconstructor::~DCAKernelReconstructor() {

  if(output_root_file_){
    output_root_file_->Write();
    output_root_file_->Close();
    delete output_root_file_;
  }

  cout<<"completed destruction of DCAKernelReconstructor . . ."<<endl;
}


void DCAKernelReconstructor::setup_density_matrix_(){

  x_length_ = run_params_.get_float("X_LENGTH");
  bins_x_ = run_params_.get_int("X_BINS");
  phantom_volume_ptr_->x_max = x_length_ / 2.0;
  phantom_volume_ptr_->x_min = -1.0 * phantom_volume_ptr_->x_max;

  y_length_ = run_params_.get_float("Y_LENGTH");
  bins_y_ = run_params_.get_int("Y_BINS");
  phantom_volume_ptr_->y_max = y_length_ / 2.0;
  phantom_volume_ptr_->y_min = -1.0 * phantom_volume_ptr_->y_max;

  z_length_ = run_params_.get_float("Z_LENGTH");
  bins_z_ = run_params_.get_int("Z_BINS");
  phantom_volume_ptr_->z_max = z_length_ / 2.0;
  phantom_volume_ptr_->z_min = -1.0 * phantom_volume_ptr_->z_max;


  float bandwidth =  run_params_.get_float("KERNEL_BANDWIDTH");
  int num_threads = run_params_.get_int("NUMBER_OF_THREADS");

  density_estimator_ptr_ = make_shared<DCAKernelDensityEstimator>(
          &conic_sections_, bandwidth, num_threads,
          phantom_volume_ptr_->x_min, phantom_volume_ptr_->x_max, bins_x_,
          phantom_volume_ptr_->y_min, phantom_volume_ptr_->y_max, bins_y_,
          phantom_volume_ptr_->z_min, phantom_volume_ptr_->z_max, bins_z_);

  density_estimator_ptr_->print();
}

void DCAKernelReconstructor::load_events_(){
  
  int data_file_format = run_params_.get_int("DATA_FILE_FORMAT");

  double percent_cone_phantom_overlap = run_params_.get_double("CONE_PHANTOM_OVERLAP_PERCENTAGE");
  int number_tries_per_random_point = static_cast<int>(100 * (1.0/percent_cone_phantom_overlap));

  auto eventsLoaderPtr = EventsLoaderFactory::create(data_file_format, run_params_["EVENT_FILE_PATH"], &run_params_,
                                                     phantom_volume_ptr_);
  vector< shared_ptr<ConicSection>> cs;
  eventsLoaderPtr->LoadEvents(cs, number_tries_per_random_point);
  
  //copy the pointers to ConicSections into an array of conicSections
  // to reduce the number of cache misses.
  conic_sections_.reserve(cs.size());
  for(size_t i = 0; i < cs.size(); ++i){
      conic_sections_.push_back(*cs[i]);
  }

  if(max_number_cones_ < conic_sections_.size()){
    cout<<"WARNING: Only using "<< max_number_cones_ <<" of " << conic_sections_.size() << " possible cones . . ." << endl;
  }
  
}


// Draws function graphics in randomly choosen interval
void DCAKernelReconstructor::StartExperiment() {
  save_results_(0);
}

void DCAKernelReconstructor::open_output_file_( int number_iterations) { 

  //////// Open file to store data
  vector<string> fileNameRootParts;
  string fileNameRoot;
  
  vector<string> parts;
  StringUtils::split( event_file_path_, "/", parts );
  if(parts.size() > 0){
    StringUtils::split( parts.back(), ".root", fileNameRootParts );
  }
  if(fileNameRootParts.size() > 0){
    fileNameRoot = fileNameRootParts.front();
  }
  sprintf(buffer_,"%s%s%s_%s_%d.root",
	  output_folder_path_.c_str(),
	  run_params_["FOLDER_DELIMITER"].c_str(),
	  run_params_["RUN_NAME"].c_str(),
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

void DCAKernelReconstructor::save_results_( int numIterationsPerformed) {
    open_output_file_( numIterationsPerformed);
    assert(output_root_file_->IsOpen());
    output_root_file_->cd();

    TTree output_tree("promptGammaOutput", "Gamma Reconstruction Information");

    float final_point[3], origin_true[3];
    float p1[3], p2[3], p3[3], p1_true[3], p2_true[3], p3_true[3];

    //when reconstructing Monte Carlo data, we know the true incident energy and we store
    // in incident_energy_true.
    float incident_energy[3], incident_energy_true[3], initial_energy_true, energy_lost;
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

        output_tree.Branch("energy_deposit", &energy_deposit, "deposit_1/F:deposit_2/F:deposit_3/F");
        output_tree.Branch("energy_deposit_true", &energy_deposit_true, "deposit_1_true/F:deposit_2_true/F");
        output_tree.Branch("energy_lost", &energy_lost, "energy_lost/F");

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
    for(auto iter = conic_sections_.begin(); iter != conic_sections_.end(); ++iter){
        //populate variables mapped to tree branches
        ComptonScatter tmpScatter = iter->getComptonScatter();
        MC_Truth mc_truth = iter->getMCTruth();
        shared_ptr<Scatter> scatter_info = iter->getScatterInfo();

        origin_true[0] = mc_truth.origin_true[0];
        origin_true[1] = mc_truth.origin_true[1];
        origin_true[2] = mc_truth.origin_true[2];

        final_point[0] = (float)iter->getLikelyOrigin().x;
        final_point[1] = (float)iter->getLikelyOrigin().y;
        final_point[2] = (float)iter->getLikelyOrigin().z;

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
        energy_deposit[2] = scatter_info->getScatter3EnergyDeposit();
        copy_vector_to_array(mc_truth.energy_deposition, energy_deposit_true);

        energy_lost = scatter_info->getEnergyLost();

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

        distance_closest_approach = iter->getDistanceToPoint( PGVector3(origin_true[0], origin_true[1], origin_true[2]) );
        pca = iter->getPointOfClosestApproach( PGVector3(origin_true[0], origin_true[1], origin_true[2]) );
        pca.set_xyz_array(point_closest_approach);

        output_tree.Fill();
    }

    printf("--- Storing data to root file: %s (%f)---\n", output_root_file_path_.c_str(), getRunningTime());
    output_tree.Write();
    // Write data to file 
  
   char buffer[1000];
   string web_folder = run_params_["WEB_FOLDER"].c_str();
   double z_min = strtod(run_params_["PLOT_Z_MIN"].c_str(),0);
   double z_max = strtod(run_params_["PLOT_Z_MAX"].c_str(),0);
   double tank_min = strtod(run_params_["TANK_MIN"].c_str(),0);
   double tank_max = strtod(run_params_["TANK_MAX"].c_str(),0);
   double tank_mid = strtod(run_params_["TANK_MID"].c_str(),0);
   double y_plane_min = strtod(run_params_["Y_PLANE_MIN"].c_str(),0);
   double y_plane_max = strtod(run_params_["Y_PLANE_MAX"].c_str(),0);
   double y_plane_step = strtod(run_params_["Y_PLANE_STEP"].c_str(),0);
   int z_bins = atoi(run_params_["PLOT_Z_BINS"].c_str());
   
  //save density matrix
   TCanvas *canvas = new TCanvas("myCan", "tuCan",600,600);

    for(int y = y_plane_min; y <= y_plane_max; y+=y_plane_step){
         auto h2D = density_estimator_ptr_->getRootHist_xz((float)y);
         h2D.GetMaximum();
         h2D.SetStats(0);
         if(h2D.GetMaximum() > 0.0) {
             h2D.Scale(1.0/h2D.GetMaximum());
         }
         h2D.Draw("COLZ");
         h2D.Write();
         sprintf(buffer,"%s/%s_%d.png", web_folder.c_str(), h2D.GetName(),y);
         canvas->SaveAs(buffer);   
    }   
   
   TH2* hist2D = density_estimator_ptr_->getRootHist_xy();
   hist2D->Draw("CONTZ");
   sprintf(buffer,"%s/%s.png", web_folder.c_str(), hist2D->GetName());
   canvas->SaveAs(buffer);
   hist2D = density_estimator_ptr_->getRootHist_yz();
   hist2D->Draw("CONTZ"); 
   sprintf(buffer,"%s/%s.png", web_folder.c_str(), hist2D->GetName());
   canvas->SaveAs(buffer);
   

  
   TLine * line = new TLine();
   TH1* hist = density_estimator_ptr_->getRootHist_z(z_bins,z_min,z_max, 0.1);
   sprintf(buffer,"%s/%s.png", web_folder.c_str(), hist->GetName());
   cout<<"Saving "<<buffer<<" . . ."<<endl;
   hist->SetMinimum(0.0);
   double y_max = hist->GetMaximum()*1.05;
   hist->SetMaximum(y_max);
   hist->Draw("C");
   line->DrawLine(tank_min,0,tank_min,hist->GetMaximum());
   line->DrawLine(tank_mid,0,tank_mid,hist->GetMaximum());
   line->DrawLine(tank_max,0,tank_max,hist->GetMaximum());
   canvas->SaveAs(buffer);

   hist->Write(); delete hist;
   hist = density_estimator_ptr_->getRootHist_z(z_bins,z_min,z_max, 1);

   sprintf(buffer,"%s/%s.png", web_folder.c_str(), hist->GetName());
   cout<<"Saving "<<buffer<<" . . ."<<endl;
   hist->SetMinimum(0.0);
   y_max = hist->GetMaximum()*1.05;
   hist->SetMaximum(y_max);  
   hist->Draw("C");
   line->DrawLine(tank_min,0,tank_min,hist->GetMaximum());
   line->DrawLine(tank_mid,0,tank_mid,hist->GetMaximum());
   line->DrawLine(tank_max,0,tank_max,hist->GetMaximum());
   canvas->SaveAs(buffer);  
   hist->Write(); delete hist;
   hist = density_estimator_ptr_->getRootHist_z(z_bins,z_min,z_max, 2);
   sprintf(buffer,"%s/%s.png", web_folder.c_str(), hist->GetName());
   cout<<"Saving "<<buffer<<" . . ."<<endl;
   hist->SetMinimum(0.0);
   y_max = hist->GetMaximum()*1.05;
   hist->SetMaximum(y_max);  
   hist->Draw("C");
   line->DrawLine(tank_min,0,tank_min,hist->GetMaximum());
   line->DrawLine(tank_mid,0,tank_mid,hist->GetMaximum());
   line->DrawLine(tank_max,0,tank_max,hist->GetMaximum());
   canvas->SaveAs(buffer);  
   hist->Write(); delete hist;

    printf("--- Storing data to root file: %s ---\n", output_root_file_path_.c_str());
    output_root_file_->Close();
    printf("-------- Saved State Histograms for %d iterations (%f)--------------\n", numIterationsPerformed, getRunningTime());
}

//////////////////////////////////////////////////////////////////////////////////////
string DCAKernelReconstructor::setup_output_folder_(){
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
    ss<<run_params_["OUTPUT_FOLDER_PATH"]<<"_DCA";

    string output_folder(ss.str());
    char command[501];
    //check to see if folder exists
    sprintf(command, "%s %s", run_params_["MAKE_DIR_COMMAND"].c_str(), run_params_["OUTPUT_FOLDER_PATH"].c_str());
    cout<<"Creating directory "<<run_params_["OUTPUT_FOLDER_PATH"]<<"."<<endl;
    printf("%s\n", command);

    int return_code = system(command);
    if(return_code != 0){
      printf("failed to make folder %s. \nDoes it exist? If not, there is a problem . . .\n", run_params_["OUTPUT_FOLDER_PATH"].c_str());
    }else{
      cout<<"Created directory "<< run_params_["OUTPUT_FOLDER_PATH"]<<". . . \n";
    }

    //copy the parametersFile and the data files to the new folder
    FileUtils::fileCopy(parameters_file_path_, output_folder, run_params_["COPY_COMMAND"]);

    return output_folder;
}//end of SetUpOutputFolder
//////////////////////////////////////////////////////////////////////////////////////

