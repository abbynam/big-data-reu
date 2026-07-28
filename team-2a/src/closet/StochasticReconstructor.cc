//C++ includes
#include <memory>
//
////ROOT includes
//#include "TTimeStamp.h"
//#include "TFile.h"

#include "StochasticReconstructor.h"
#include "DensityMatrix.h"
#include "ASHDensity.h"
#include "EventsLoader.h"
#include "PGSimulatedEventsLoader.h"
#include "SimulatedCCEventsLoader.h"
#include "CSVEventsLoader.h"
#include "DetectorEffectsEventsLoader.h"
#include "PhantomVolumeBuilder.h"
#include "utilities/utilities.h"
#include "EventsLoaderFactory.h"
#include "utilities/FileUtils.h"

using namespace prompt_gamma_reconstruction;
using namespace pg_tools;


StochasticReconstructor::StochasticReconstructor(const string &parameters_file_path):
    parameters_file_path_(parameters_file_path),
    run_time_parameters_(parameters_file_path),
    num_scatters_(0){

    auto phantom_volume_ptr = PhantomVolumeBuilder::build(run_time_parameters_);
    phantom_volume_ = *phantom_volume_ptr;

    iterations_ = run_time_parameters_.get_int("ITERATIONS");
    double percent_cone_phantom_overlap = run_time_parameters_.get_double("CONE_PHANTOM_OVERLAP_PERCENTAGE");
    
    number_tries_for_random_ = static_cast<int>(100 * (1.0/percent_cone_phantom_overlap));

    event_file_path_ = run_time_parameters_["EVENT_FILE_PATH"];//path to file with gamma dections
    max_number_cones_ = run_time_parameters_.get_int("MAX_NUM_CONES"); //limit the number of cones used in reconstruction

    output_folder_path_ =  setup_output_folder_(run_time_parameters_["OUTPUT_FOLDER_PATH"]);
    output_root_file_path_ = output_folder_path_ ;
    output_root_file_path_ +=  run_time_parameters_["FOLDER_DELIMITER"];
    output_root_file_path_ += run_time_parameters_["RUN_NAME"] + ".root"; //set in _openOutputFile

    // Set starting timestamp
    start_time_ = new TTimeStamp();

    // Copy initial cone data into modifed cone data variables
    int totalVoxels = x_bins_ * y_bins_ * z_bins_;
    double eventsPerVox = (double)num_cones_ / (double)totalVoxels;
    printf("--- Number of Voxels: %d ---\n", totalVoxels);

    printf("--- Run Time: %f sec ---\n", get_running_time());
    setup_density_matrix_(1.0);
    load_events_();
    num_cones_ = this->conic_sections_.size();
    printf("--- Events per Voxel: %.2f ---\n", eventsPerVox);
    printf("processing %d cones . . .\n", num_cones_);
}


StochasticReconstructor::~StochasticReconstructor() {
    //ROOT file should be explicitly deleted
    if(output_root_file_){
      output_root_file_->Write();
      output_root_file_->Close();
      delete output_root_file_;
    }
    cout<<"completed destruction of StochasticReconstructor . . ."<<endl;
}


void StochasticReconstructor::load_events_(){
    
    shared_ptr<EventsLoader> eventsLoaderPtr;
    const int data_file_format = run_time_parameters_.get_int("DATA_FILE_FORMAT");

    eventsLoaderPtr = EventsLoaderFactory::create(data_file_format, event_file_path_, &run_time_parameters_,
                                                  std::make_shared<const PhantomVolume>(phantom_volume_));
    
    vector< shared_ptr<ConicSection>> cs;
    eventsLoaderPtr->LoadEvents(cs, number_tries_for_random_);
    for(size_t i = 0; i < cs.size(); ++i){
        conic_sections_.push_back(*cs[i]);
    }

    num_scatters_ = static_cast<int>(conic_sections_.size());
    cout<<"StochasticReconstructor::LoadEvents: processing "<<num_scatters_<<" events . . ."<<endl;

    if(max_number_cones_ > num_scatters_){
        cout<<"WARNING: Only using "<< num_scatters_ <<" cones; requested " << max_number_cones_ << " . . ." << endl;
    }
}


void StochasticReconstructor::setup_density_matrix_(float bins_scalar){

    x_length_ = run_time_parameters_.get_double("X_LENGTH");  
    x_bins_   = run_time_parameters_.get_int("X_BINS");

    /// multiply bins by scalar so that the binning can be changed based on the number of iterations.
    /// usually we make bins smaller after many iterations so that we can continually
    /// improve the resolution -- DSM 2012-04-08
    x_bins_ = x_bins_ * bins_scalar;
    if(x_bins_ % 2 == 0) x_bins_++; //use odd number of bins to avoid asymmetry around 0

//    phantom_volume_ptr_->x_max = x_length_ / 2.0;
//    phantom_volume_ptr_->x_min = -1.0 * phantom_volume_ptr_->x_max;

    y_length_ = run_time_parameters_.get_double("Y_LENGTH");
    y_bins_ = run_time_parameters_.get_int("Y_BINS");
    y_bins_ = floor(y_bins_ * bins_scalar);
    
    if(y_bins_ % 2 == 0) y_bins_++; //use odd number of bins to avoid asymmetry around 0

//    phantom_volume_ptr_->y_max = y_length_ / 2.0;
//    phantom_volume_ptr_->y_min = -1.0 * phantom_volume_ptr_->y_max;

    z_length_ = run_time_parameters_.get_double("Z_LENGTH");
    z_bins_ = run_time_parameters_.get_int("Z_BINS");
    z_bins_  = floor(z_bins_ * bins_scalar);
    if(z_bins_ % 2 == 0) z_bins_++; //use odd number of bins to avoid assymmetry around 0

//    phantom_volume_ptr_->z_max = z_length_ / 2.0;
//    phantom_volume_ptr_->z_min = -1.0 * phantom_volume_ptr_->z_max;

    //Create the density estimator based on the DENSITY_ESTIMATOR_TYPE specified
    // in the parameters file. 
    
    /// @TODO Break this code out into a factory

    int density_estimator_type = run_time_parameters_.get_int("DENSITY_ESTIMATOR_TYPE");
    if(2 == density_estimator_type){//using averaged shifted histograms
        cout<<"Using Averaged Shifted Histograms for density estimation . . ."<<endl;
        int number_of_shifts = run_time_parameters_.get_int("NUMBER_OF_SHIFTS");
        density_estimator_ptr_ = shared_ptr<DensityEstimator>(
            new ASHDensity(number_of_shifts,
                           phantom_volume_.x_min, phantom_volume_.x_max, x_bins_,
                           phantom_volume_.y_min, phantom_volume_.y_max, y_bins_,
                           phantom_volume_.z_min, phantom_volume_.z_max, z_bins_) );
    }else{
        cout<<"Using a standard 3D histogram for density estimation . . ."<<endl;
        density_estimator_ptr_ = shared_ptr<DensityEstimator>(
          new DensityMatrix(
                  phantom_volume_.x_min, phantom_volume_.x_max, x_bins_,
                  phantom_volume_.y_min, phantom_volume_.y_max, y_bins_,
                  phantom_volume_.z_min, phantom_volume_.z_max, z_bins_) );
    }

    //phantom_volume_.print();
    //density_estimator_ptr_->print();
} 

/*! \brief Load the Compton Scatter representative points into the density matrix
 */
void StochasticReconstructor::populate_density_matrix_(const vector<ConicSection> &conicSections, DensityEstimator &density_estimator){
  density_estimator.clear();
  density_estimator = 0.0;
  
  //vector<ConicSection>::iterator conic_iter = conicSections.begin();
  auto conic_iter = conicSections.begin();
  for(/* */; conic_iter != conicSections.end(); ++conic_iter){
      auto origin = conic_iter->getLikelyOrigin();
      auto weight = conic_iter->getWeight();
	  density_estimator.fill(origin, weight);
  }
}

/*! \brief Controls the flow of the reconstruction algorithm
 */
void StochasticReconstructor::start_calc() {

    /// @TODO Create bin_width and temperatures run time parameters 
    
    //Hardcoded options for testing multiple bin widths and annealing 
    float bin_width_scalars[] = {1.0}; 
    float temperatures[] = {1.0};

    size_t num_test_points = static_cast<size_t>(sizeof(bin_width_scalars)/sizeof(bin_width_scalars[0]));
    
    for(size_t iBin_width_scalar=0; iBin_width_scalar < num_test_points; ++iBin_width_scalar){
        setup_density_matrix_(bin_width_scalars[iBin_width_scalar]);
        populate_density_matrix_(conic_sections_, *density_estimator_ptr_);
        printf("completed construction of StochasticReconstructor width scalar %.1f (%0.3f sec) . . .\n", bin_width_scalars[iBin_width_scalar], get_running_time());

        calculate_density(iterations_, iBin_width_scalar, temperatures[iBin_width_scalar]);
    }
}

/*! \brief Stochastic origins ensembles implementation function.
 * Loops over the cones attempting to improve the representative 
 * point for each cone.
 */
void StochasticReconstructor::calculate_density( int iterations, int width_num, float temperature ) {

    // Variable Declaration
    int changeCount = 0;

    PGVector3 oldEventPos, randomPosition;
    double oldDensity, newDensity;

    // Loop through iterations
    num_cones_ = conic_sections_.size();
    printf("--- Number of MCMC iterations: %d ---\n", iterations);
    printf("--- Number of events: %d ---\n", num_cones_);

    save_density_info_(width_num*iterations, true);

    int NUM_THREADS = run_time_parameters_.get_int("NUMBER_OF_THREADS");
    int i = 0;
    double old_prior = 1.0; //Gaussian likelihood weight of current representative point
    double new_prior = 1.0; //Gaussian likelihood weight of random test point
    double event_weight = 0.0; //Better (logistic regression) events weight more

    double gaussian_width = run_time_parameters_.get_double("GAUSS_WIDTH");
    double offset_x = run_time_parameters_.get_double("OFFSET_X");
    double offset_y = run_time_parameters_.get_double("OFFSET_Y");
    double C1 = -0.5/(gaussian_width*gaussian_width);

    for (; i < iterations; i++) {
        changeCount = 0;

        long long index = i;
        #pragma omp parallel for reduction(+:changeCount), \
                private(index, oldEventPos, randomPosition, oldDensity, newDensity, event_weight, old_prior, new_prior ) num_threads(NUM_THREADS)
        for(int jEventNum=0; jEventNum< num_cones_; ++jEventNum){

            oldEventPos = conic_sections_[jEventNum].getLikelyOrigin();
            old_prior = exp(C1*( (offset_x - oldEventPos.x)*(offset_x-oldEventPos.x)+(offset_y - oldEventPos.y)*(offset_y - oldEventPos.y)));

            event_weight = conic_sections_[jEventNum].getWeight();
            oldDensity = (density_estimator_ptr_->getDensity(oldEventPos) - event_weight)*old_prior;

            int number_tries = conic_sections_[jEventNum].getRandomPointInPhantom(randomPosition, number_tries_for_random_);
            if( -1 == number_tries) continue; //did not find random point so continue to next point
            new_prior = exp(C1*((offset_x - randomPosition.x)*(offset_x - randomPosition.x)+(offset_y - randomPosition.y)*(offset_y - randomPosition.y)));

            newDensity = density_estimator_ptr_->getDensity(randomPosition)*new_prior;

            index = static_cast<long long>(num_cones_)*i + jEventNum;
            float rand = RandomSingleton::Instance()->getRandIndex(index);

            //Acquire the mutex and make the change if the new density
            // is greater than the old density times a random raised to 
            // power temperature.
            if( newDensity >= pow(rand, temperature)*( oldDensity )){
                conic_sections_[jEventNum].setLikelyOrigin(randomPosition);
                #pragma omp critical
                {
                    density_estimator_ptr_->updateMatrix(oldEventPos,randomPosition,event_weight);
                }
                ++changeCount;
            }
        }//end for

        //Give status update to the terminal
        if ( ( i < 1000 &&  (i+1) % 100 == 0) ||
             ( i < 10000 &&  (i+1) % 1000 == 0) ||
             ( i < 25000 &&  (i+1) % 5000 == 0) ||
             ( i < 100000 &&  (i+1) % 10000 == 0)  ){
            printf("Iteration: %d, time %.1f, Number of Position Changes: %d, ratio: %.3f",
                    i+1, get_running_time(), changeCount, static_cast<double>(changeCount)/static_cast<double>(num_cones_) );
            cout<<endl;//print it out now!

            //Save data for analysis
            save_density_info_(width_num*iterations + i+1, true);
        }
    }
    save_density_info_(width_num*iterations+iterations, true);
    printf("--- Total Iterations: %d, time %.1f, Number of Position Changes: %d, %.1f, ratio: %.3f\n",
                iterations, get_running_time(), changeCount, (double)(num_cones_), (double)changeCount/(double)(num_cones_));
}


void StochasticReconstructor::open_output_file_( int number_iterations) {

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
    
    sprintf(str_buffer_,"%s%s%s_%s_%d.root",
            output_folder_path_.c_str(),
            run_time_parameters_["FOLDER_DELIMITER"].c_str(),
            run_time_parameters_["RUN_NAME"].c_str(),
            fileNameRoot.c_str(),
            number_iterations);

    output_root_file_ = new TFile(str_buffer_,"RECREATE");
    
    if(! output_root_file_->IsOpen()){
        string error_message("StochasticReconstructor::open_output_file_: failed to open ");
        error_message += str_buffer_;
        error_message += ".\nABORTING\n\n";
        throw runtime_error(error_message);
    };
};


string StochasticReconstructor::setup_output_folder_(const string &outputFolderPath){
//////////////////////////////////////////////////////////////////////////////////////
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

    string output_folder = run_time_parameters_["OUTPUT_FOLDER_PATH"];
    output_folder += "_";
    output_folder += run_time_parameters_["ITERATIONS"];
    output_folder += "_estimator";
    output_folder += run_time_parameters_["DENSITY_ESTIMATOR_TYPE"];
    output_folder += "_move";
    output_folder += run_time_parameters_["MOVE_TYPE"];
    char command[501];


    //check to see if folder exists
    sprintf(command, "%s %s",
            run_time_parameters_["MAKE_DIR_COMMAND"].c_str(),
            output_folder.c_str());

    printf("%s\n", command);

    if(system(command) != 0){
      printf("failed to make %s. \nDoes it exist? If not, there is a problem . . .\n",
        output_folder.c_str());
    }else{
      cout<<"Created directory "<< output_folder<<". . . \n";
    }

    //copy the parametersFile and the data files to the new folder
    FileUtils::fileCopy(parameters_file_path_, output_folder, run_time_parameters_["COPY_COMMAND"]);

    return output_folder;
}//end of SetUpOutputFolder


void StochasticReconstructor::save_density_info_( int numIterationsPerformed, Bool_t saveEventData) {
    //////// Store start data and create image for GUI and to store to file
    open_output_file_( numIterationsPerformed);
    assert(output_root_file_->IsOpen());
    output_root_file_->cd();

    /// Declare ntuples to store cone and final event data
    sprintf(str_buffer_,"eventData_iter_%05d", numIterationsPerformed);

    TTree output_tree("promptGammaOutput", "Gamma Reconstruction Information");

    float final_point[3], origin_true[3];
    float p1[3], p2[3], p3[3], p1_true[3], p2_true[3], p3_true[3];

    float incident_energy[3], incident_energy_true[3], energy_lost;
    float energy_deposit[3], energy_deposit_true[3];

    float  apex[3], axis[3];
    float angle[2], angle_true[2];
    float dca=0.0, dca_weight=0.0;
    PGVector3 pca;
    float point_closest_approach[3];

    output_tree.Branch("incident_energy", &incident_energy, "E0/F:E1/F:E2/F");
    output_tree.Branch("final_point", &final_point, "final_x/F:final_y/F:final_z/F");
    output_tree.Branch("angle", &angle, "angle_1/F:angle_2/F");
    output_tree.Branch("angle_true", &angle_true, "angle_1_true/F:angle_2_true/F");
    output_tree.Branch("distance_closest_approach", &dca, "dca/F");
    output_tree.Branch("dca_weight", &dca_weight, "dca_weight/F");
    output_tree.Branch("point_closest_approach", &point_closest_approach, "pca_x/F:pca_y/F:pca_z/F");
    
    //Write the initial state data only if this is the first iteration.
    // the initial information does not need to saved more than once.
    if( 0 == numIterationsPerformed){

        // Set the ROOT tree branch addresses to the addresses of the 
        // appropriate array.
        
        ///POSITION
        output_tree.Branch("origin", &origin_true, "origin_x/F:origin_y/F:origin_z/F");
        output_tree.Branch("position1", &p1, "p1_x/F:p1_y/F:p1_z/F");
        output_tree.Branch("position2", &p2, "p2_x/F:p2_y/F:p2_z/F");
        output_tree.Branch("position3", &p3, "p3_x/F:p3_y/F:p3_z/F");
        output_tree.Branch("position1true", &p1_true, "p1_x_true/F:p1_y_true/F:p1_z_true/F");
        output_tree.Branch("position2true", &p2_true, "p2_x_true/F:p2_y_true/F:p2_z_true/F");
        output_tree.Branch("position3true", &p3_true, "p3_x_true/F:p3_y_true/F:p3_z_true/F");

        ///ENERGY
        output_tree.Branch("incident_energy", &incident_energy, "E0/F:E1/F:E2/F");
        output_tree.Branch("incident_energy_true", &incident_energy_true, "E0_true/F:E1_true/F:E2_true/F");
        output_tree.Branch("energy_deposit", &energy_deposit, "deposit_1/F:deposit_2/F:deposit_3/F:deposit_3/F");
        output_tree.Branch("energy_deposit_true", &energy_deposit_true, "deposit_1_true/F:deposit_2_true/F:deposit_true_3/F");
        output_tree.Branch("energy_lost", &energy_lost, "energy_lost/F");

        ///SCATTERING ANGLES
        output_tree.Branch("angle", &angle, "angle_1/F:angle_2/F");
        output_tree.Branch("angle_true", &angle_true, "angle_1_true/F:angle_2_true/F");

        ///CALCULATED VALUES
          output_tree.Branch("apex", &apex, "apex_x/F:apex_y/F:apex_z/F");
          output_tree.Branch("axis", &axis, "axis_x/F:axis_y/F:axis_z/F");
    };

    // Create matrices to store final density data
    sprintf(str_buffer_,"density_iter_%05d", numIterationsPerformed);
    density_hist_ = new TH3F(str_buffer_, "3D histogram of density matrix;z (mm);z (mm);y (mm)",
                             z_bins_, phantom_volume_.z_min, phantom_volume_.z_max,
                             z_bins_, phantom_volume_.z_min, phantom_volume_.z_max,
                             y_bins_, phantom_volume_.y_min, phantom_volume_.y_max);
    sprintf(str_buffer_,"Event Origin Density (%d iterations);x (mm);y (mm);z (mm)", numIterationsPerformed);
    density_hist_->SetTitle(str_buffer_);

    PGVector3 origin(0.0,0.0,0.0);
    
    // Loop through cones to store data (apexPos, scatVect, scatAng, finalEventPos) into ntuples
    
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

        PGVector3 origin_vec(origin_true[0], origin_true[1], origin_true[2]);
        dca = iter->getDistanceToPoint(origin_vec);
        dca_weight = iter->getWeight();

        pca = iter->getPointOfClosestApproach(origin_vec);
        pca.set_xyz_array(point_closest_approach);

        //Copy the data from the local variables to the ROOT tree branches
        output_tree.Fill();
    }

    printf("--- Storing data to root file: %s (%f)---\n", output_root_file_path_.c_str(), get_running_time());
    output_tree.Write();

    if(! saveEventData) {
        delete density_hist_;
    }

    output_root_file_->Close();
    printf("-------- Saved State Histograms for %d iterations (%f)--------------\n", numIterationsPerformed, get_running_time());
};

