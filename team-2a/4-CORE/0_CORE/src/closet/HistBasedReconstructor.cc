#include <sstream>
#include "omp.h"
#include <memory>

//ROOT Includes
#include "TH1F.h"
#include "TH3F.h"


//PromptGamma includes
#include "HistBasedReconstructor.h"
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
#include "utilities/Random.h"

using namespace prompt_gamma_reconstruction;
using namespace pg_tools;


HistBasedReconstructor::HistBasedReconstructor(const string &parameters_file_path):
    start_time_(TTimeStamp()),
    parameters_file_path_(parameters_file_path),
    run_time_parameters_(parameters_file_path),
    ptr_random_(RandomSingleton::Instance()),
    ptr_random_sqrt_(RandomSqrtSingleton::Instance()){

  event_file_path_ = run_time_parameters_["EVENT_FILE_PATH"];
  max_number_cones_ = std::atoi(run_time_parameters_["MAX_NUM_CONES"].c_str());
  density_points_per_cone_ = std::atoi(run_time_parameters_["DENSITY_POINTS_PER_CONE"].c_str());

  move_type_ = std::atoi(run_time_parameters_["MOVE_TYPE"].c_str());
  density_estimator_type_ = StringUtils::strtoi(run_time_parameters_["DENSITY_ESTIMATOR_TYPE"]);

  setup_density_matrix_();

  output_folder_path_ =  setup_output_folder_(run_time_parameters_["OUTPUT_FOLDER_PATH"]);
  output_root_file_path_ = output_folder_path_ ;
  output_root_file_path_ +=  run_time_parameters_["FOLDER_DELIMITER"];
  output_root_file_path_ += run_time_parameters_["RUN_NAME"] + ".root"; //set in _openOutputFile


  load_events_();
  number_cones_ = this->conic_sections_.size();

  printf("processing %d cones . . .\n",  number_cones_);

  int totalVoxels = bins_x_ * bins_y_ * bins_z_;
  double eventsPerVox = (double)number_cones_ / (double)totalVoxels;
  printf("--- Number of Voxels: %d ---\n", totalVoxels);
  printf("--- Events per Voxel: %.2f ---\n", eventsPerVox);
  printf("--- Run Time: %f sec ---\n", getRunningTime());

  cout<<"completed construction of HistBasedReconstructor . . ."<<endl;
}


HistBasedReconstructor::~HistBasedReconstructor() {

  if(output_root_file_){
    output_root_file_->Write();
    output_root_file_->Close();
    delete output_root_file_;
  }

  cout<<"completed destruction of HistBasedReconstructor . . ."<<endl;
}


void HistBasedReconstructor::setup_density_matrix_(){

  x_length_ = strtod(run_time_parameters_["X_LENGTH"].c_str(), 0);
  bins_x_ = StringUtils::strtoi(run_time_parameters_["X_BINS"].c_str());

  phantom_volume_->x_max = x_length_ / 2.0;
  phantom_volume_->x_min = -1.0 * phantom_volume_->x_max;

  y_length_ = strtod(run_time_parameters_["Y_LENGTH"].c_str(), 0);
  bins_y_ = StringUtils::strtoi(run_time_parameters_["Y_BINS"].c_str());

  phantom_volume_->y_max = y_length_ / 2.0;
  phantom_volume_->y_min = -1.0 * phantom_volume_->y_max;

  z_length_ = strtod(run_time_parameters_["Z_LENGTH"].c_str(), 0);
  bins_z_ = StringUtils::strtoi(run_time_parameters_["Z_BINS"].c_str());

  phantom_volume_->z_max = z_length_ / 2.0;
  phantom_volume_->z_min = -1.0 * phantom_volume_->z_max;

  density_matrix1_ptr_ = shared_ptr<DensityMatrix>(
        new DensityMatrix(
          phantom_volume_->x_min, phantom_volume_->x_max, bins_x_,
          phantom_volume_->y_min, phantom_volume_->y_max, bins_y_,
          phantom_volume_->z_min, phantom_volume_->z_max, bins_z_) );
  density_matrix2_ptr_ = shared_ptr<DensityMatrix>(
        new DensityMatrix(
          phantom_volume_->x_min, phantom_volume_->x_max, bins_x_,
          phantom_volume_->y_min, phantom_volume_->y_max, bins_y_,
          phantom_volume_->z_min, phantom_volume_->z_max, bins_z_) );
  density_matrix1_ptr_->print();
}

void HistBasedReconstructor::load_events_(){
  shared_ptr<EventsLoader> eventLoaderPtr;
  int data_file_format = atoi(run_time_parameters_["DATA_FILE_FORMAT"].c_str());

  string gamma_tree_name = run_time_parameters_["GAMMA_TREE_NAME"];
  double percent_cone_phantom_overlap = strtod(run_time_parameters_["CONE_PHANTOM_OVERLAP_PERCENTAGE"].c_str(), 0);
  int number_tries_per_random_point = static_cast<int>(100 * (1.0/percent_cone_phantom_overlap));

  eventLoaderPtr = EventsLoaderFactory::create(data_file_format, event_file_path_, &run_time_parameters_, phantom_volume_);
  eventLoaderPtr->LoadEvents(conic_sections_, number_tries_per_random_point);
  number_cones_ = static_cast<int>(conic_sections_.size());

  if(max_number_cones_ < number_cones_){
    conic_sections_.erase(conic_sections_.begin() + max_number_cones_, conic_sections_.end());
    random_shuffle(conic_sections_.begin(), conic_sections_.end());
    cout<<"WARNING: Only using "<< max_number_cones_ <<" of " << conic_sections_.size() << " possible cones . . ." << endl;
    number_cones_ = max_number_cones_;
  }

}


void HistBasedReconstructor::populate_density_matrix_(const vector<shared_ptr<ConicSection> > &conicSections, DensityEstimator &density_estimator){

  density_estimator.clear();

  cout<<"populating density matrix with "<<density_points_per_cone_<<" points per cone."<<endl;
  int total_entries  = 0, total_misses = 0;
  int num_cones = conicSections.size();

  int NUM_THREADS = std::atoi(run_time_parameters_["NUMBER_OF_THREADS"].c_str());
  #pragma omp parallel for reduction(+:total_entries, total_misses), num_threads(NUM_THREADS)
  for(int iCone=0; iCone < num_cones; ++iCone){
    vector<PGVector3> points_in_phantom;
    PGVector3 random_point;
    for(int i=0; i< density_points_per_cone_; i++){
        if( -1 != conicSections[iCone]->getRandomPointInPhantom(random_point, 1)){ //each attempt only gets 1 try to avoid creating local maxima
            points_in_phantom.push_back(random_point);
            ++total_entries;
        }else{
            ++total_misses;
        }
    }
    std::vector<PGVector3>::iterator iter = points_in_phantom.begin();
    float weight = 1.0/static_cast<float>(points_in_phantom.size());
    for(/* */; iter != points_in_phantom.end(); ++iter){
        #pragma omp critical
        { density_estimator.fill(*iter, weight); }
    }

    points_in_phantom.clear();
    if( (total_entries + total_misses) % (density_points_per_cone_*conicSections.size()/100) == 0 ) {
        cout<< "DENSITY_MATRIX "<< getRunningTime()<< "(s):"<< total_entries <<" entries, "<< total_misses<< " misses."<<endl;
    }
  }
  cout<<"DENSITY_MATRIX completed "<<getRunningTime()<<"(s):"
      << total_entries <<" entries, "<< total_misses<< " misses."<<endl;
}


void HistBasedReconstructor::repopulate_density_matrix_parallel(
              const vector<shared_ptr<ConicSection> > &conicSections,
              const DensityMatrix &old_de,
              DensityMatrix &new_de,
              int number_threads){
  float step_size = strtod(run_time_parameters_["RANDOM_STEP_SIZE"].c_str(),0);
  vector<DensityMatrix> de_array(number_threads);
  int number_cones = conicSections.size();
  int remaining_cones = number_cones;
  size_t range_min_index = 0;
  size_t range_max_index = 0;
  int number_cones_for_thread;

  vector<pair<size_t,size_t> > ranges(number_threads);
  for(int i=0; i<number_threads; ++i){
      number_cones_for_thread = ceil(float(remaining_cones)/float(number_threads - i));
      range_min_index = range_max_index;
      range_max_index = range_max_index + number_cones_for_thread;
      remaining_cones -= number_cones_for_thread;
      ranges[i] = pair<size_t,size_t>(range_min_index,range_max_index);
      cout<<"range, "<<range_min_index<<", "<<range_max_index<<endl;
  }

  #pragma omp parallel num_threads(number_threads)
  {
      int thread_id = omp_get_thread_num();
      cout<<"thread id: "<<thread_id<<endl;
      random_walk_algo(conicSections, ranges[thread_id].first, ranges[thread_id].second, old_de, de_array[thread_id], step_size);
  }

  for(int i=0; i<number_threads; ++i){
    new_de += de_array[i];
  }

}


void HistBasedReconstructor::random_walk_algo(const vector<shared_ptr<ConicSection> > &conicSections,
        size_t first_cone_index, size_t last_cone_index, const DensityMatrix &old_de, DensityMatrix &new_de, float step_size){

    TTimeStamp time;
    pg_tools::Random rand(last_cone_index + time.GetNanoSec());
    new_de = old_de;
    new_de.clear();
    cout<<"populating density matrix for cones["<<first_cone_index<<","<<last_cone_index<< "], points/cone = "<< density_points_per_cone_<<endl;

    int total_entries  = 0, total_misses = 0, density_zero=0;
    double gaussian_width = strtod(run_time_parameters_["GAUSS_WIDTH"].c_str(), 0);
    double offset_x = strtod(run_time_parameters_["OFFSET_X"].c_str(), 0);
    double offset_y = strtod(run_time_parameters_["OFFSET_Y"].c_str(), 0);
    double C1 = -0.5/(gaussian_width*gaussian_width);

    for(size_t jCone=first_cone_index; jCone < last_cone_index; ++jCone){

    PGVector3 current_point, next_random_point;
    vector< std::pair<PGVector3,float> > points_in_phantom;
    float weight = 0.0, current_weight = 0.0, sum_of_weights=0.0;
    float weight_max;
    PGVector3 most_dense_point;

    int number_tries_for_random = 0;

    weight_max = 0.0;

    current_point = conicSections[jCone]->getLikelyOrigin();
    current_weight = old_de.getDensity(current_point);
    int i;
    for(i=0; i< density_points_per_cone_; i++){
      number_tries_for_random = conicSections[jCone]->getRandomStepInPhantom(current_point, next_random_point, step_size, 100);

      if(-1 != number_tries_for_random){ //each attempt only gets 1 try to avoid creating local maxima
        double prior = exp(C1*((offset_x - next_random_point.x)*(offset_x - next_random_point.x)+(offset_y - next_random_point.y)*(offset_y - next_random_point.y)));
        weight = old_de.getDensity(next_random_point)*prior;
        if(weight_max < weight){ most_dense_point = next_random_point; }
        if( weight != weight) throw std::runtime_error("HistBasedReconstructor::repopulate_density_matrix_()\nweight is NAN!");

        float random_num = rand.Rndm();
        if( weight > random_num* current_weight){
          /// DSM - to soften the convergence, move the following line outside of the if block
          ///  so that all points will be included in the histogram even if the current point is not moved
          ///  to the location.
          auto pair = std::pair<PGVector3, float>(next_random_point, weight);
          points_in_phantom.push_back(pair);
          current_point = next_random_point;
          current_weight = weight;
          sum_of_weights += weight;
          ++total_entries;
        }
      }else{
        ++total_misses;
      }
    }

    conicSections[jCone]->setLikelyOrigin(current_point);

    std::vector<std::pair<PGVector3,float> >::iterator iter = points_in_phantom.begin();

    //     float weight = 1.0/static_cast<float>(points_in_phantom.size());
    ///weight all the points equally when the points are generated according to Markov chain

    if(0.0 < sum_of_weights){

      conicSections[jCone]->setLikelyOrigin(most_dense_point);
      float cone_weight = conicSections[jCone]->getWeight();
      for(/* */; iter != points_in_phantom.end(); ++iter){
        weight = (*iter).second/sum_of_weights;
        if( sum_of_weights != sum_of_weights) throw std::runtime_error("HistBasedReconstructor::repopulate_density_matrix_()\nsum_of_weights is NAN!");
        if(weight != weight) {
          cout<<(*iter).first.print()<<", "<<"weight, "<<weight<<", sum_of_weights, "<<sum_of_weights<<endl;
          throw std::runtime_error("HistBasedReconstructor::repopulate_density_matrix_()\npoint weight is NAN!");
        }
        new_de.fill((*iter).first, weight*cone_weight);
      }
    }else{
      ++density_zero;
    }

    points_in_phantom.clear();

    }
    cout<<"DENSITY_MATRIX completed "<<getRunningTime()<<"(s):"<< total_entries<<" entries, "
                                   << total_misses << " misses, "<< density_zero << " zero contributions."<<endl;

}

void HistBasedReconstructor::repopulate_density_matrix_(const vector<shared_ptr<ConicSection> > &conicSections,
							const DensityMatrix &old_de,
							DensityMatrix &new_de){
    new_de.clear();

    cout<<"populating density matrix with "<<density_points_per_cone_<<" points per cone."<<endl;


    int num_conic_sections = conicSections.size();
    long long random_index;
  
    int total_entries  = 0, total_misses = 0, density_zero=0;
    int NUM_THREADS = std::atoi(run_time_parameters_["NUMBER_OF_THREADS"].c_str());
    float step_size = strtod(run_time_parameters_["RANDOM_STEP_SIZE"].c_str(),0);

    DensityMatrix de_private(new_de);
    double sum_of_weights = 0.0;
    #pragma omp parallel for private(sum_of_weights, de_private, random_index),\
      reduction(+:total_entries, total_misses, density_zero), num_threads(NUM_THREADS)
    for(int jCone=0; jCone < num_conic_sections; ++jCone){
        printf("generating points for cone %d . . .\n", jCone);

        PGVector3 current_point, next_random_point;
        vector< std::pair<PGVector3,float> > points_in_phantom;
        float weight = 0.0, current_weight = 0.0;
        float weight_max;
        PGVector3 most_dense_point;

        int number_tries_for_random = 0;

        weight_max = 0.0;

        current_point = conicSections[jCone]->getLikelyOrigin();
        current_weight = old_de.getDensity(next_random_point);
        int i;
        for(i=0; i< density_points_per_cone_; i++){
            number_tries_for_random = conicSections[jCone]->getRandomStepInPhantom(current_point, next_random_point, step_size, 100);

            if(-1 != number_tries_for_random){ //each attempt only gets 1 try to avoid creating local maxima
                weight = old_de.getDensity(next_random_point);
                if(weight_max < weight){ most_dense_point = next_random_point; }
                if( weight != weight) throw std::runtime_error("HistBasedReconstructor::repopulate_density_matrix_()\nweight is NAN!");

                random_index = static_cast<long long>(i)*jCone + i;
                if( weight > ptr_random_->Instance()->getRandIndex(random_index) * current_weight){
					auto pair = std::pair<PGVector3, float>(next_random_point, weight);
					points_in_phantom.push_back(pair);
					current_point = next_random_point;
					current_weight = weight;
					sum_of_weights += weight;
					++total_entries;
                }
            }else{
                ++total_misses;
            }
        }

        conicSections[jCone]->setLikelyOrigin(current_point);

        std::vector<std::pair<PGVector3,float> >::iterator iter = points_in_phantom.begin();
        
        if(0.0 < sum_of_weights){
            de_private = new_de;
            de_private.clear();
            conicSections[jCone]->setLikelyOrigin(most_dense_point);
            for(/* */; iter != points_in_phantom.end(); ++iter){
                weight = (*iter).second/sum_of_weights;
                if( sum_of_weights != sum_of_weights) throw std::runtime_error("HistBasedReconstructor::repopulate_density_matrix_()\nsum_of_weights is NAN!");
                if(weight != weight) {
                    cout<<(*iter).first.print()<<", "<<"weight, "<<weight<<", sum_of_weights, "<<sum_of_weights<<endl;
                    throw std::runtime_error("HistBasedReconstructor::repopulate_density_matrix_()\npoint weight is NAN!");
                }
                de_private.fill((*iter).first, weight);
            }
            #pragma omp critical
            {new_de += de_private;}
        }else{
            ++density_zero;
        }

        points_in_phantom.clear();
        if( (total_entries + total_misses) % 100000000 == 0 ) {
            cout<< "DENSITY_MATRIX "<< getRunningTime()<< "(s):"<< total_entries <<" entries, "<< total_misses<< " misses."<<endl;
        }
    }
    cout<<"DENSITY_MATRIX completed "<<getRunningTime()<<"(s):"<< total_entries<<" entries, "
        << total_misses << " misses, "<< density_zero << " zero contributions."<<endl;
}

// Draws function graphics in randomly choosen interval
void HistBasedReconstructor::StartExperiment() {

  int iterations = atoi(run_time_parameters_["ITERATIONS"].c_str());
  int NUM_THREADS = std::atoi(run_time_parameters_["NUMBER_OF_THREADS"].c_str());
  populate_density_matrix_(conic_sections_,*density_matrix1_ptr_);
  save_density_info_(0);

  for(int i=1; i<iterations+1; ++i){
    cout<<"Starting iteration " << i << " . . ."<<endl;
    repopulate_density_matrix_parallel(conic_sections_,*density_matrix1_ptr_,*density_matrix2_ptr_, NUM_THREADS);

    density_matrix1_ptr_.swap(density_matrix2_ptr_);
    if ( ( i < 10 ) ||
         ( i < 100 &&  i % 10 == 0) ||
         ( i < 1000 &&  i % 100 == 0) ||
         ( i < 10000 &&  i % 1000 == 0) ||
         ( i < 25000 &&  i % 5000 == 0) ||
         ( i < 100000 &&  i % 10000 == 0)  ){
        save_density_info_(i);
    }
  }
}

void HistBasedReconstructor::open_output_file_( int number_iterations) {

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

void HistBasedReconstructor::save_density_info_( int numIterationsPerformed) {
  open_output_file_( numIterationsPerformed);
  assert(output_root_file_->IsOpen());
  output_root_file_->cd();

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
                    bins_z_, phantom_volume_->z_min, phantom_volume_->z_max,
                    bins_x_, phantom_volume_->z_min, phantom_volume_->z_max,
                    bins_y_, phantom_volume_->y_min, phantom_volume_->y_max) ;
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

  //save density matrix
  TH3* hist3D = density_matrix1_ptr_->getRootHist();
  hist3D->SetName("density3D");
  hist3D->Write();
  TH1* hist = hist3D->Project3D("x");
  hist->SetName("density_x");
  hist->Write(); delete hist;
  hist = hist3D->Project3D("y");
  hist->SetName("density_y");
  hist->Write(); delete hist;
  hist = hist3D->Project3D("z");
  hist->SetName("density_z");
  hist->Write(); delete hist;
  hist = hist3D->Project3D("xy");
  hist->SetName("density_xy");
  hist->Write(); delete hist;
  hist = hist3D->Project3D("yz");
  hist->SetName("density_yz");
  hist->Write(); delete hist;
  hist = hist3D->Project3D("xz");
  hist->SetName("density_xz");
  hist->Write(); delete hist;
  delete hist3D;

  printf("--- Storing data to root file: %s ---\n", output_root_file_path_.c_str());
  output_root_file_->Close();
  printf("-------- Saved State Histograms for %d iterations (%f)--------------\n", numIterationsPerformed, getRunningTime());
}



void HistBasedReconstructor::choose_random_points_() {

    cout<<"WARNING! method HistBasedReconstructor::choose_random_points_ no longer supported."<<endl;

    // Variable Declaration
    int changeCount = 0;

    PGVector3 current_most_dense_point, new_pos, random_point;
    double best_density, new_density;

    for(int jConeNum=0; jConeNum < number_cones_; ++jConeNum){
        current_most_dense_point = conic_sections_[jConeNum]->getLikelyOrigin();
        best_density= static_cast<double>(density_matrix1_ptr_->getDensity(current_most_dense_point));
        for (int i=0; i < test_points_per_cone_; i++) {

          int number_tries = conic_sections_[jConeNum]->getRandomPointInPhantom(random_point,1000);
          if( -1 == number_tries) continue; //did not find random point so continue to next point

          new_density = static_cast<double>(density_matrix1_ptr_->getDensity(random_point));

          //double ratio = new_density/(new_density + best_density);
          if( new_density > best_density ||
              (2 == move_type_ && new_density > ptr_random_->getRand()*(new_density + best_density)) ||
              (3 == move_type_ && new_density > ptr_random_->getRand()*best_density ) ){

            conic_sections_[jConeNum]->setLikelyOrigin(random_point);
            best_density = new_density;
            changeCount++;
          }
        }
        if(jConeNum % 10000 == 0) {
            printf("--- Test Points per cone: %d, time %.1f, Number of Position Changes: %d, ratio: %.3f\n",
              test_points_per_cone_, getRunningTime(), changeCount, (double)changeCount/(double)(test_points_per_cone_ * number_cones_));
        }
    }
    save_density_info_(test_points_per_cone_);
    printf("--- Test Points per cone: %d, time %.1f, Number of Position Changes: %d, ratio: %.3f\n",
           test_points_per_cone_, getRunningTime(), changeCount, (double)changeCount/(double)(test_points_per_cone_ * number_cones_));
}



//////////////////////////////////////////////////////////////////////////////////////
string HistBasedReconstructor::setup_output_folder_(const string &outputFolderPath){
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
    ss<<run_time_parameters_["OUTPUT_FOLDER_PATH"]<<"_"<<density_points_per_cone_;
    ss<<"_"<<test_points_per_cone_;
    ss<<"_"<<bins_x_<<"x"<<bins_y_<<"x"<<bins_z_;
    ss<<"_move" << move_type_;
    ss<<"_densityEstimator" << density_estimator_type_;

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

