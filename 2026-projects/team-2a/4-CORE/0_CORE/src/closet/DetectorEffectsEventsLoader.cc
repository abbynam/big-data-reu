#include "DetectorEffectsEventsLoader.h"

// Standard C++ Includes
#include <cmath>
#include <valarray>
#include <vector>
#include <stdexcept>
#include <iostream>


// ROOT Includes
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"

//PG Includes
#include "PhantomVolume.h"
#include "ReconstructionParabola.h"
#include "ReconstructionEllipse.h"
#include "ConicSection.h"
#include "EventsLoader.h"
#include "TripleScatter.h"


using namespace std;
using namespace prompt_gamma_reconstruction;

DetectorEffectsEventsLoader::DetectorEffectsEventsLoader(const string &data_file_path, const pg_tools::RunTimeParameters *params, shared_ptr<const PhantomVolume>phantom):
   EventsLoader(data_file_path, params, phantom){

    gamma_tree_name_ = (*params)["GAMMA_TREE_NAME"];

}


void DetectorEffectsEventsLoader::AddDetectorEffects(TripleScatter &ts) const{

  static bool have_set_parameters = false;
  static PGVector3 positionUnc1, positionUnc2, positionUnc3;

  if( ! have_set_parameters){

    positionUnc1.x =  strtod((*params_ptr_)["DETECTOR1_X_UNC"].c_str(),0);
    positionUnc1.y =  strtod((*params_ptr_)["DETECTOR1_Y_UNC"].c_str(),0);
    positionUnc1.z =  strtod((*params_ptr_)["DETECTOR1_Z_UNC"].c_str(),0);
    ts.applyPositionUncertainty(positionUnc1, 0);

    positionUnc2.x =  strtod((*params_ptr_)["DETECTOR2_X_UNC"].c_str(),0);
    positionUnc2.y =  strtod((*params_ptr_)["DETECTOR2_Y_UNC"].c_str(),0);
    positionUnc2.z =  strtod((*params_ptr_)["DETECTOR2_Z_UNC"].c_str(),0);
    ts.applyPositionUncertainty(positionUnc2, 1);

    positionUnc3.x =  strtod((*params_ptr_)["DETECTOR3_X_UNC"].c_str(),0);
    positionUnc3.y =  strtod((*params_ptr_)["DETECTOR3_Y_UNC"].c_str(),0);
    positionUnc3.z =  strtod((*params_ptr_)["DETECTOR3_Z_UNC"].c_str(),0);

    have_set_parameters = true;
  }
  ts.applyPositionUncertainty(positionUnc3, 2);

  //energy uncertainty uses the formula dE(E) = sqrt(alpha +beta*E)
  string detector_type = (*params_ptr_)["DETECTOR1_MATERIAL"];
  float e_scalar =  strtod((*params_ptr_)["DETECTOR1_E_SCALAR"].c_str(),0);

  ts.applyEnergyUncertainty(detector_type, e_scalar, 0);

  detector_type = (*params_ptr_)["DETECTOR2_MATERIAL"];
  e_scalar =  strtod((*params_ptr_)["DETECTOR2_E_SCALAR"].c_str(),0);
  ts.applyEnergyUncertainty(detector_type, e_scalar, 1);

}


void DetectorEffectsEventsLoader::LoadEvents(vector<shared_ptr<ConicSection> > &conics, int number_tries_per_random_point){

  cout<<"Reading events from "<<data_file_path_<<" . . ."<<endl;
  TFile *f = open_root_file_();
  TTree *gamma_tree = get_gamma_tree_(*f);

  cout<<gamma_tree->GetEntries()<<" events in tree . . ."<<endl;

  // calculate cone data, returning cone data and number of cones calculated
  int num_cones = read_tree_into_vector_(*gamma_tree, conics, number_tries_per_random_point);
  cout<<"number cones processed "<<num_cones<<" . . ."<<endl;
  delete gamma_tree;
  f->Close();
  delete f;
}


float DetectorEffectsEventsLoader::get_weight(const Scatter &sc) const{
  static float c_E0=0.0, c_E1=0.0, c_E2=0.0, c_theta1=0.0, c_theta2 =0.0, c_dE1=0, c_dE2=0.0, c_dE3=0.0, c_intercept=0.0;
  static bool have_set_constants = false;

  if(! have_set_constants){
    c_E0 = strtof((*params_ptr_)["C_E0"].c_str(),0);
    c_E1 = strtof((*params_ptr_)["C_E1"].c_str(),0);
    c_E2 = strtof((*params_ptr_)["C_E2"].c_str(),0);
    c_theta1 = strtof((*params_ptr_)["C_THETA1"].c_str(),0);
    c_theta2 = strtof((*params_ptr_)["C_THETA2"].c_str(),0);
    c_dE1 = strtof((*params_ptr_)["C_dE1"].c_str(),0);
    c_dE2 = strtof((*params_ptr_)["C_dE2"].c_str(),0);
    c_dE3 = strtof((*params_ptr_)["C_dE3"].c_str(),0);
    c_intercept = strtof((*params_ptr_)["C_INTERCEPT"].c_str(),0);
  }
  float E1 = sc.getGammaEnergy() - sc.getScatter1EnergyDeposit();
  float E2 = sc.getGammaEnergy() - sc.getScatter1EnergyDeposit() - sc.getScatter2EnergyDeposit();

  float eta = c_intercept + c_E0*sc.getGammaEnergy() + c_E1*E1 + c_E2*E2
              + c_theta1*sc.getTheta1Degrees() + c_theta2*sc.getTheta2Degrees()
              + c_dE1*sc.getScatter1EnergyDeposit() + c_dE2*sc.getScatter2EnergyDeposit() + c_dE3*sc.getScatter3EnergyDeposit();

  float logistic_probability = 1.0/(1.0 + exp(-1.0 * eta));

  return logistic_probability;
}


size_t DetectorEffectsEventsLoader::read_tree_into_vector_(TTree &tree, vector< shared_ptr<ConicSection> > &conics, int number_tries_per_random_point) const {

  // Get number of TTree entries
  int rows = (int)tree.GetEntries();
  cout<<"read_tree_into_vector_: number entries in tree "<<tree.GetEntries()<<" . . ."<<endl;

  struct EventData{
    int event_num;
    PGVector3 cone_apex;
    PGVector3 cone_axis;
    PGVector3 true_origin;
    float x[3];
    float y[3];
    float z[3];
    float origin_x;
    float origin_y;
    float origin_z;
    float scatAng[3];
    float energy_deposited[3];
    float energy_incident[3];
    float energy_initial;
  } ed;

  tree.SetBranchAddress("event", &ed.event_num);

  tree.SetBranchAddress("pos_x", ed.x);
  tree.SetBranchAddress("pos_y", ed.y);
  tree.SetBranchAddress("pos_z", ed.z);
//   tree.SetBranchAddress("gammaEnergy", &ed.energy_true);
  tree.SetBranchAddress("gammaIncidentEnergy", &ed.energy_incident);
  tree.SetBranchAddress("gammaInitialEnergy", &ed.energy_initial);

  tree.SetBranchAddress("scattering_angle", ed.scatAng);
  tree.SetBranchAddress("energyDeposited", ed.energy_deposited);

  tree.SetBranchAddress("origin_x", &ed.origin_x);
  tree.SetBranchAddress("origin_y", &ed.origin_y);
  tree.SetBranchAddress("origin_z", &ed.origin_z);
  ed.true_origin = PGVector3(ed.origin_x,ed.origin_y,ed.origin_z);

//   tree.SetBranchAddress("origin_energy", ed.energy_true);

  PGVector3 p[3];
  PGVector3 random_point(phantom_volume_->x_max+1000.0,phantom_volume_->y_max+1000.0,phantom_volume_->z_max+1000.0);
  int reconstructable_count = 0; //counter for number of reconstructable
  int number_ellipses = 0; //number cones intersecting as ellipse
  int number_parabolas = 0; //number cones intersecting as parabolas
  int number_cones_skipped = 0; //count of number of unreconstructable conics
  int nan_skipped = 0;
  int randoms_not_in_cone_skipped = 0;
  int low_energy_skipped = 0;
  int high_energy_skipped = 0;
  int DCA_max_skipped = 0;
  int DCA_min_skipped = 0;
  int dca_cut_skipped = 0;

  double min_energy = strtod((*params_ptr_)["MIN_GAMMA_ENERGY"].c_str(),0);
  double max_energy = strtod((*params_ptr_)["MAX_GAMMA_ENERGY"].c_str(),0);
  double DCA_max_cutoff = strtod((*params_ptr_)["DCA_MAX"].c_str(),0);
  double DCA_min_cutoff = strtod((*params_ptr_)["DCA_MIN"].c_str(),0);


  PGVector3 dca_cut_point;
  dca_cut_point.x = params_ptr_->get_float("DCA_CUT_X");
  dca_cut_point.y = params_ptr_->get_float("DCA_CUT_Y");
  dca_cut_point.z = params_ptr_->get_float("DCA_CUT_Z");
  float dca_cut = params_ptr_->get_float("DCA_CUT");


  //number of attempts needed to find random point in phantom
  // for all of the cones
  int number_random_points_tried = 0;
  int num_cones_requested  = atoi((*params_ptr_)["MAX_NUM_CONES"].c_str());
  if( 1 > num_cones_requested){
    string error_message("ERROR: invalid limit on number of cones(");
    error_message += (*params_ptr_)["MAX_NUM_CONES"].c_str();
    error_message += ")\nCheck the value for parameter MAX_NUM_CONES in the parameters file.\n";
    error_message += ".\nABORTING\n\n";
    throw runtime_error(error_message);
  };

  int offset  = atoi((*params_ptr_)["NUM_CONES_OFFSET"].c_str());
  cout<<"GAMMA OFFSET: "<<offset<<", "<<num_cones_requested<<" cones requested, "<<tree.GetEntries()<<" cones total . . ."<<endl;

  for(int i=offset; i< rows; ++i){
    tree.GetEntry(i);
    p[0] = PGVector3(ed.x[0],ed.y[0], ed.z[0]);
    p[1] = PGVector3(ed.x[1], ed.y[1], ed.z[1]);
    p[2] = PGVector3(ed.x[2], ed.y[2], ed.z[2]);
    ed.true_origin = PGVector3(ed.origin_x,ed.origin_y,ed.origin_z);

    shared_ptr<TripleScatter> ts(new TripleScatter(ed.energy_deposited[0], ed.energy_deposited[1], ed.energy_deposited[2], p[0], p[1], p[2]));

    while(1) {
      AddDetectorEffects(*ts);
      ///if energy is NAN than add effects again.
      if( ts->getGammaEnergy() == ts->getGammaEnergy() ) break;
    }

    if(ts->getGammaEnergy() <min_energy){
      ++number_cones_skipped;
      ++low_energy_skipped;
      continue;
    }

    if(ts->getGammaEnergy() >max_energy){
      ++number_cones_skipped;
      ++high_energy_skipped;
//       printf("high %d of %d, %d . . .\n", i, num_cones_requested,reconstructable_count);
      continue;
    }

    double scattering_angle = ts->getConeOpeningAngle();

    //scattering_angle in not a number, skip it
    if( scattering_angle != scattering_angle){
      ++number_cones_skipped;
      ++nan_skipped;
      continue;
    }
    shared_ptr<ComptonScatter> tmpComptonScatter;
    try{
      tmpComptonScatter = shared_ptr<ComptonScatter> (new ComptonScatter(
                                                        ts->getConeApex(),
                                                        ts->getConeAxis(),
                                                        ed.true_origin,
                                                        ts->getConeOpeningAngle(),
                                                        ts->getScatter1EnergyDeposit(),
                                                        ts->getScatter2EnergyDeposit(),
                                                        ts->getGammaEnergy()
                                                      ));
    }catch( runtime_error ){
      ++number_cones_skipped;
//       printf("runtime_error %d of %d, %d . . .\n", i, num_cones_requested,reconstructable_count);
      continue;
    }
    double alpha = tmpComptonScatter->getAlpha();
    double DCA = 0.0;
    shared_ptr<ConicSection> ptrConicSection;
    //phantom_volume_ptr_->print();

    if( alpha + ts->getConeOpeningAngle() > M_PI/2.0){//then conic section is parabola
      ++number_parabolas;
      ptrConicSection = shared_ptr<ConicSection>( new ReconstructionParabola( *tmpComptonScatter, phantom_volume_,i));
    }else{//then conic section is an ellipse
      ptrConicSection = shared_ptr<ConicSection>( new ReconstructionEllipse( *tmpComptonScatter, phantom_volume_,i));
      ++number_ellipses;
    }

    if(ptrConicSection->getDistanceToPoint(dca_cut_point) > dca_cut){
      ++dca_cut_skipped;
      continue;
    }

    PGVector3 random_point;
    int number_tries = ptrConicSection->getRandomPointInPhantom(random_point, number_tries_per_random_point);
    if( -1 == number_tries ){
      ++number_cones_skipped;
      ++randoms_not_in_cone_skipped;
//       printf("number_tries %d of %d, %d . . .\n", i, num_cones_requested,reconstructable_count);
      continue;
    }
    ptrConicSection->setMCTruth(ed.scatAng, p, ed.energy_deposited, ed.energy_incident, ed.energy_initial, &ed.origin_x);
//     ts->print();
    ptrConicSection->setScatter(ts);

//    double weight = get_weight(*ts);
////     printf("WEIGHT, %.5f, %.5f\n", weight, DCA);
//    double weight_power = strtod((*params_ptr_)["WEIGHT_POWER"].c_str(), 0);
//
////    ptrConicSection->setWeight(pow(weight,weight_power));
    ptrConicSection->setWeight(1.0);
    conics.push_back(ptrConicSection);

    ++reconstructable_count;

    if(num_cones_requested <= reconstructable_count){
      printf("%d cones, last cone loaded is %d . . .\n", reconstructable_count, i);
      break;
    }

    /// Monitor progress
    int increment = (10000 < num_cones_requested)? num_cones_requested/ 10: 1000;
    if (i % increment == 0) {
      printf("    Row Number: %d, Event Number: %i, Reconstructable: %d\n", i, ed.event_num, reconstructable_count);
    }
  }
  printf("--- Number of parabolas in TTree: %d ---\n", number_parabolas);
  printf("--- Number of ellipses in TTree: %d ---\n", number_ellipses);
  printf("--- Number of skipped due to nan angle: %d ---\n", nan_skipped);
  printf("--- Number of skipped due to energy below %f: %d ---\n", min_energy, low_energy_skipped);
  printf("--- Number of skipped due to energy above %f: %d ---\n", max_energy, high_energy_skipped);
  printf("--- Number of skipped due to DCA CUT: %d ---\n", dca_cut_skipped);
  printf("--- Number of skipped due to DCA above %f: %d ---\n", DCA_max_cutoff, DCA_max_skipped);
  printf("--- Number of skipped due to DCA below %f: %d ---\n", DCA_min_cutoff, DCA_min_skipped);
//   printf("--- Number of skipped due to no cone phantom intercept: %d ---\n", not_in_phantom_skipped);

  printf("--- Number of skipped due to no randoms in the phantom: %d ---\n", randoms_not_in_cone_skipped);


  printf("--- Number of useable triple scatters in TTree: %d ---\n", reconstructable_count);

  if(0 < reconstructable_count){
    printf("--- Avg. randoms points tested per conic section: %.2f ---\n",
     static_cast<float>(number_random_points_tried)/static_cast<float>(reconstructable_count));
  }

  return reconstructable_count;
}

