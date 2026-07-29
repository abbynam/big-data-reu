#include "PGSimulatedEventsLoader.h"


// Standard C++ Includes
#include <cmath>
#include <valarray>
#include <vector>
#include <stdexcept>
#include <iostream>
//

//PG Includes
#include "PhantomVolume.h"
#include "ReconstructionParabola.h"
#include "ReconstructionEllipse.h"
#include "ConicSection.h"
#include "EventsLoader.h"


using namespace std;
using namespace prompt_gamma_reconstruction;

PGSimulatedEventsLoader::PGSimulatedEventsLoader(const string &data_file_path, const pg_tools::RunTimeParameters *params, shared_ptr<const PhantomVolume>phantom):
   EventsLoader(data_file_path, params, phantom){

    gamma_tree_name_ = (*params)["GAMMA_TREE_NAME"];
}


void PGSimulatedEventsLoader::LoadEvents(vector<shared_ptr<ConicSection> > &conics, int number_tries_per_random_point){

  cout<<"Reading events from "<<data_file_path_<<" . . ."<<endl;
  TFile *f = open_root_file_();
  TTree *gamma_tree = get_gamma_tree_(*f);

  cout<<"number entries in tree "<<gamma_tree->GetEntries()<<" . . ."<<endl;

  // calculate cone data, returning cone data and number of cones calculated
  int num_cones = read_tree_into_vector_(*gamma_tree, conics, number_tries_per_random_point);
  cout<<"number entries in tree after loading "<<gamma_tree->GetEntries()<<" . . ."<<endl;
  printf("%d cones calculated . . .\n", num_cones);

  delete gamma_tree;
  f->Close();
  delete f;
}

// function accepts TTree contains triple scatter data and returns cone data (cone apex position, vector between steps 1 & 0, scatter angle)
size_t PGSimulatedEventsLoader::read_tree_into_vector_(TTree &tree, vector<shared_ptr<ConicSection> > &conics, int number_tries_per_random_point) {

  // Variable declaration
  Double_t threePos[3][3], direction[3];
  Int_t event, step;
  Double_t pos_x, pos_y, pos_z, scatAng;
  char detector[300], process[300];
  string detectorName, processName;

  // declare tree and branches
  tree.SetBranchAddress("event", &event);
  tree.SetBranchAddress("step", &step);
  tree.SetBranchAddress("detector", &detector);
//  tree.SetBranchAddress("detector_num", &detectorNum);
  tree.SetBranchAddress("process", &process);
  tree.SetBranchAddress("pos_x", &pos_x);
  tree.SetBranchAddress("pos_y", &pos_y);
  tree.SetBranchAddress("pos_z", &pos_z);
  tree.SetBranchAddress("scatAng", &scatAng);


  // Get number of TTree entries
  int rows = (int)tree.GetEntries();
  printf("--- Number of Rows in TTree: %d ---\n", rows);

  // Initialize triple scatter count
  int count = 0;
  int num_parabolas = 0;
  int num_ellipses = 0;
  int num_skipped = 0;
  int num_random_points_tried = 0;

  //initialize the random_point to somewhere outside the phantom.
  PGVector3 random_point(phantom_volume_->x_max+1000.0,phantom_volume_->y_max+1000.0,phantom_volume_->z_max+1000.0);
  PGVector3 axis, apex, true_origin;

  // Cycle through triple scatter data
  for (int i = 0; i < rows; i++) {
    // Load data
    tree.GetEntry(i);
    // Convert char arrays into strings
    detectorName = detector;
    processName = process;
/*
    // Check that first step is Comptom event with scatter angle between 0 and 90 degrees
    if ((step == 0 || detectorName.Contains("One"))
      && processName.Contains("Compton") && scatAng > 0.0 && scatAng < 90.0) {
*/
    // Check that first step is Comptom event with scatter angle between 0 and 85 degrees
    if (  scatAng < 85.0 && scatAng > 0.0
    && (step == 0 || detectorName.find("One") != string::npos)
    && (processName.find("Compton") != string::npos || processName.find("compt") != string::npos) ) {

      // Store position
      for (int j = 0; j < 3; j++) {
        tree.GetEntry(i + j);
        threePos[j][0] = pos_x;
        threePos[j][1] = pos_y;
        threePos[j][2] = pos_z;
      }

      // Calculate vector between step 0 and step 1 (subtract step 0 from step 1 at each dimension / flips direction)
      //DSM don't flip the direction
      for (int k = 0; k < 3; k++) {
        //direction[k] = threePos[0][k] - threePos[1][k];
        direction[k] =  threePos[1][k] - threePos[0][k];
      }
      //reload the first event in the scatter
      tree.GetEntry(i);

      // Store data, reverse axis direction to down from up

      axis.x = -direction[0];
      axis.y = -direction[1];
      axis.z = -direction[2];

      apex.x = pos_x;
      apex.y = pos_y; //make y a positive value (I don't know why it is negative)
      apex.z = pos_z;

      //ComptonScatter tmpComptonScatter( ComptonScatter(axis, apex, scatAng, 0.0, 0.0, 0.0) );

      ///@TODO add the scattering energies to the comptonScatters
      ComptonScatter tmpComptonScatter = ComptonScatter(apex, axis, true_origin,
                                                        scatAng*M_PI/180.0, 0.0, 0.0, 0.0 );

      double alpha = tmpComptonScatter.getAlpha();
      shared_ptr<ConicSection> ptrConicSection;

      int number_tries = 0; //number of randoms thrown before getting a point in phantom
      if( alpha + scatAng*M_PI/180.0 > M_PI/2.0){//then conic section is parabola
        //ptrConicSection = new ReconstructionParabola( tmpComptonScatter, phantom_volume_ptr_);
        ptrConicSection = shared_ptr<ConicSection>( new ReconstructionParabola( tmpComptonScatter, phantom_volume_,i));
        ++num_parabolas;
      }else{//then conic section is an ellipse
        //ptrConicSection = new ReconstructionEllipse( tmpComptonScatter, phantom_volume_ptr_);
        ptrConicSection = shared_ptr<ConicSection>( new ReconstructionEllipse( tmpComptonScatter, phantom_volume_,i));
        ++num_ellipses;
      }
      number_tries = ptrConicSection->getRandomPointInPhantom(random_point, number_tries_per_random_point);
      if( -1 == number_tries || ! ptrConicSection->doesConicIntersectPhantom()){
        ++num_skipped;
//         delete ptrConicSection;
        continue;
      }

      num_random_points_tried += number_tries;
      conics.push_back(ptrConicSection);
      ptrConicSection->setLikelyOrigin(random_point);

      ++count;
    }

    // Monitor progress
    int increment = rows / 10;
    if (i % increment == 0) {
      printf("    Row Number: %d, Event Number: %i\n", i, event);
    }
  }
  printf("--- Number of parabolas in TTree: %d ---\n", num_parabolas);
  printf("--- Number of ellipses in TTree: %d ---\n", num_ellipses);
  printf("--- Number of useable triple scatters in TTree: %d ---\n", count);

  if(0 < count){
    printf("--- Avg. randoms points tested per conic section: %.2f ---\n",
     static_cast<float>(num_random_points_tried)/static_cast<float>(count));
  }
  return count;

}


