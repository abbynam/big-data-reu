#include "TFile.h"
#include "ConicEnsemble.h"
#include "RandomSingleton.h"


using namespace prompt_gamma_reconstruction;


/// CONSTRUCTOR and DESTRUCTOR
ConicEnsemble::ConicEnsemble(const string &root_input_file_path, double smoothing, double width): h_(smoothing), sigma_(width){

  readInConicSections_(root_input_file_path);
  setCoefficient();
}
 
ConicEnsemble::~ConicEnsemble() { /*  */  }


// opens root file (contains tripleGammmas) and calculates cone data (apexPos, scatVect, scatAng)
void ConicEnsemble::readInConicSections_(const string &root_input_file_path){
  cout<<"Reading gammas from "<<root_input_file_path<<" . . ."<<endl;
  TFile *f = TFile::Open(root_input_file_path.c_str(), "READ");
  if( !(f->IsOpen()) ){
    printf("ERROR: root file \n%s\ndid not open properly. ABORTING . . .\n",root_input_file_path.c_str());
    exit(-111);
  }
  TTree *scatterTree;
  f->GetObject("tripleGammas", scatterTree);

  ///////// Creating cone data (position, vector and angle) from triple scatter data / return number of cones
  // calculate cone data, returning cone data and number of cones calculated
  number_cones_ = loadConicSections_(scatterTree);
  printf("%d cones calculated . . .\n", number_cones_);

  delete scatterTree;
  f->Close();
  delete f;
}

// function accepts TTree contains triple scatter data and returns cone data (cone apex position, vector between steps 1 & 0, scatter angle)
int ConicEnsemble::loadConicSections_(TTree *tree) {

  // Variable declaration
  Double_t threePos[3][3], direction[3];
  Int_t event, step;
  Double_t pos_x, pos_y, pos_z, scatAng;
  char detector[300], process[300];
  string detectorName, processName;

  // declare tree and branches
  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("step", &step);
  tree->SetBranchAddress("detector", &detector);
//  tree->SetBranchAddress("detector_num", &detectorNum);
  tree->SetBranchAddress("process", &process);
  tree->SetBranchAddress("pos_x", &pos_x);
  tree->SetBranchAddress("pos_y", &pos_y);
  tree->SetBranchAddress("pos_z", &pos_z);
  tree->SetBranchAddress("scatAng", &scatAng);


  // Get number of TTree entries
  int rows = (int)tree->GetEntries();
  printf("--- Number of Rows in TTree: %d ---\n", rows);

  // Initialize triple scatter count
  int count = 0;

  PGVector3 axis, apex, true_origin;
  // Cycle through triple scatter data
  for (int i = 0; i < rows; i++) {
    // Load data
    tree->GetEntry(i);
    // Convert char arrays into strings
    detectorName = detector;
    processName = process;

    // Check that first step is Comptom event with scatter angle between 0 and 85 degrees
    if (  scatAng < 85.0 && scatAng > 0.0
	  && (step == 0 || detectorName.find("One") != string::npos)
	  && processName.find("Compton") != string::npos ) {

      // Store position
      for (int j = 0; j < 3; j++) {
        tree->GetEntry(i + j);
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
      tree->GetEntry(i);

      // Store data, reverse axis direction to down from up

      axis.x = -direction[0];
      axis.y = -direction[1];
      axis.z = -direction[2];

      apex.x = pos_x;
      apex.y = pos_y; //make y a positive value (I don't know why it is negative)
      apex.z = pos_z;

      ComptonScatter tmpComptonScatter = ComptonScatter( apex, axis, true_origin,
                                                         scatAng, 0.0, 0.0, 0.0 );

      //KernelConic *ptrConic = new KernelConic( tmpComptonScatter, &phantom_volume_ptr_, 0.31, 10.0);
      KernelConic conic = KernelConic( tmpComptonScatter, phantom_volume_);

      kernel_conic_sections_.push_back(conic);
      count++;
    }
    // Monitor progress
    int increment = rows / 10;
    if (i % increment == 0) {
      printf("    Row Number: %d, Event Number: %i\n", i, event);
    }
  }
  printf("--- Number of useable triple scatters in TTree: %d ---\n", count);
  std::random_shuffle(kernel_conic_sections_.begin(), kernel_conic_sections_.end());
  printf("--- Finished shuffling ---\n");
  setSmoothingParameterForEnsemble(count);
  return count;
}

//////////////////////////////////////////////////////////////////////////////////////
void ConicEnsemble::setSmoothingParameterForEnsemble(double h, int number_cones){
//////////////////////////////////////////////////////////////////////////////////////
/// Sets the smoothing parameter (h) for kernels. If h is 0 or < 0, then
/// a plugin estimator for h is used based on the number of dimensions and
/// the number of conic sections used in the estimatation.
///
/// @param h The smoothing parameter.
/// @param number_cones The number of cones used to estimate the density.
/// @returns void
//////////////////////////////////////////////////////////////////////////////////////

    if(0 >= h) {
        h = pow(4.0/3.0, 0.2)*pow(number_cones,-0.20);
    }

    cout<<"Setting smoothing paramaters to "<< h<<" . . ."<<endl;
    h_ = h;
    setExponentialDenominator_();
    return;
}//end of setSmoothingParameterForEnsemble
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////
void ConicEnsemble::setWidthForEnsemble(double width){
//////////////////////////////////////////////////////////////////////////////////////
///
/// Sets the width of the Gaussian kernels.
///
/// @param point The point in 3 space where density is to be tested.
/// @returns normalized density
//////////////////////////////////////////////////////////////////////////////////////

    cout<<"Setting width for kernels to "<< width<<" . . ."<<endl;
    sigma_ = width;

    setExponentialDenominator_();
    return;
}//end of setSmoothingParameterForEnsemble
//////////////////////////////////////////////////////////////////////////////////////


double ConicEnsemble::getDensityForPoint_(const double &distance) const{

    double ratio = distance*exponential_denominator_;
    if(8.0 < ratio) return 0.0;
    double density = coefficient_ * exp( -0.5*ratio*ratio);

    return density;
}

//////////////////////////////////////////////////////////////////////////////////////
double ConicEnsemble::getDensity(const PGVector3 &point, const int number_conics = 1E9) const{
//////////////////////////////////////////////////////////////////////////////////////
///
/// Sums the density contributioins from alll conic sections in the ensemble and
/// divides by the number of conics to normalize the result.
///
/// @param point The point in 3 space where density is to be tested.
/// @returns normailzed density
//////////////////////////////////////////////////////////////////////////////////////

    vector<KernelConic>::const_iterator iter = kernel_conic_sections_.begin();
    vector<KernelConic>::const_iterator end_ptr = kernel_conic_sections_.end();

    double density = 0.0;
    int counter = 0;

    for(/* */; iter != end_ptr; ++iter){
        double distance = iter->getDistanceToPoint(point);
        density += getDensityForPoint_(distance);
        ++counter;
        if( number_conics <= counter) break;
    }
    density /= static_cast<double>(counter);

    return density;
}//end of getDensity
//////////////////////////////////////////////////////////////////////////////////////

