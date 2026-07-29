// //
// DensityReconstructor
//
// v1: created reconstruction GUI using stochastic methods
// v2: modified field-of-view specifications from number of voxels to image dimensions
// v3: added timer (not necessary, found right refresh commands) to periodically update GUI (textbox [finished], graphs [not finished])
//
//  ---- added by Dennis Mackin 2010-09-20
// v4: designed to compile to executable
// v5: replaced GUI interface with a config file

// Automatically compiles code / correct format: root reconGUIv3.C+
//  - must be in same directory (would like to change this)
#include <memory>
#include "TStyle.h"

#include "DensityReconstructor.h"
#include "RandomSingleton.h"
#include "ConicEnsemble.h"



using namespace prompt_gamma_reconstruction;
using namespace pg_tools;

// Create the main frame
//DensityReconstructor::DensityReconstructor(const TGWindow *p, UInt_t w, UInt_t h) {
DensityReconstructor::DensityReconstructor(const string &parameters_file_path):
    parameters_file_path_(parameters_file_path),
    run_time_parameters_(parameters_file_path),
    num_cones_(0){

  // Set starting timestamp
  start_time_ = shared_ptr<TTimeStamp>(new TTimeStamp());

  output_folder_path_ =  setup_output_folder_(run_time_parameters_["OUTPUT_FOLDER_PATH"]);
  output_root_file_path_ = output_folder_path_ ;
  output_root_file_path_ +=  run_time_parameters_["FOLDER_DELIMITER"];
  output_root_file_path_ += run_time_parameters_["RUN_NAME"] + ".root"; //set in _openOutputFile

  smoothing_parameter_ = std::strtod(run_time_parameters_["SMOOTHING_PARAMETER"].c_str(),0);
  number_kernels_ = StringUtils::strtoi(run_time_parameters_["NUMBER_KERNELS"]);
  events_file_path_ = run_time_parameters_["EVENT_FILE_PATH"]; //path to file with gamma dections
  width_ = std::strtod(run_time_parameters_["KERNEL_WIDTH"].c_str(),0);

  ensemble_ = shared_ptr<ConicEnsemble>(new ConicEnsemble(events_file_path_,smoothing_parameter_,width_));
  ensemble_->setSmoothingParameterForEnsemble(smoothing_parameter_, number_kernels_);

  printf("processing %d cones . . .\n", static_cast<int>(ensemble_->size()));
  cout<<"completed construction of DensityReconstructor . . ."<<endl;
  printf("--- Run Time (start up): %f sec ---\n", getRunningTime());
}


// Draws function graphics in randomly choosen interval
void DensityReconstructor::StartExperiment() {
  cout<<"starting the experiment . . ."<<endl;
  produce_output_(run_time_parameters_["RUN_NAME"]);
}

void DensityReconstructor::open_output_file_(const string &name){
  cout<<"opening output file "<< output_root_file_path_ <<" . . ."<<endl;
  output_root_file_ = new TFile(output_root_file_path_.c_str(),"RECREATE");
  if(! output_root_file_->IsOpen()){
    string error_message("DensityReconstructor::open_output_file_: failed to open ");
    error_message += name;
    error_message += ".\nABORTING\n\n";
    cout<<error_message<<endl;
    throw runtime_error(error_message);
  };
}

void DensityReconstructor::produce_1D_plot_(TH1D &hist){

  cout<<"making 1D plot . . ."<<endl;
  //produce 1D histogram

  int bins = run_time_parameters_.get_int("BINS_1D");
  double min = run_time_parameters_.get_double("MIN_1D");
  double max = strtod(run_time_parameters_["MAX_1D"].c_str(),0);
  int number_kernels_1D = StringUtils::strtoi(run_time_parameters_["NUMBER_KERNELS_1D"].c_str());
  double smoothing_parameter_1D = strtod(run_time_parameters_["SMOOTHING_PARAMETER_1D"].c_str(),0);
  ensemble_->setSmoothingParameterForEnsemble(smoothing_parameter_1D, number_kernels_1D);

  TH1D myHist("myHist","Prompt Gamma Density", bins, min, max);

  int dimension = StringUtils::strtoi(run_time_parameters_["DIMENSION_1D"]);
  dimension %= 3;
  PGVector3 point(0.0,0.0,0.0);
  double xyz[] = {0.0,0.0,0.0};
  for(int i = 1; i <= bins; ++i){
    xyz[dimension] = myHist.GetBinCenter(i);
    point.x = xyz[0];
    point.y = xyz[1];
    point.z = xyz[2];
    double density = ensemble_->getDensity(point, number_kernels_1D);
    myHist.SetBinContent(i,density);
    if(i%100 == 0 ) std::cout<<"1D hist bin "<<i<<" density: "<< density<<" . . ."<<std::endl;
  }
  myHist.Scale(1.0/myHist.GetMaximum());

  hist = myHist;
}

void DensityReconstructor::produce_2D_plot_(TH2D &hist){

  cout<<"making 2D plot . . ."<<endl;
  int bins_x = StringUtils::strtoi(run_time_parameters_["BINS_X"]);
  double min_x = strtod(run_time_parameters_["MIN_X"].c_str(),0);
  double max_x = strtod(run_time_parameters_["MAX_X"].c_str(),0);
  double bin_size_x = (max_x - min_x)/static_cast<double>(bins_x);

  int bins_z = StringUtils::strtoi(run_time_parameters_["BINS_Z"]);
  double min_z = strtod(run_time_parameters_["MIN_Z"].c_str(),0);
  double max_z = strtod(run_time_parameters_["MAX_Z"].c_str(),0);
  double bin_size_z = (max_z - min_z)/static_cast<double>(bins_z);

  int number_kernels_2D = StringUtils::strtoi(run_time_parameters_["NUMBER_KERNELS_2D"].c_str());
  ensemble_->setSmoothingParameterForEnsemble(smoothing_parameter_, number_kernels_2D);

  TH2D myHist2D("myHist 2D","Prompt Gamma Density XZ", bins_z, min_z, max_z, bins_x, min_x, max_x);
  PGVector3 point(0.0,0.0,0.0);
  double sum_of_bins = 0.0;
  double density = 0.0;
  for(int i = 1; i <= bins_x; ++i){
    point.y = (static_cast<double>(i) - 0.5)*bin_size_x + min_x;
    for(int j = 1; j <= bins_z; ++j){
      point.z = (static_cast<double>(j) - 0.5)*bin_size_z + min_z;
      density = ensemble_->getDensity(point, number_kernels_2D);
      myHist2D.SetBinContent(j, i, density);
      sum_of_bins += density;

    }
    if(i%10 == 0) std::cout<<"bin("<<i<<",0); "<<point.print()<<" density: "<< density<<" . . ."<<std::endl;
  }
  //myHist2D.Scale(1.0/sum_of_bins);
  myHist2D.Scale(1.0/myHist2D.GetMaximum());

  gStyle->SetPalette(1);

  hist = myHist2D;
}


TCanvas *DensityReconstructor::setup_canvas_(){
  cout<<"setting up canvas . . ."<<endl;
  int width = StringUtils::strtoi(run_time_parameters_["CANVAS_WIDTH"]);
  int height = StringUtils::strtoi(run_time_parameters_["CANVAS_HEIGHT"]);
  TCanvas *myCanvas = new TCanvas("myTest","Density Canvas", width, height);
  myCanvas->SetRightMargin(myCanvas->GetRightMargin()*1.5);
  myCanvas->SetLeftMargin(myCanvas->GetLeftMargin()*1.1);
  myCanvas->SetBottomMargin(myCanvas->GetBottomMargin()*1.1);

    cout<<"completed setting up cavas . . ."<<endl;
  return myCanvas;
}


void DensityReconstructor::produce_output_(const string &output_name){
  //////// Store start data and create image for GUI and to store to file

  open_output_file_(output_name);
  assert(output_root_file_->IsOpen());
  output_root_file_->cd();

  TH1D hist1D;
  produce_1D_plot_(hist1D);

  auto_ptr<TCanvas> myCanvas(setup_canvas_());

  myCanvas->Draw();
  hist1D.Draw("L");
  string hist1D_file = output_folder_path_;
  hist1D_file += run_time_parameters_["FOLDER_DELIMITER"];
  hist1D_file += run_time_parameters_["RUN_NAME"];
  hist1D_file += "_" + run_time_parameters_["NAME_HIST_1D"];

  string hist1D_file_C = hist1D_file + ".C";
  hist1D_file += ".gif";
  myCanvas->SaveAs(hist1D_file.c_str());
  myCanvas->SaveAs(hist1D_file_C.c_str());

  TH2D hist2D;
  produce_2D_plot_(hist2D);

  string hist2D_file(output_folder_path_);
  hist2D_file += run_time_parameters_["FOLDER_DELIMITER"] + run_time_parameters_["RUN_NAME"];
  hist2D_file += "_" + run_time_parameters_["NAME_HIST_2D"];
  string hist2D_file_C(hist2D_file);
  hist2D_file_C += ".C";
  hist2D_file += ".gif";

  hist2D.SetStats(0);
  hist2D.Draw("CONTZ");
  gStyle->SetPalette(1);//set the color palette to the visible light spectrum
  myCanvas->SaveAs(hist2D_file.c_str());
  myCanvas->SaveAs(hist2D_file_C.c_str());

  printf("--- Storing data to root file: %s ---\n", output_root_file_path_.c_str());
  hist1D.Write();
  hist2D.Write();
  output_root_file_->Close();
}

DensityReconstructor::~DensityReconstructor() {
  if(output_root_file_){
    output_root_file_->Write();
    output_root_file_->Close();
    delete output_root_file_;
  }
  cout<<"completed destruction of DensityReconstructor . . ."<<endl;
}


//////////////////////////////////////////////////////////////////////////////////////
string DensityReconstructor::setup_output_folder_(const string &outputFolderPath){
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


  string output_folder = run_time_parameters_["OUTPUT_FOLDER_PATH"];
  output_folder += "_";
  output_folder += run_time_parameters_["SMOOTHING_PARAMETER"];
  output_folder += "_";
  output_folder += run_time_parameters_["NUMBER_KERNELS"];

  char command[501];
  //check to see if folder exists
  sprintf(command, "%s %s",
	  run_time_parameters_["MAKE_DIR_COMMAND"].c_str(),
	  output_folder.c_str());

  cout<<"running command {"<<command<<"} . . ."<<endl;

  if(system(command) != 0){
    printf("failed to make %s. \nDoes it exist? If not, there is a problem . . .\n",
      run_time_parameters_["OUTPUT_FOLDER_PATH"].c_str());
  }else{
    cout<<"Created directory "<< run_time_parameters_["OUTPUT_FOLDER_PATH"]<<". . . \n";
  }

  cout<<"Copying the parameters file "<< parameters_file_path_ <<" . . ."<<endl;
  //copy the parametersFile and the data files to the new folder
  FileUtils::fileCopy(parameters_file_path_, output_folder, run_time_parameters_["COPY_COMMAND"]);

  cout<<"completed setting up output folder . . ."<<endl;
  return output_folder;
}//end of SetUpOutputFolder
//////////////////////////////////////////////////////////////////////////////////////

