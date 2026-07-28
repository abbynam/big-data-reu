/*! \mainpage ConicRecon - A Compton Gamma Camera Reconstruction Algorithm
 *
 * \section intro_sec Overview
 *
 * ConicRecon reconstructs two-stage and three-stage Compton cameras. The 
 * intended use in vivo imaging during proton therapy. The algorithm is 
 * based on the paper 
 * <i>Fast image reconstruction for Compton camera using stochastic origin ensemble approach</i>
 * by Andreyev, Sitek, and Celler, Med. Phys. 2011. However, we adapted the
 * algorithm where appropriate to improve performance, both processing speed and
 * image quality.
 * 
 * This code relies heavily on the <a href="http://root.cern.ch">ROOT</a> 
 * analysis framework developed originally for processing high energy
 * physics data. ROOT has some built in memory management so code using 
 * ROOT objects will have naked news and deletes.
 * 
 * One aspect of the style is to prefer descriptive name to comments answering 
 * "what?". Comments should be provided when an explanation of "why?" is needed.
 * 
 * This code is intended for research purposes only. 
 * 
 * \section building Building the Code
 * This code based has been developed and tested on Linux Mint with kernel 
 * <i>Linux 3.2.0-23-generic x86_64</i>. Because this application relies on 
 * compiled ROOT libraries, it may not link properly on other Linux systems. 
 * 
 * Expand sample code tar file to the directory <kbd>[some path]/sample</kbd>
 * and the sample packages to <nobr><kbd>[some path]/sample_packages</kbd></nobr>.
 * 
 * Then, run <kbd>make</kbd> from a command prompt with the current 
 * working directory set to  <kbd>[some path]/sample</kbd>. To build the 
 * tests and the documentation enter <kbd>make test</kbd> 
 * and <kbd>doxygen Doxyfile.cfg</kbd>, respectively. Doxygen must be
 * installed on the system for the documentation to build correctly.  
 * 
 * 
 * \section Instructions Running the Code
 * To add the ROOT libraries in <kbd>sample_packages</kbd> to 
 * the library path, first run <kbd>source setup.sh</kbd>.
 * To run the application, type <kbd>./conicrecon soe sample.cfg</kbd>. 
 * To run the tests, type <kbd>./unit_tests</kbd>. Both commands require
 * that the current working directory be set to <kbd>[some path]/sample</kbd>.
 * 
 * @author Dennis Mackin
 */

//standard includes
#include <iostream>
#include <string>
#include <memory>
#include <cmath>


//local includes
#include "ReconstructorBuilder.h"
#include "ReconstructorTemplate.h"
//#include "StochasticReconstructor.h"
//#include "ConicEnsemble.h"
//;#include "DensityReconstructor.h"
//#include "HistBasedReconstructor.h"
//#include "DCAKernelReconstructor.h"
//#include "NearestNeighborReconstructor.h"
#include "RunTimeParameters.h"

using namespace std;
using namespace prompt_gamma_reconstruction;
//
//int runHistBased(const string &parameters_file){
//    std::shared_ptr<HistBasedReconstructor> myRecoPtr(new HistBasedReconstructor(parameters_file));
//    myRecoPtr->StartExperiment();
//    return 0;
//}
//
//int runStochastic(const string &parameters_file){
//    std::shared_ptr<StochasticReconstructor> myRecoPtr(new StochasticReconstructor(parameters_file));
//    myRecoPtr->start_calc();
//    return 0;
//}
//
//int runDensityReconstructor(const string &parameters_file){
//
//    std::shared_ptr<DensityReconstructor> myRecoPtr(new DensityReconstructor(parameters_file));
//    cout<<"Completed creating DensityReconstructor object . . ."<<endl;
//    myRecoPtr->StartExperiment();
//
//    return 0;
//}
//
//int runDCAKernelReconstructor(const string &parameters_file){
//
//  std::shared_ptr<DCAKernelReconstructor> myRecoPtr(new DCAKernelReconstructor(parameters_file));
//  cout<<"Completed creating DCAKernelReconstructor object . . ."<<endl;
//  myRecoPtr->StartExperiment();
//
//  return 0;
//}
//
//int runNearestNeighborlReconstructor(const string &parameters_file){
//
//    std::shared_ptr<NearestNeighborReconstructor> myRecoPtr(new NearestNeighborReconstructor(parameters_file));
//    cout<<"Completed creating DCAKernelReconstructor object . . ."<<endl;
//    myRecoPtr->StartExperiment();
//
//    return 0;
//}
//

size_t runReconstructor(const string &parameters_file){
    auto builder = ReconstructorBuilder();
    auto reconstructor = builder.build(parameters_file);
    reconstructor->run();

    return 0;
}


size_t usage(char *name){
    string USAGE("usage: ");
    USAGE.append(name);
    //USAGE.append(" [soe | kernel | hist | dca | NN | dyn] [parameters file]\n");
    USAGE.append(" [parameters file]\n");
    cerr<<USAGE;
    return -1;
}


int main(int argc, char *argv[]){

    if( argc != 2){
      return usage(argv[0]);
    }

    //std::shared_ptr<TStopwatch> myTimer(TStopwatch());
//    auto myTimer = std::make_shared<TStopwatch>(TStopwatch());
//    myTimer->Start();

//    string run_type(argv[1]);
    string parameters_file(argv[1]);

//    if("dyn" == run_type){
//        runReconstructor(parameters_file);
//    }else{
//        return usage(argv[0]);
//    }
    runReconstructor(parameters_file);

//    myTimer->Stop();
  //  Double_t real_time = myTimer->RealTime();

    //cout<<"Run Time "<<run_type<<": "<<real_time<<"(s) . . ."<<endl;
    
    return 0;
}
