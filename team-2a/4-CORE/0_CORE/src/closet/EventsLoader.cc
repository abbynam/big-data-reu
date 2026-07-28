#include "PGSimulatedEventsLoader.h"

// Standard C++ Includes
#include <cmath>
#include <valarray>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <sstream>

// ROOT Includes
#include "TFile.h"
#include "TTree.h"

//PG Includes
#include "PhantomVolume.h"
#include "ReconstructionParabola.h"
#include "ReconstructionEllipse.h"
#include "ConicSection.h"
#include "EventsLoader.h"


using namespace std;
using namespace prompt_gamma_reconstruction;

TFile *EventsLoader::open_root_file_(){

    TFile *f = TFile::Open(data_file_path_.c_str(), "READ");
    if( 0 == f || !(f->IsOpen()) ){
        stringstream ss;
        ss << "PGEventsLoader::open_root_file_: failed to open " 
                <<  data_file_path_
                << ".\nAborting\n\n";
        cout<<ss.str()<<endl;
        throw runtime_error(ss.str());
    }

    return f;
}

TTree *EventsLoader::get_gamma_tree_(TFile &file){

    TTree *t;
    file.GetObject(gamma_tree_name_.c_str(), t);
    if( 0 == t){
        string error_message("ERROR: failed to get tree ");
        error_message += gamma_tree_name_;
        error_message += ".\nABORTING\n\n";
        cout<<error_message<<endl;
        throw runtime_error(error_message);
    }
    if( 0 == t->GetEntries()){
        string error_message("ERROR: tree ");
        error_message += gamma_tree_name_;
        error_message += " is empty.\nABORTING\n\n";
        cout<<error_message<<endl;
        throw runtime_error(error_message);
    }

  return t;
}
