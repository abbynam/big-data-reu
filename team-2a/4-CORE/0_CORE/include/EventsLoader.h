#ifndef EVENTS_LOADER_H_
#define EVENTS_LOADER_H_
#define _USE_MATH_DEFINES

//standard C++ includes
#include <memory>

//ROOT includes
//#include "TTree.h"
//PromptGamma includes
#include "ConicSection.h"
#include "PGVector3.h"
#include "RunTimeParameters.h"


using namespace std;
namespace prompt_gamma_reconstruction{

    
/*! \brief Base class for reading in triple scatter events,
 * 
 * Concrete subclasses should include at least a class 
 * to  read in events from the Geant4 simulation and 
 * one to load in events from a real detector.
 * 
 * @author Dennis Mackin
 */    
class EventsLoader{

  public:

    EventsLoader(const string &data_file_path, const pg_tools::RunTimeParameters *params, shared_ptr<const PhantomVolume>phantom):
        data_file_path_(data_file_path), params_ptr_(params), phantom_volume_(phantom) { /* DO NOTHING */ };

    virtual ~EventsLoader(){
        cout<<"Destroying ~PGEventsLoader . . ."<<endl;
    }
    virtual void LoadEvents(vector<shared_ptr<ConicSection> > &conics, size_t tries_per_random)=0;

  protected:
    string data_file_path_;
    string gamma_tree_name_;
    const pg_tools::RunTimeParameters *params_ptr_;
    shared_ptr<const PhantomVolume> phantom_volume_;

//    virtual TTree *get_gamma_tree_(TFile &file);
//    virtual TFile *open_root_file_();
};

};
#endif //EVENTS_LOADER_
