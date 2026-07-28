#ifndef SIMULATED_CC_EVENTS_LOADER_H_
#define SIMULATED_CC_EVENTS_LOADER_H_
#define _USE_MATH_DEFINES

//standard C++ includes
#include <memory>


//PromptGamma includes
#include "ConicSection.h"
#include "PGVector3.h"
#include "EventsLoader.h"
#include "RunTimeParameters.h"
#include "TripleScatter.h"


using namespace std;
namespace prompt_gamma_reconstruction{

////////////////////////////////////////////////////////////
/// SimulatedCCEventsLoader
///
/// Reads in triple scatter events from Geant4 based compton camera simulator
/// and copies them into ConicSection or daughter objects.
/// A vector of pointers is used to track the conic sections.
///
/// This is a concrete class which implement the abstract class PGEventsLoader.
///
///
////////////////////////////////////////////////////////////
class SimulatedCCEventsLoader: public EventsLoader {

  public:

    SimulatedCCEventsLoader(const string &data_file_path, const pg_tools::RunTimeParameters *params, shared_ptr<const PhantomVolume> phantom);
    ~SimulatedCCEventsLoader(){
      cout<<"Destroying ~SimulatedCCEventsLoader . . ."<<endl;
    }
    void LoadEvents(vector<shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point);

  private:
//     pg_tools::RunTimeParameters* run_time_parameters_ptr_;

    size_t read_tree_into_vector_(TTree &tree, vector< shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point);
    void AddDetectorEffects(TripleScatter &ts);

};

}
#endif //SIMULATED_CC_EVENTS_LOADER_
