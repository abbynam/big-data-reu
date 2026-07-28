#ifndef PGSIMULATED_EVENTS_LOADER_H_
#define PGSIMULATED_EVENTS_LOADER_H_

//standard C++ includes
#include <memory>

//PromptGamma includes
#include "ConicSection.h"
#include "PGVector3.h"
#include "EventsLoader.h"


using namespace std;
namespace prompt_gamma_reconstruction{


/*! \brief Reads in triple scatter events from PGC Geant4 based simulations.
 * 
 *  Reads in triple scatter events from PGC Geant4 based simulations.
 *  simulated and copies them into ConicSection or daughter objects.
 *  A vector of pointers is used to track
 *  the conic sections.
 * 
 * @author Dennis Mackin
 */
    
class PGSimulatedEventsLoader: public EventsLoader {

  public:

    PGSimulatedEventsLoader(const string &data_file_path, const pg_tools::RunTimeParameters *params, shared_ptr<const PhantomVolume>phantom);
    ~PGSimulatedEventsLoader(){
      cout<<"Destroying ~PGSimulatedEventsLoader . . ."<<endl;
    }
    void LoadEvents(vector<shared_ptr<ConicSection> > &conics, int tries_per_random);

  private:

    size_t read_tree_into_vector_(TTree &tree, vector<shared_ptr<ConicSection> > &conics, int number_tries_per_random_point);

};

}
#endif //PGSIMULATED_EVENTS_LOADER_
