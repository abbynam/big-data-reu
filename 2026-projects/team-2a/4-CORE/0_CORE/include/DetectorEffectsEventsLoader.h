#ifndef DETECTOR_EFFECTS_EVENTS_LOADER_H_
#define DETECTOR_EFFECTS_EVENTS_LOADER_H_
#define _USE_MATH_DEFINES

//standard C++ includes
#include <memory>

//ROOT includes
#include "TTree.h"

//PromptGamma includes
#include "ConicSection.h"
#include "PGVector3.h"
#include "EventsLoader.h"
#include "RunTimeParameters.h"
#include "TripleScatter.h"


using namespace std;
namespace prompt_gamma_reconstruction{
    
/*! \brief Reads in events from the Geant4 simulations designed
 * to simulate detector finite resolution effects.
 * 
 * The events are read in from an ROOT file. Then smeared based on the 
 * parameters read in from the runtime configuration file.
 * 
 * @author Dennis Mackin
 */    
class DetectorEffectsEventsLoader: public EventsLoader {

    public:

        DetectorEffectsEventsLoader(const string &data_file_path, const pg_tools::RunTimeParameters *params, shared_ptr<const PhantomVolume> phantom);
        ~DetectorEffectsEventsLoader(){
            cout<<"Destroying ~DetectorEffectsEventsLoader . . ."<<endl;
        }
        void LoadEvents(vector<shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point);

    private:
        size_t read_tree_into_vector_(TTree &tree, vector< shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point) const;
        void AddDetectorEffects(TripleScatter &ts) const;
        float get_weight(const Scatter &sc) const;

};

}
#endif //DETECTOR_EFFECTS_EVENTS_LOADER_
