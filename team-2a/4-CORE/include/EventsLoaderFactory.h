#ifndef EVENTS_LOADER_FACTORY_H_
#define EVENTS_LOADER_FACTORY_H_

//standard C++ includes
#include <valarray>
#include <vector>
#include <algorithm>
#include <memory>

//local package includes
#include "EventsLoader.h"
//#include "PGSimulatedEventsLoader.h"
//#include "SimulatedCCEventsLoader.h"
#include "CSVEventsLoader.h"
//#include "DetectorEffectsEventsLoader.h"

/*! \brief Creates the EventsLoader based on the type of file
 * used for input. They is specified by the <kbd>DATA_FILE_FORMAT</kbd> 
 * parameter.
 * 
 * @author Dennis Mackin
 */
namespace prompt_gamma_reconstruction{
    class EventsLoaderFactory{
        
    public:
        static shared_ptr<EventsLoader> create(const size_t event_file_type, const string &data_file_path,
                                               const pg_tools::RunTimeParameters *params,
                                               shared_ptr<const PhantomVolume> phantom_ptr);
    };
}

#endif //EVENTS_LOADER_FACTORY_H_