#include "EventsLoaderFactory.h"
#include "CSVEventsLoader.h"
#include "DCALineEventsLoader.h"

using namespace prompt_gamma_reconstruction;

shared_ptr<EventsLoader> EventsLoaderFactory::create(const size_t event_file_type,
                                                     const string &data_file_path,
                                                     const pg_tools::RunTimeParameters *params,
                                                     shared_ptr<const PhantomVolume> phantom_ptr) {

    stringstream ss;
    switch (event_file_type) {
        case (1):
//        return shared_ptr<EventsLoader>(new PGSimulatedEventsLoader(data_file_path, params, phantom_ptr));
            ss << "ERROR: PGSimulatedEventsLoader is not supported in this version of CORE.";
            throw runtime_error(ss.str());
            break;
        case (2):
//        return shared_ptr<EventsLoader>(new SimulatedCCEventsLoader(data_file_path, params, phantom_ptr));
            ss << "ERROR: SimulatedCCEventsLoader is not supported in this version of CORE. ";
            throw runtime_error(ss.str());
            break;
        case (3):
            return shared_ptr<EventsLoader>(new CSVEventsLoader(data_file_path, params, phantom_ptr));
        case (4):
//        return shared_ptr<EventsLoader>(new DetectorEffectsEventsLoader(data_file_path, params, phantom_ptr))
            ss << "ERROR: DetectorEffectsEventsLoader is not supported in this version of CORE.";
            throw runtime_error(ss.str());
            break;
        case (5):
            return shared_ptr<EventsLoader>(new DCALineEventsLoader(data_file_path, params, phantom_ptr));
        default:
            ss << "ERROR: Invalid event file type " << event_file_type << ".";
            throw runtime_error(ss.str());
    };

};