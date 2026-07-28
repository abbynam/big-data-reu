#ifndef EVENT_FILTERS_H_
#define EVENT_FILTERS_H_
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

    pair<int, string> scatter_distance_filter(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const *params);
    pair<int, string> nan_scattering_angle_filter(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const *params);
    pair<int, string> energy_window_filter(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const *params);
    pair<int, string> energy_lost_filter(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const *params);
    pair<int, string> scattering_angle_filter(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const *params);
    pair<int, string> DCA_filter(shared_ptr<ConicSection> cs, pg_tools::RunTimeParameters const *params);
    pair<int, string> find_random_in_phantom_filter(shared_ptr<ConicSection> cs, pg_tools::RunTimeParameters const *params);

}
#endif //EVENT_FILTERS_H_
