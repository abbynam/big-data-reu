#ifndef CSV_EVENTS_LOADER_H_
#define CSV_EVENTS_LOADER_H_
#define _USE_MATH_DEFINES

//standard C++ includes
#include <memory>

//PromptGamma includes
#include "ConicSection.h"
#include "PGVector3.h"
#include "EventsLoader.h"
#include "RunTimeParameters.h"
#include "TripleScatter.h"
#include "EventFilters.h"


using namespace std;
using namespace prompt_gamma_reconstruction;
namespace prompt_gamma_reconstruction{


/*! \brief Read in list mode gamma camera events
 *
 * Reads in triple scatter events from CSV file
 * and copies them into ConicSection or daughter objects.
 * A vector of pointers is used to track the conic sections.
 *
 * @author Dennis Mackin
 * @TODO set the parameters using a builder. Do not pass in the parameters file.
 */
class CSVEventsLoader: public EventsLoader {

  public:
    CSVEventsLoader(const string &root_file_path, const pg_tools::RunTimeParameters *params, shared_ptr<const PhantomVolume> phantom);
    ~CSVEventsLoader(){ /* */ };
    void LoadEvents(vector<shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point);

  private:
    pg_tools::RunTimeParameters const * run_time_parameters_ptr_;
    vector< pair<int, string> (*)(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const * const params)> filters_;
    size_t read_CSV_file_into_vector_(const string &file_path, vector< shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point);
    shared_ptr<ConicSection> compareScatterOrderings_(shared_ptr<ConicSection> &cs1, shared_ptr<ConicSection> &cs2);

    void set_filters(){
        filters_.push_back(scatter_distance_filter);
        filters_.push_back(nan_scattering_angle_filter);
        filters_.push_back(energy_window_filter);
        filters_.push_back(energy_lost_filter);
        filters_.push_back(scattering_angle_filter);
    }

    pair<int, string> filter_double_scatter_(shared_ptr<Scatter> tmpScatter){
        if(filters_.size() == 0) set_filters();

        for(auto const &f: filters_){
            auto result = f(tmpScatter, run_time_parameters_ptr_);
            if(result.first != 0) return result;
        }
        return std::make_pair(0, "double scatter passed filters");
    };

};



}
#endif //CSV_EVENTS_LOADER_H_
