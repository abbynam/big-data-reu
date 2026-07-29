#ifndef DCA_LINE_EVENTS_LOADER_H_
#define DCA_LINE_EVENTS_LOADER_H_
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
namespace prompt_gamma_reconstruction{

/*! \brief Read in list mode gamma camera events
 *
 * Reads in triple scatter events from CSV file and copies them into ConicSection or daughter objects.
 * Measured energies are replaced with energies from the know gamma spectra that produce the
 * smalles DCA to the line corresponding to the proton beam path.
 *
 * @author Dennis Mackin
 */
class DCALineEventsLoader: public EventsLoader {

  public:
    DCALineEventsLoader(const string &root_file_path, const pg_tools::RunTimeParameters *params, shared_ptr<const PhantomVolume> phantom);
    ~DCALineEventsLoader(){ /* */ };
    void LoadEvents(vector<shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point);

  private:
    pg_tools::RunTimeParameters const * run_time_parameters_ptr_;
    size_t read_CSV_file_into_vector_(const vector<vector<float> > &scattering_data, vector< shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point);
    shared_ptr<ConicSection> compareScatterOrderings_(shared_ptr<ConicSection> cs1, shared_ptr<ConicSection> cs2);
    void put_scatters_in_order_(vector<vector<float> > &data);
    vector<shared_ptr<ConicSection> > load_events_(vector<vector<float> > &data, const vector<float> &energies);
    pair<int, string> get_conic_section_(vector<float> &event, shared_ptr<ConicSection> &cs, long event_number);
    pair<int, string> dca_line_filter_(vector<shared_ptr<ConicSection> > &cs_vec, shared_ptr<ConicSection> &cs_selected);

    vector<vector<float> > apply_energy_filter_(vector<vector<float> >  &data);

    pair<int, string> filter_double_scatter_(shared_ptr<Scatter> &tmpScatter){
        vector< pair<int, string> (*)(shared_ptr<Scatter> tmpScatter, pg_tools::RunTimeParameters const * const params)> filters;
        filters.push_back(scatter_distance_filter);
        filters.push_back(nan_scattering_angle_filter);
        filters.push_back(energy_window_filter);
        filters.push_back(energy_lost_filter);
        filters.push_back(scattering_angle_filter);

        for(auto const &f: filters){
            auto result = f(tmpScatter, run_time_parameters_ptr_);
            if(result.first != 0) return result;
        }
        return std::make_pair(0, "double scatter passed filters");
    };

};



}
#endif //DCA_LINE_EVENTS_LOADER_H_
