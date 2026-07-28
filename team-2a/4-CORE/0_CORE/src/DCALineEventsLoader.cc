#include "DCALineEventsLoader.h"

// Standard C++ Includes
#include <cmath>
#include <valarray>
#include <vector>
#include <stdexcept>
#include <iostream>


//PG Includes
#include "PhantomVolume.h"
#include "ReconstructionParabola.h"
#include "ReconstructionEllipse.h"
#include "ConicSection.h"
#include "EventsLoader.h"
#include "DoubleScatter.h"
#include "TripleScatter.h"
#include "utilities/FileUtils.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

DCALineEventsLoader::DCALineEventsLoader(const string &data_file_path, const pg_tools::RunTimeParameters *params,
                                         shared_ptr<const PhantomVolume> phantom) :
        EventsLoader(data_file_path, params, phantom), run_time_parameters_ptr_(params) {
    /* DO NOTHING */
}


void DCALineEventsLoader::LoadEvents(vector<shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point) {
    vector<vector<float> > scatter_data;

    cout<<"Reading in the CSV file . . ."<< endl;
    pg_tools::FileUtils::ReadDataFromFile(scatter_data, data_file_path_, ",");

    cout<<"Remove scatters based on energy deposits . . ." <<endl;
    scatter_data = apply_energy_filter_(scatter_data);

    cout<<"Put scatters in order . . ."<<endl;
    put_scatters_in_order_(scatter_data);

    cout<<"Build conic section objects . . ."<<endl;
    auto gamma_energies = params_ptr_->get_csv_values("KNOWN_GAMMA_ENERGIES");
    conics = load_events_(scatter_data, gamma_energies);

    cout<<"Shuffling the events . . . "<<endl;
    std::srand(42);
    std::random_shuffle(conics.begin(), conics.end());

    auto max_number_events = this->params_ptr_->get_int("MAX_NUM_CONES");
    cout<<"Selecting first "<< max_number_events <<" events . . ."<< endl;
    if(conics.size() > (size_t)max_number_events) conics.resize(max_number_events);
//    cout<<"Creating the list of conic sections from the scatters . . ."<<endl;
//    read_CSV_file_into_vector_(scatter_data, conics, number_tries_per_random_point);
}

vector<vector<float> > DCALineEventsLoader::apply_energy_filter_(vector<vector<float> >  &data){

    const float MIN_ENERGY_SCATTER = this->params_ptr_->get_float("MIN_ENERGY_SCATTER");
    const float MIN_ENERGY_EVENT = this->params_ptr_->get_float("MIN_ENERGY_EVENT");

    vector<vector<float> > selected_events;
    selected_events.reserve(data.size());
    for(size_t i=0; i<data.size(); ++i){
        auto total_energy = data[i][0] + data[i][4];
        auto min_scatter = min(data[i][0], data[i][4]);
        if(min_scatter > MIN_ENERGY_SCATTER && total_energy > MIN_ENERGY_EVENT){
            selected_events.push_back(data[i]);
        }
    }

    cout<<"Energy threshold (" << MIN_ENERGY_SCATTER <<", "<< MIN_ENERGY_EVENT << "): " << selected_events.size() << " of " << data.size() << " passed . . ."<<endl;
    return selected_events;
}

vector<shared_ptr<ConicSection> > DCALineEventsLoader::load_events_(vector<vector<float> > &data, const vector<float> &energies){

    map<string, long> msg_map;
    auto add_message = [&](float e, string m){ stringstream ss; ss << e <<": "<< m; msg_map[ss.str()]++; };

    auto adjusted_event = data[0];
    vector<shared_ptr<ConicSection> > conic_sections;
    vector<shared_ptr<ConicSection> > conic_sections_to_compare;
    shared_ptr<ConicSection> cs;
    pair<int, string> msg;

    for(size_t i=0; i<data.size(); ++i){
        for(size_t j=0; j<energies.size(); ++j){
            adjusted_event = data[i];
            adjusted_event[4] = energies[j] - adjusted_event[0];

            msg = get_conic_section_(adjusted_event, cs, i);
            add_message(energies[j], msg.second);
            if(msg.first == 0) {
                conic_sections_to_compare.push_back(cs);
            }
        }
        if(conic_sections_to_compare.size() > 0){
            msg = dca_line_filter_(conic_sections_to_compare, cs);
            if(msg.first == 0) {
//                #pragma omp critical
                { conic_sections.push_back(cs); }
                add_message(cs->getComptonScatter().getInitialEnergy(), "selected");
            }else{
                add_message(0.0, "All tested energies failed DCA line cut.");
            }
        }
        conic_sections_to_compare.clear();
    }

    printf("\n-----------Results-----------\n");
    for (auto msg : msg_map) {
        std::cout << msg.first << ": " << msg.second << "\n";
    }

    return conic_sections;
}

/// dca_line_filter_ - Selects the ConicSection with known energy that comes closest to the line corresponding to the
///                 proton beam path.
/// @param cs_vec - vector of conic sections for each known energy that survived the initial filters
/// @param cs_select - shared pointer to the ConicSection that has the smallest DCA to the line of the proton beam.
///
/// If multiple energies produce ConicSections that intersect the line producing DCA 0, then the largest energy
/// will be selected since it will have the smallest scattering angle and, therefore, the smalllest delta theta.
pair<int, string> DCALineEventsLoader::dca_line_filter_(vector<shared_ptr<ConicSection> > &cs_vec, shared_ptr<ConicSection> &cs_selected){

    auto p = params_ptr_->get_csv_values("BEAM_LINE_POINT1");
    PGVector3 p1(p[0], p[1], p[2]);
    assert(p.size() == 3);
    p = params_ptr_->get_csv_values("BEAM_LINE_POINT2");
    PGVector3 p2(p[0], p[1], p[2]);
    assert(p.size() == 3);
    const float DCA_LINE_CUT = params_ptr_->get_float("DCA_LINE_CUT");

    pair<float, PGVector3> results;
    float min_dca = 100.0*DCA_LINE_CUT;

    cs_selected = nullptr;
    for(size_t i=0; i < cs_vec.size(); ++i){
        results = cs_vec[i]->getDistanceToLine(p1, p2);

        //If DCA is less then 1 mm, set it equal to 0.0 for comparison purposes
        if(results.first < 0.5) results.first = 0.0;
        if(results.first <= min_dca){
            cs_selected = cs_vec[i];
            min_dca = results.first;
        }
    }

    if(min_dca < DCA_LINE_CUT){
        return make_pair(0, "Passed DCA line filter.");
    }
//    cout<<"Failed DCA line filter . . ."<<endl;
    return make_pair(1000, "Failed DCA line filter.");
}


pair<int, string> DCALineEventsLoader::get_conic_section_(vector<float> &event, shared_ptr<ConicSection> &cs, long event_number){

    auto p1 = PGVector3(event[1], event[2], event[3]);
    auto p2 = PGVector3(event[5], event[6], event[7]);

    auto tmpScatter = shared_ptr<Scatter>(new DoubleScatter(event[0], event[4], p1, p2));

    auto msg = filter_double_scatter_(tmpScatter);
    if(msg.first != 0) return msg;

    auto compton_scatter = shared_ptr<ComptonScatter>(new ComptonScatter(*tmpScatter));
    cs = shared_ptr<ConicSection>(new ConicSection(*compton_scatter, phantom_volume_, event_number));
    cs->setScatter(tmpScatter);

    msg = DCA_filter(cs, run_time_parameters_ptr_);
    if (msg.first != 0) return msg;

    msg = find_random_in_phantom_filter(cs, run_time_parameters_ptr_);
    if (msg.first != 0) return msg;

    return make_pair(0, "event accepted.");
};


void DCALineEventsLoader::put_scatters_in_order_(vector<vector<float> > &data){
    vector<float> tmp; //buffer for reordering event data
    auto num_records = data.size();
    for(size_t i = 0; i < num_records; ++i){
        tmp = data[i];
        //validate input
        if(data[i].size() % 4 != 0){
            cerr<<"ERROR: Invalid event line " << i <<".\n"<<endl;
            for_each(data[i].begin(), data[i].end(), [&data](float x){cerr<<" "<<x<<",";});
            cerr<<endl;
        }

        //set pointers to start of each scatter
        vector<float const *> scatter_ptrs(data[i].size() / 4);
        for(size_t j=0; j< scatter_ptrs.size(); ++j){
            scatter_ptrs[j] = &tmp[4*j];
        }

        //sort the points based on distance from isocenter
        sort(scatter_ptrs.begin(), scatter_ptrs.end(), [&](float const a[], float const b[]){
            return (a[1]*a[1] + a[2]*a[2] + a[3]*a[3] < b[1]*b[1] + b[2]*b[2] + b[3]*b[3]);
        });

        //point the data into order in temporary array
        for(size_t j=0; j< scatter_ptrs.size(); ++j){
            for(size_t k = 0; k < 4; ++k) data[i][j*4 + k] = scatter_ptrs[j][k];
        }

    }
}


size_t DCALineEventsLoader::read_CSV_file_into_vector_(const vector<vector<float> > &scattering_data,
                                                       vector<shared_ptr<ConicSection> > &conics,
                                                       size_t number_tries_per_random_point) {

    cout << "read_CSV_file_into_vector_: number entries in file " << scattering_data.size() << " . . ." << endl;

    size_t reconstructable_count = 0; //counter for number of reconstructable cones
    size_t num_cones_requested = params_ptr_->get_int("MAX_NUM_CONES");

    float energy_deposited[3];
    size_t numbers_records = scattering_data.size();
    shared_ptr<ConicSection> tmpConicSection1;
    shared_ptr<ConicSection> tmpConicSection2;
    shared_ptr<ConicSection> tmpConicSection;
    shared_ptr<ComptonScatter> tmpComptonScatter;
    shared_ptr<Scatter> tmpScatter;
    map<string, long> msg_map;
    for (size_t i = 0; i < numbers_records; ++i) {

        if (scattering_data[i].size() == 12) { //Triple Scatter
            energy_deposited[0] = scattering_data[i][0];
            energy_deposited[1] = scattering_data[i][4];
            energy_deposited[2] = scattering_data[i][8];

            auto p1 = PGVector3(scattering_data[i][1], scattering_data[i][2], scattering_data[i][3]);
            auto p2 = PGVector3(scattering_data[i][5], scattering_data[i][6], scattering_data[i][7]);
            auto p3 = PGVector3(scattering_data[i][9], scattering_data[i][10], scattering_data[i][11]);

            tmpScatter = shared_ptr<Scatter>(
                    new TripleScatter(energy_deposited[0], energy_deposited[1], energy_deposited[2], p1, p2, p3));
            auto cs = shared_ptr<ComptonScatter>(new ComptonScatter(*tmpScatter));
            tmpConicSection = shared_ptr<ConicSection>(new ConicSection(*cs, phantom_volume_, i));
            tmpConicSection->setInitialOrigin(run_time_parameters_ptr_->get_float("INVERSE_SQUARE_PARAM"));
            tmpConicSection->setScatter(tmpScatter);

        } else if (scattering_data[i].size() == 8) { //Double Scatter

            energy_deposited[0] = scattering_data[i][0];
            energy_deposited[1] = scattering_data[i][4];
            energy_deposited[2] = 0.0;

            auto p1 = PGVector3(scattering_data[i][1], scattering_data[i][2], scattering_data[i][3]);
            auto p2 = PGVector3(scattering_data[i][5], scattering_data[i][6], scattering_data[i][7]);

            auto tmpScatter1 = shared_ptr<Scatter>(new DoubleScatter(energy_deposited[0], energy_deposited[1], p1, p2));
            auto tmpScatter2 = shared_ptr<Scatter>(new DoubleScatter(energy_deposited[1], energy_deposited[0], p2, p1));

            auto msg1 = filter_double_scatter_(tmpScatter1);
            if (msg1.first == 0) {
                auto cs = shared_ptr<ComptonScatter>(new ComptonScatter(*tmpScatter1));
                tmpConicSection1 = shared_ptr<ConicSection>(new ConicSection(*cs, phantom_volume_, i));
                tmpConicSection1->setScatter(tmpScatter1);
                msg1 = DCA_filter(tmpConicSection1, run_time_parameters_ptr_);
                if (msg1.first == 0) msg1 = find_random_in_phantom_filter(tmpConicSection1, run_time_parameters_ptr_);
            }

            auto msg2 = filter_double_scatter_(tmpScatter2);
            if (msg2.first == 0) {
                auto cs = shared_ptr<ComptonScatter>(new ComptonScatter(*tmpScatter2));

                tmpConicSection2 = shared_ptr<ConicSection>(new ConicSection(*cs, phantom_volume_, i));
                tmpConicSection2->setScatter(tmpScatter2);

                msg2 = DCA_filter(tmpConicSection2, run_time_parameters_ptr_);
                if (msg2.first == 0) msg2 = find_random_in_phantom_filter(tmpConicSection2, run_time_parameters_ptr_);
            }

            if (msg1.first == 0 && msg2.first == 0) {
                msg_map["Both orderings work"]++;
//                tmpConicSection = compareScatterOrderings_(tmpConicSection1, tmpConicSection2);
                tmpConicSection = tmpConicSection1;
            } else if (msg1.first == 0) {
                tmpConicSection = tmpConicSection1;
            } else if (msg2.first == 0) {
                tmpConicSection = tmpConicSection2;
            } else {
                msg_map[msg1.second]++;
                msg_map[msg2.second]++;
                msg_map["Both orderings failed"]++;
                continue;
            }


        } else { //Invalid record

            stringstream ss;
            ss << "ERROR: invalid record length of " << scattering_data[i].size() << "\n";
            cout << ss.str() << endl;
            throw runtime_error(ss.str());
        }

        conics.push_back(tmpConicSection);
        msg_map["Events kept"]++;
        reconstructable_count++;

        if (num_cones_requested <= reconstructable_count) {
            printf("%d cones, last cone loaded is %lu . . .\n", (int) reconstructable_count, i);
            break;
        }

        // Monitor progress
        size_t increment = (numbers_records / 10);
        if (0 == increment) ++increment;
        if (i % increment == 0) {
            printf("    Row Number: %lu\n", i);
        }
    }

    printf("\n-----------Results-----------\n");
    for (auto msg : msg_map) {
        std::cout << msg.first << ": " << msg.second << "\n";
    }

    return reconstructable_count;
}


shared_ptr<ConicSection> DCALineEventsLoader::compareScatterOrderings_(shared_ptr<ConicSection> cs1,
                                                                       shared_ptr<ConicSection> cs2) {

    PGVector3 dca_center_point;
    dca_center_point.x = run_time_parameters_ptr_->get_float("DCA_CENTER_X");
    dca_center_point.y = run_time_parameters_ptr_->get_float("DCA_CENTER_Y");
    dca_center_point.z = run_time_parameters_ptr_->get_float("DCA_CENTER_Z");

    auto dca1 = cs1->getDistanceToPoint(dca_center_point);
    auto dca2 = cs2->getDistanceToPoint(dca_center_point);
    if (dca1 < dca2) {
        return cs1;
    } else {
        return cs2;
    }
};

