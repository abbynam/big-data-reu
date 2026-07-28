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
#include "CSVEventsLoader.h"
#include "TripleScatter.h"
#include "utilities/FileUtils.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

CSVEventsLoader::CSVEventsLoader(const string &data_file_path, const pg_tools::RunTimeParameters *params, shared_ptr<const PhantomVolume>phantom):
                    EventsLoader(data_file_path, params, phantom), run_time_parameters_ptr_(params){
        /* DO NOTHING */
}


void CSVEventsLoader::LoadEvents(vector<shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point){
    read_CSV_file_into_vector_(data_file_path_, conics, number_tries_per_random_point);
}


size_t CSVEventsLoader::read_CSV_file_into_vector_(const string &file_path, vector< shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point) {

    vector< vector<float> > scattering_data;
    pg_tools::FileUtils::ReadDataFromFile(scattering_data, file_path,",");
    std::srand(23);
    std::random_shuffle(scattering_data.begin(), scattering_data.end());

    cout<<"read_CSV_file_into_vector_: number entries in file "<<scattering_data.size()<<" . . ."<<endl;

    //Counters for Scatting Events Characteristics
    size_t reconstructable_count = 0; //counter for number of reconstructable cones

    //number of attempts needed to find random point in phantom for all of the cones
    size_t num_cones_requested  = params_ptr_->get_int("MAX_NUM_CONES");

    float energy_deposited[3];
    size_t numbers_records = scattering_data.size();

    shared_ptr<ConicSection> tmpConicSection;
    shared_ptr<ComptonScatter> tmpComptonScatter;
    shared_ptr<Scatter> tmpScatter;
    map<string, long> msg_map;
    for(size_t i=0; i < numbers_records; ++i){
        
        if(scattering_data[i].size() == 12){ //Triple Scatter
            energy_deposited[0] = scattering_data[i][0];
            energy_deposited[1] = scattering_data[i][4];
            energy_deposited[2] = scattering_data[i][8];

            auto p1 = PGVector3(scattering_data[i][1], scattering_data[i][2], scattering_data[i][3]);
            auto p2 = PGVector3(scattering_data[i][5], scattering_data[i][6], scattering_data[i][7]);
            auto p3 = PGVector3(scattering_data[i][9], scattering_data[i][10], scattering_data[i][11]);

            tmpScatter = shared_ptr<Scatter> (new TripleScatter(energy_deposited[0], energy_deposited[1], energy_deposited[2], p1, p2, p3));
            auto msg = filter_double_scatter_(tmpScatter);

            if(msg.first != 0) {
                msg_map[msg.second]++;
                continue;
            }

            auto cs = shared_ptr<ComptonScatter> (new ComptonScatter(*tmpScatter));
            tmpConicSection = shared_ptr<ConicSection>( new ConicSection( *cs, phantom_volume_, i));
            tmpConicSection->setInverseSquareParam(run_time_parameters_ptr_->get_float("INVERSE_SQUARE_PARAM"));
            tmpConicSection->setScatter(tmpScatter);
            msg = find_random_in_phantom_filter(tmpConicSection, run_time_parameters_ptr_);
            if(msg.first != 0) {
                msg_map[msg.second]++;
                continue;
            }

        }else if(scattering_data[i].size() == 8 || scattering_data[i].size() == 9) { //Double Scatter

            shared_ptr<ConicSection> tmpConicSection1;
            shared_ptr<ConicSection> tmpConicSection2;

            energy_deposited[0] = scattering_data[i][0];
            energy_deposited[1] = scattering_data[i][4];
            energy_deposited[2] = 0.0;

            auto p1 = PGVector3(scattering_data[i][1], scattering_data[i][2], scattering_data[i][3]);
            auto p2 = PGVector3(scattering_data[i][5], scattering_data[i][6], scattering_data[i][7]);

            auto tmpScatter1 = shared_ptr<Scatter>(new DoubleScatter(energy_deposited[0], energy_deposited[1], p1, p2));
            auto tmpScatter2 = shared_ptr<Scatter>(new DoubleScatter(energy_deposited[1], energy_deposited[0], p2, p1));

            auto msg1 = filter_double_scatter_(tmpScatter1);
            if(msg1.first == 0){
                auto cs = shared_ptr<ComptonScatter> (new ComptonScatter(*tmpScatter1));
                tmpConicSection1 = shared_ptr<ConicSection>( new ConicSection( *cs, phantom_volume_, i+1));
                tmpConicSection1->setScatter(tmpScatter1);
                msg1 = DCA_filter(tmpConicSection1, run_time_parameters_ptr_);
//                if(msg1.first == 0) msg1 = find_random_in_phantom_filter(tmpConicSection1, run_time_parameters_ptr_);
            }

            auto msg2 = filter_double_scatter_(tmpScatter2);
            if(msg2.first == 0){
                auto cs = shared_ptr<ComptonScatter> (new ComptonScatter(*tmpScatter2));

                tmpConicSection2 = shared_ptr<ConicSection>( new ConicSection(*cs, phantom_volume_, i+1));
                tmpConicSection2->setScatter(tmpScatter2);

                msg2 = DCA_filter(tmpConicSection2, run_time_parameters_ptr_);
//                if(msg2.first == 0) msg2 = find_random_in_phantom_filter(tmpConicSection2, run_time_parameters_ptr_);
            }

            if (msg1.first == 0 && msg2.first == 0) {
                msg_map["Both orderings work"]++;

                //DSM I changed it so that CORE does not override the ordering it is given unless
                // the ordering is non-physical. This ways we can test ordering from outside of CORE.
                //DSM I changed it back 2021-01-17
//                tmpConicSection = tmpConicSection1;
                tmpConicSection = compareScatterOrderings_(tmpConicSection1, tmpConicSection2);
//                cout<<"dca,both,"<<tmpConicSection->getDistanceToPoint(PGVector3(0,0,0))<<endl;
//                continue;
            }else if(msg1.first == 0){
                tmpConicSection = tmpConicSection1;
//                cout<<"dca,first,"<<tmpConicSection->getDistanceToPoint(PGVector3(0,0,0))<<endl;
//                continue;
            }else if(msg2.first == 0){
                tmpConicSection = tmpConicSection2;
//                cout<<"dca,second,"<<tmpConicSection->getDistanceToPoint(PGVector3(0,0,0))<<endl;
            }else{
                msg_map[msg1.second]++;
                msg_map[msg2.second]++;
                msg_map["Both orderings failed"]++;
                continue;
            }


        }else{ //Invalid record
            stringstream ss;
            ss<<"ERROR: invalid record length of "<< scattering_data[i].size() <<"\n";
            cout<<ss.str()<<endl;
            throw runtime_error(ss.str());
        }

        conics.push_back(tmpConicSection);

        msg_map["Events kept"]++;
        reconstructable_count++;

        if(num_cones_requested <= reconstructable_count){
            printf("%d cones, last cone loaded is %lu . . .\n", (int)reconstructable_count, i);
            break;
        }

        // Monitor progress
        size_t increment = (numbers_records / 10);
        if( 0 == increment) ++increment;
        if (i % increment == 0) {
            printf("    Row: %lu\n", i);
        }
    }

    printf("\n-----------Results-----------\n");
    for(auto msg : msg_map)
    {
        std::cout << msg.first << ": " << msg.second << "\n";
    }

    return reconstructable_count;
}


shared_ptr<ConicSection> CSVEventsLoader::compareScatterOrderings_(shared_ptr<ConicSection> &cs1, shared_ptr<ConicSection> &cs2){
    ///If the total energy is less than 0.7 MeV, then make the the highest energy scatter the first scatter.
    // If the total energy is greater than 0.7MeV, then choose the scatter closest to the DCA Center as the first scatter.
    // The DCA center should usually be isocenter.

    auto E0 = cs2->getScatterInfo()->getGammaEnergy();
    PGVector3 p0;
    p0.x = run_time_parameters_ptr_->get_float("DCA_CENTER_X");
    p0.y = run_time_parameters_ptr_->get_float("DCA_CENTER_Y");
    p0.z = run_time_parameters_ptr_->get_float("DCA_CENTER_Z");

    auto distance_calc = [](PGVector3 a, PGVector3 b){
        return fabs(sqrt( (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) + (a.z - b.z)* (a.z - b.z)));
    };
    if(E0 > 0.7){
        if(distance_calc(p0, cs1->getScatterInfo()->getConeApex()) < distance_calc(p0, cs2->getScatterInfo()->getConeApex())){
            return cs1;
        } else {
            // cout<<"Switching the event ordering (LittleD). . ."<< endl;
            return cs2;
        }
    }

    if(cs1->getScatterInfo()->getScatter1EnergyDeposit() > cs2->getScatterInfo()->getScatter1EnergyDeposit()){
        return cs1;
    } else {
        // cout<<"Switching the event ordering (Big E) . . ."<< endl;
        return cs2;
    }

};

