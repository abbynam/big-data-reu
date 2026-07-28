#include "SimulatedCCEventsLoader.h"

// Standard C++ Includes
#include <cmath>
#include <valarray>
#include <vector>
#include <stdexcept>
#include <iostream>


// ROOT Includes
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"

//PG Includes
#include "PhantomVolume.h"
#include "ReconstructionParabola.h"
#include "ReconstructionEllipse.h"
#include "ConicSection.h"
#include "EventsLoader.h"
#include "TripleScatter.h"


using namespace std;
using namespace prompt_gamma_reconstruction;

SimulatedCCEventsLoader::SimulatedCCEventsLoader(const string &data_file_path, const pg_tools::RunTimeParameters *params, shared_ptr<const PhantomVolume>phantom):
    EventsLoader(data_file_path, params, phantom){
    gamma_tree_name_ = (*params)["GAMMA_TREE_NAME"];
}


void SimulatedCCEventsLoader::AddDetectorEffects(TripleScatter &ts){
    PGVector3 positionUnc1, positionUnc2, positionUnc3;
    positionUnc1.x =  strtod((*params_ptr_)["DETECTOR1_X_UNC"].c_str(),0);
    positionUnc1.y =  strtod((*params_ptr_)["DETECTOR1_Y_UNC"].c_str(),0);
    positionUnc1.z =  strtod((*params_ptr_)["DETECTOR1_Z_UNC"].c_str(),0);
    ts.applyPositionUncertainty(positionUnc1, 0);

    positionUnc2.x =  strtod((*params_ptr_)["DETECTOR2_X_UNC"].c_str(),0);
    positionUnc2.y =  strtod((*params_ptr_)["DETECTOR2_Y_UNC"].c_str(),0);
    positionUnc2.z =  strtod((*params_ptr_)["DETECTOR2_Z_UNC"].c_str(),0);
    ts.applyPositionUncertainty(positionUnc2, 1);

    positionUnc3.x =  strtod((*params_ptr_)["DETECTOR3_X_UNC"].c_str(),0);
    positionUnc3.y =  strtod((*params_ptr_)["DETECTOR3_Y_UNC"].c_str(),0);
    positionUnc3.z =  strtod((*params_ptr_)["DETECTOR3_Z_UNC"].c_str(),0);
    ts.applyPositionUncertainty(positionUnc3, 2);

    //energy uncertainty uses the formula dE(E) = sqrt(alpha +beta*E)
    string detector_type = (*params_ptr_)["DETECTOR1_MATERIAL"];
    float e_scalar =  strtod((*params_ptr_)["DETECTOR1_E_SCALAR"].c_str(),0);

    ts.applyEnergyUncertainty(detector_type, e_scalar, 0);

    detector_type = (*params_ptr_)["DETECTOR2_MATERIAL"];
    e_scalar =  strtod((*params_ptr_)["DETECTOR2_E_SCALAR"].c_str(),0);
    ts.applyEnergyUncertainty(detector_type, e_scalar, 1);

}

void SimulatedCCEventsLoader::LoadEvents(vector<shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point){

    cout<<"Reading events from "<<data_file_path_<<" . . ."<<endl;
    TFile *f = open_root_file_();
    TTree *gamma_tree = get_gamma_tree_(*f);

    cout<<"number entries in tree "<<gamma_tree->GetEntries()<<" . . ."<<endl;

    // calculate cone data, returning cone data and number of cones calculated
    size_t num_cones = read_tree_into_vector_(*gamma_tree, conics, number_tries_per_random_point);
    cout<<"read  "<<num_cones<<" from tree . . ."<<endl;
    delete gamma_tree;
    f->Close();
    delete f;
}


size_t SimulatedCCEventsLoader::read_tree_into_vector_(TTree &tree, vector< shared_ptr<ConicSection> > &conics, size_t number_tries_per_random_point) {

    // Get number of TTree entries
    size_t rows = (int)tree.GetEntries();
    cout<<"read_tree_into_vector_: number entries in tree "<<tree.GetEntries()<<" . . ."<<endl;

    struct EventData{
        size_t event_num;
        PGVector3 cone_apex;
        PGVector3 cone_axis;
        PGVector3 true_origin;
        float x[3];
        float y[3];
        float z[3];
        float origin_x;
        float origin_y;
        float origin_z;
        float scatAng[3];
        float energy_deposited[3];
        float energy_incident[3];
        float energy_initial;
    } ed;

    tree.SetBranchAddress("event", &ed.event_num);

    tree.SetBranchAddress("pos_x", ed.x);
    tree.SetBranchAddress("pos_y", ed.y);
    tree.SetBranchAddress("pos_z", ed.z);
    tree.SetBranchAddress("incident_energy", &ed.energy_incident);
    ed.energy_initial = ed.energy_incident[0];

    tree.SetBranchAddress("scatAng", ed.scatAng);
    tree.SetBranchAddress("energyDeposited", ed.energy_deposited);

    tree.SetBranchAddress("origin_x", &ed.origin_x);
    tree.SetBranchAddress("origin_y", &ed.origin_y);
    tree.SetBranchAddress("origin_z", &ed.origin_z);
    ed.true_origin = PGVector3(ed.origin_x,ed.origin_y,ed.origin_z);

    //   tree.SetBranchAddress("origin_energy", ed.energy_true);

    PGVector3 p[3];

    PGVector3 random_point(phantom_volume_->x_max+1000.0,phantom_volume_->y_max+1000.0,phantom_volume_->z_max+1000.0);
    size_t reconstructable_count = 0; //counter for number of reconstructable
    size_t number_ellipses = 0; //number cones intersecting as ellipse
    size_t number_parabolas = 0; //number cones intersecting as parabolas
    size_t number_cones_skipped = 0; //count of number of un-reconstructable conics
    size_t nan_skipped=0;
    size_t dca_cut_skipped=0;
    size_t not_in_phantom_skipped=0;
    size_t randoms_not_in_cone_skipped=0;

    //number of attempts needed to find random point in phantom
    // for all of the cones
    size_t number_random_points_tried = 0;
    size_t max_num_cones  = params_ptr_->get_int("MAX_NUM_CONES");
    size_t offset  = params_ptr_->get_int("NUM_CONES_OFFSET");
    cout<<"GAMMA OFFSET: "<<offset<<", "<<max_num_cones<<" cones requested, "<<tree.GetEntries()<<" cones total . . ."<<endl;
    PGVector3 dca_cut_point;
    dca_cut_point.x = params_ptr_->get_float("DCA_CUT_X");
    dca_cut_point.y = params_ptr_->get_float("DCA_CUT_Y");
    dca_cut_point.z = params_ptr_->get_float("DCA_CUT_Z");
    float dca_cut =   params_ptr_->get_float("DCA_CUT");

    for(size_t i=offset; i< rows ; ++i){
        tree.GetEntry(i);
        p[0] = PGVector3(ed.x[0],ed.y[0], ed.z[0]);
        p[1] = PGVector3(ed.x[1], ed.y[1], ed.z[1]);
        p[2] = PGVector3(ed.x[2], ed.y[2], ed.z[2]);
        ed.true_origin = PGVector3(ed.origin_x,ed.origin_y,ed.origin_z);

        shared_ptr<TripleScatter> ts(new TripleScatter(ed.energy_deposited[0], ed.energy_deposited[1], ed.energy_deposited[2], p[0], p[1], p[2]));

        AddDetectorEffects(*ts);
        double scattering_angle = ts->getConeOpeningAngle();
        if( scattering_angle != scattering_angle){ //scattering_angle is not a number, skip it
            ++number_cones_skipped;
            ++nan_skipped;
            continue;
        }
        shared_ptr<ComptonScatter> tmpComptonScatter;
        try{
          tmpComptonScatter = shared_ptr<ComptonScatter> (new ComptonScatter(
                                                            ts->getConeApex(),
                                                            ts->getConeAxis(),
                                                            ed.true_origin,
                                                            ts->getConeOpeningAngle(),
                                                            ts->getScatter1EnergyDeposit(),
                                                            ts->getScatter2EnergyDeposit(),
                                                            ts->getGammaEnergy()
                                                          ));
        }catch( runtime_error e ){
            std::cout<<e.what()<<std::endl;
            ++number_cones_skipped;
            continue;
        }


        double alpha = tmpComptonScatter->getAlpha();
        shared_ptr<ConicSection> ptrConicSection;

        if( alpha + ts->getConeOpeningAngle() > M_PI/2.0){//then conic section is parabola
            ptrConicSection = shared_ptr<ConicSection>( new ReconstructionParabola( *tmpComptonScatter, phantom_volume_,i));
            ++number_parabolas;
        }else{//then conic section is an ellipse
            ptrConicSection = shared_ptr<ConicSection>( new ReconstructionEllipse( *tmpComptonScatter, phantom_volume_,i));
            ++number_ellipses;
        }

        if(ptrConicSection->getDistanceToPoint(dca_cut_point) > dca_cut){
            ++dca_cut_skipped;
            continue;
        }   

        size_t number_tries = 0; //number of randoms thrown before getting a point in phantom
        number_random_points_tried += number_tries;

        ptrConicSection->setLikelyOrigin(random_point);
        ptrConicSection->setMCTruth(ed.scatAng, p, ed.energy_deposited, ed.energy_incident, ed.energy_initial, &ed.origin_x);
        ptrConicSection->setScatter(ts);
        conics.push_back(ptrConicSection);

        ++reconstructable_count;
        if(max_num_cones <= reconstructable_count){
            printf("%d cones, last cone loaded is %d . . .\n", reconstructable_count, i);
            break;
        }

        // Monitor progress
        size_t increment = rows / 10;
        if (i % increment == 0) {
            printf("    Row Number: %d, Event Number: %i, Reconstructable: %d\n", i, ed.event_num, reconstructable_count);
        }
    }
    printf("--- Number of parabolas in TTree: %d ---\n", number_parabolas);
    printf("--- Number of ellipses in TTree: %d ---\n", number_ellipses);
    printf("--- Number of skipped due to DCA CUT: %d ---\n", dca_cut_skipped);
    printf("--- Number of skipped due to nan angle: %d ---\n", nan_skipped);
    printf("--- Number of skipped due to no cone phantom intercept: %d ---\n", not_in_phantom_skipped);

    printf("--- Number of skipped due to no randoms in the phantom: %d ---\n", randoms_not_in_cone_skipped); 


    printf("--- Number of useable triple scatters in TTree: %d ---\n", reconstructable_count);

    if(0 < reconstructable_count){
        printf("--- Avg. randoms points tested per conic section: %.2f ---\n",
        static_cast<float>(number_random_points_tried)/static_cast<float>(reconstructable_count));
    }
    return reconstructable_count;
}
