#ifndef _RECONSTRUCTOR_BUILDER
#define _RECONSTRUCTOR_BUILDER

//Standard includes
#include <vector>
#include <memory>
//
////Custom includes
#include "RunTimeParameters.h"
#include "EventsLoaderFactory.h"
#include "ReconstructorTemplate.h"
//#include "ImageAlgorithmFactory.h"
#include "PhantomVolumeBuilder.h"
#include "OctaneAlgorithmBuilder.h"
#include "OctaneEMAlgorithmBuilder.h"
#include "KEMAlgorithmBuilder.h"
#include "VectorAlgorithmBuilder.h"
#include "SBPAlgorithmBuilder.h"
#include "PointAlgorithmBuilder.h"
#include "SOEAlgorithmBuilder.h"
#include "EventDataAlgorithmBuilder.h"
#include "DCALineAlgorithmBuilder.h"
#include "Co60AlgorithmBuilder.h"
#include "ResultsFormatterFactory.h"
#include "StringUtils.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{        
    
class ReconstructorBuilder{ 
        
    public:
        ReconstructorBuilder(){ /* DO NOTHING*/ };
        ~ReconstructorBuilder(){};
        void  run();
        
        static shared_ptr<ReconstructorTemplate> build(const string &parameters_file_path){
            
            RunTimeParameters params(parameters_file_path);

            RandomSingleton::Instance()->setSeed(params.get_float("RANDOM_SEED"));
            
            auto reconstructor_ptr = std::make_shared<ReconstructorTemplate>();
            auto conic_sections = loadConicSections(params);
            
            auto IA = buildImageAlgorithm(conic_sections, params);
            
            reconstructor_ptr->setImageAlgorithm(IA);
            
            auto results_formatter = ResultsFormatterFactory::create(parameters_file_path, IA);
            reconstructor_ptr->setResultsFormatter(results_formatter);
            
            return reconstructor_ptr;            
        }
             
        
    private:
        static shared_ptr<ImageAlgorithm> buildImageAlgorithm(const vector<ConicSection> &conics, const RunTimeParameters &params){
            auto algo_str = params["IMAGE_ALGORITHM"];
            auto algorithm = pg_tools::StringUtils::to_upper(pg_tools::StringUtils::strip(algo_str));

            if (algorithm == "OCTANE") {
                cout << "Running octane algorithm . . ." << endl;
                OctaneAlgorithmBuilder b;
                auto algo = b.build(conics, params);
                return algo;
            }else if (algorithm == "OCTANE_EM") {
                    cout << "Running octane EM algorithm . . ." << endl;
                    OctaneEMAlgorithmBuilder b;
                    auto algo = b.build(conics, params);
                    return algo;
            }else if (algorithm == "VECTOR"){
                    cout<<"Running vector algorithm . . ."<<endl;
                    VectorAlgorithmBuilder b;
                    auto algo = b.build(conics, params);
                    return algo;
            }else if (algorithm == "SBP"){
                cout<<"Running SBP algorithm . . ."<<endl;
                SBPAlgorithmBuilder sbp;
                return sbp.build(conics, params);
            }else if (algorithm == "KEM"){
                cout<<"Running KEM algorithm . . ."<<endl;
                KEMAlgorithmBuilder a;
                return a.build(conics, params);
            }else if(algorithm == "POINT"){
                PointAlgorithmBuilder pab;
                return pab.build(conics, params);
            }else if(algorithm == "SOE"){
                SOEAlgorithmBuilder b;
                return b.build(conics, params);
            }else if(algorithm == "EVENT_DATA"){
                EventDataAlgorithmBuilder builder;
                return builder.build(conics, params);
            }else if(algorithm == "DCA_LINE"){
                DCALineAlgorithmBuilder builder;
                return builder.build(conics, params);
            }else if(algorithm == "CO60"){
                Co60AlgorithmBuilder builder;
                return builder.build(conics, params);
            }else{
                string err_msg("ERROR: Invalid IMAGE_ALGORITHM ");
                err_msg += params["IMAGE_ALGORITHM"];
                err_msg += ". Valid algorithms are SOE, POINT, and OCTANE.";
                throw std::runtime_error(err_msg);
            }
            
        }
        
        static vector<ConicSection> loadConicSections(const RunTimeParameters &params){
            /// \todo Most of this method should be in EventsLoaderFactory.

            auto data_file_format = params.get_int("DATA_FILE_FORMAT");
            auto phantom_volume = PhantomVolumeBuilder::build(params);
            
            auto eventsLoaderPtr = EventsLoaderFactory::create(data_file_format, params["EVENT_FILE_PATH"], &params, phantom_volume);
            vector<shared_ptr<ConicSection> > cs_ptrs;

//            auto percent_cone_phantom_overlap = params.get_double("CONE_PHANTOM_OVERLAP_PERCENTAGE");
//            auto number_tries_for_random = static_cast<int>(100 * (1.0/percent_con    e_phantom_overlap));
            auto number_tries_for_random = params.get_double("NUM_TRIES_FOR_RANDOM");
            eventsLoaderPtr->LoadEvents(cs_ptrs, number_tries_for_random);  
            
            
            /// \todo Get rid of this extra loop
            vector<ConicSection> cs;
            cs.reserve(cs_ptrs.size());
            for(size_t i = 0; i < cs_ptrs.size(); ++i){
                cs.push_back(*cs_ptrs[i]);;
            }    
            
            return cs;
        };

        //Image reconstruction steps 
        void load_compton_scatter_events_();
        void prepare_density_estimator_();
        void calculate_image_();
        void save_results_();


         std::time_t start_time_;
    };
};

#endif //RECONSTRUCTOR_BUILDER
