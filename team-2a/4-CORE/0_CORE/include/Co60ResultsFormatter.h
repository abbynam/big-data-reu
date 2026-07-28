#ifndef _CO60_RESULTS_FORMATTER
#define _CO60_RESULTS_FORMATTER

//C++ standard includes
#include <memory>
#include <string>

// Custom includes
#include "ImageAlgorithm.h"
#include "AbstractResultsFormatter.h"
#include "StringUtils.h"

using namespace std;
using namespace pg_tools;
namespace prompt_gamma_reconstruction{
    
    
/*! A basic results formattor class for use on Linux Machines */    
class Co60ResultsFormatter: public AbstractResultsFormatter{
    
    public:
        Co60ResultsFormatter(const shared_ptr<const ImageAlgorithm> image_algo_ptr, const string &output_folder, const string &parameters_file):
            AbstractResultsFormatter(image_algo_ptr, output_folder), output_folder_(output_folder), parameters_file_path_(parameters_file){ 
        
                cout<<"Created Co60ResultsFormatter  . . ."<<endl;
        };
            
        ~Co60ResultsFormatter() {cout<<"Destroying Co60ResultsFormatter  . . ."<<endl; };
        void saveOutput() const {
            setup_output_folder_();

            auto results = image_algorithm_ptr_->getDataAsString();
            cout<<"results size: "<<results.size()<<endl;

            auto parts = StringUtils::split(results, '~');
            assert(2 == parts.size());
            write_results(StringUtils::strip(parts[0]), "output.dat");
            write_results(StringUtils::strip(parts[1]), "output_corrected.dat");

            write_results(image_algorithm_ptr_->getConicInformationAsString(), "events.dat");
        };
        
    private:
        string setup_output_folder_() const;//end of SetUpOutputFolder
        void write_results(const string &results, const string &filename) const{
            string ofile = output_folder_ + "/" + filename;
            std::ofstream out(ofile);
            out<< results;
        }
        string output_folder_;
        string parameters_file_path_;
    };
};
#endif // _CO60_RESULTS_FORMATTER