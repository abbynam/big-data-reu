#ifndef _SIMPLE_RESULTS_FORMATTER
#define _SIMPLE_RESULTS_FORMATTER

//C++ standard includes
#include <memory>
#include <string>

// Custom includes
#include "ImageAlgorithm.h"
#include "AbstractResultsFormatter.h"
#include "RunTimeParameters.h"

using namespace std;
namespace prompt_gamma_reconstruction{
    
    
/*! A basic results formattor class for use on Linux Machines */    
class SimpleResultsFormatter: public AbstractResultsFormatter{
    
    public:
        SimpleResultsFormatter(const shared_ptr<const ImageAlgorithm> image_algo_ptr, const string &output_folder, const string &parameters_file):
            AbstractResultsFormatter(image_algo_ptr, output_folder), output_folder_(output_folder), parameters_file_path_(parameters_file){

            pg_tools::RunTimeParameters params(parameters_file);

            output_x_bins_ = params.get_int("OUTPUT_BINS_X");
            output_x_min_ = params.get_float("OUTPUT_X_MIN");
            output_x_max_ = params.get_int("OUTPUT_X_MAX");
            output_y_bins_ = params.get_int("OUTPUT_BINS_Y");
            output_y_min_ = params.get_int("OUTPUT_Y_MIN");
            output_y_max_ = params.get_int("OUTPUT_Y_MAX");
            output_z_bins_ = params.get_int("OUTPUT_BINS_Z");
            output_z_min_ = params.get_int("OUTPUT_Z_MIN");
            output_z_max_ = params.get_int("OUTPUT_Z_MAX");

            stringstream ss;
            if(output_x_min_ >= output_x_max_){
                ss <<"\n\nINVALID X OUTPUT RANGE: OUTPUT_X_MIN="<< output_x_min_ <<" , "<< "OUTPUT_X_MAX=" << output_x_max_ << " . . .\n"<<endl;
                throw runtime_error(ss.str());
            }
            if(output_z_min_ >= output_z_max_){
                ss <<"\n\nINVALID Y OUTPUT RANGE: OUTPUT_Y_MIN="<< output_y_min_ <<" , "<< "OUTPUT_Y_MAX=" << output_y_max_ << " . . .\n"<<endl;
                throw runtime_error(ss.str());
            }
            if(output_z_min_ >= output_z_max_){
                ss <<"\n\nINVALID Z OUTPUT RANGE: OUTPUT_Z_MIN="<< output_z_min_ <<" , "<< "OUTPUT_Z_MAX=" << output_z_max_ << " . . .\n"<<endl;
                throw runtime_error(ss.str());
            }

        };
            
        ~SimpleResultsFormatter() {cout<<"Destroying SimpleResultsFormatter  . . ."<<endl; };        
        void saveOutput() const {
            cout << "SimpleResultsFormatter: Saving the results . . ." << endl;
            setup_output_folder_();
            cout << "SimpleResultsFormatter: calling image_algorithm_ptr_->getDataAsString . . ." << endl;
            auto results = image_algorithm_ptr_->getDataAsString(output_x_bins_, output_x_min_, output_x_max_,
                                                                 output_y_bins_, output_y_min_, output_y_max_,
                                                                 output_z_bins_, output_z_min_, output_z_max_);
            write_results(results, "output.dat");
            auto sysmat_results = image_algorithm_ptr_->getSystemMatrixAsString(output_x_bins_, output_x_min_, output_x_max_,
                                                                 output_y_bins_, output_y_min_, output_y_max_,
                                                                 output_z_bins_, output_z_min_, output_z_max_);
            write_results(sysmat_results, "sysmat.dat");
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

        size_t output_x_bins_;
        float output_x_min_;
        float output_x_max_;
        size_t output_y_bins_;
        float output_y_min_;
        float output_y_max_;
        size_t output_z_bins_;
        float output_z_min_;
        float output_z_max_;
    };
};
#endif // _SIMPLE_RESULTS_FORMATTER