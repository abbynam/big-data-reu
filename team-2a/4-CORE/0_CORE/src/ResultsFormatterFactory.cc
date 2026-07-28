#include <sstream>
#include <cctype>
#include "RunTimeParameters.h"
#include "ResultsFormatterFactory.h"
#include "SimpleResultsFormatter.h"
#include "Co60ResultsFormatter.h"
#include "StringUtils.h"


using namespace prompt_gamma_reconstruction;

shared_ptr<AbstractResultsFormatter> ResultsFormatterFactory::create(const string &parameters_file_path, const shared_ptr<const ImageAlgorithm> image_algo_ptr)
{
    pg_tools::RunTimeParameters params(parameters_file_path);
    string formatter_type = pg_tools::StringUtils::to_upper(params["RESULTS_FORMAT"]);
    
    if(formatter_type == "SIMPLE") {
        return make_shared<SimpleResultsFormatter>(image_algo_ptr, params["OUTPUT_FOLDER_PATH"], parameters_file_path);
    }else if(formatter_type == "CO60"){
        return make_shared<Co60ResultsFormatter>(image_algo_ptr, params["OUTPUT_FOLDER_PATH"], parameters_file_path);
    }else{
        stringstream ss;
        ss<<"Invalid RESULTS_FORMAT: " << formatter_type <<". Valid types are [SIMPLE].\n";
        throw std::runtime_error(ss.str());
    };
    
    return nullptr;
};

