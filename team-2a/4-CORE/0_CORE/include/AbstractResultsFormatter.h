#ifndef _ABSTRACT_RESULTS_FORMATTER
#define _ABSTRACT_RESULTS_FORMATTER

//C++ standard includes
#include <memory>
#include <string>

// Custom includes
#include "ImageAlgorithm.h"

using namespace std;
namespace prompt_gamma_reconstruction{
    
    
class AbstractResultsFormatter{
    
public:
    AbstractResultsFormatter(const shared_ptr<const ImageAlgorithm> image_algo_ptr, const string &output_folder):
        image_algorithm_ptr_(image_algo_ptr), output_folder_(output_folder) { };

    virtual ~AbstractResultsFormatter() { };        
    virtual void saveOutput() const = 0;
              
protected:
    const shared_ptr<const ImageAlgorithm> image_algorithm_ptr_;
    string output_folder_;
    };
};
#endif // _ABSTRACT_RESULTS_FORMATTER