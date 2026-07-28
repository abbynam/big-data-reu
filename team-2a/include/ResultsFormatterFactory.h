#ifndef RESULTS_FORMATTER_FACTORY_H_
#define RESULTS_FORMATTER_FACTORY_H_

//standard C++ includes
#include <memory>

//local package includes
#include "AbstractResultsFormatter.h"

/*! \brief Creates the EventsLoader based on the type of file
 * used for input. They is specified by the <kbd>DATA_FILE_FORMAT</kbd> 
 * parameter.
 * 
 * @author Dennis Mackin
 */
namespace prompt_gamma_reconstruction{
    class ResultsFormatterFactory{
        
    public:
        static shared_ptr<AbstractResultsFormatter> create(const string &parameters_file_path, const shared_ptr<const ImageAlgorithm> image_algo_ptr);
    };
}

#endif //RESULTS_FORMATTER_FACTORY_H_