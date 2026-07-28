#ifndef IMAGE_ALGORITHM_FACTORY_H_
#define IMAGE_ALGORITHM_FACTORY_H_

//standard C++ includes
#include <memory>
#include <iostream>

//local package includes
#include "AbstractImageAlgorithmBuilder.h"
#include "OctaneAlgorithmBuilder.h"
//#include "SOEAlgorithmBuilder.h"

/*! \brief Creates an ImageAlgorithm object based on the <KBD>IMAGE_ALGORITHM</KBD> 
 * parameter in the parameters file. 
 * 
 * @author Dennis Mackin
 * @data Nov. 29, 2015
 */
namespace prompt_gamma_reconstruction{
    class ImageAlgorithmFactory{
        
    public:
        static shared_ptr<ImageAlgorithm> createImageAlgorithm(const vector<ConicSection> &conics, const pg_tools::RunTimeParameters &params){
            auto image_algo_builder_ptr = getImageAlgorithmBuilder(params);
            auto image_algorithm_ptr = image_algo_builder_ptr->build(conics, params);
            
            return image_algorithm_ptr;
        };
        
    private:
        static shared_ptr<AbstractImageAlgorithmBuilder> getImageAlgorithmBuilder(const RunTimeParameters &params){
            if(params["IMAGE_ALGORITH"] == "octane") {
                return std::make_shared<OctaneAlgorithmBuilder>();
            }else if (params["IMAGE_ALGORITH"] == "soe"){
                return std::make_shared<SOEAlgorithmBuilder>();
            } else{
                stringstream ss;
                ss<<"ERROR: Invalid image algorithm {"<<params["IMAGE_ALGORITH"]<<"} in the parameters file."<<endl;
                throw std::runtime_error(ss.str());
            }
        };
    };
}

#endif //IMAGE_ALGORITHM_FACTORY_H_