/* ****************************************************************************
 *  SimpleResultsFormatter -
 * 
 * @author Dennis Mackin
 * @date November 21, 2015
 */

// C++ Includes
#include <sstream>

// Custom Includes
#include "SimpleResultsFormatter.h"
#include "utilities/FileUtils.h"

using namespace std;
using namespace prompt_gamma_reconstruction;

string SimpleResultsFormatter::setup_output_folder_() const{
//////////////////////////////////////////////////////////////////////////////////////
/// Creates a new directory to store inputs and outputs for a looPDE run. This is useful
/// for keeping track of runs and their inputs and outputs.
/// - Creates the directory to store run information
/// - Copies parameters file to run folder
///
/// @param
/// @returns the path to the new output file
//////////////////////////////////////////////////////////////////////////////////////
    
    char command[501];
    //Create the folder if it doesn't exist 
    sprintf(command, "mkdir -p %s", output_folder_.c_str());

    if(system(command) != 0){
        printf("WARNING: failed to make %s . . .\n", output_folder_.c_str());
    }else{
        cout<<"Created directory "<< output_folder_<<". . . \n";
    }

    string output_file = output_folder_ + "/parameters.cfg";
    pg_tools::FileUtils::fileCopy(parameters_file_path_, output_file);

    return output_folder_;
};//end of SetUpOutputFolder