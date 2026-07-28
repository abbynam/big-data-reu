#include <sstream>
#include <stdexcept>
#include "RunTimeParameters.h"

using namespace std;

namespace pg_tools {
    RunTimeParameters::RunTimeParameters(const string &parameters_file): parameters_file_(parameters_file){
        readParamFile(parameters_file);
    };
    
    RunTimeParameters::RunTimeParameters(const RunTimeParameters *params){
        parameters_file_ = params->getParametersfile();
        readParamFile(parameters_file_);
    };    

    RunTimeParameters::~RunTimeParameters(void){};

    /** reads in the configuration file and writes the key value parairs to runTimeParams
            @params configFilePath The path to the configuration file
            @params an address of a map to store to the key-value pairs
            @returns The number of run time parameters set.
    */
    size_t RunTimeParameters::readParamFile(string runTimeParamsFilePath){

        ifstream paramsFile;

        cout << "Reading parameters file " << runTimeParamsFilePath.c_str() << ".\n";
        paramsFile.open(runTimeParamsFilePath.c_str(), ios::in);

        if(! paramsFile.is_open()){
            stringstream ss;
            ss << "ERROR: Could not open parameters file " << runTimeParamsFilePath 
               << " in RunTimeParameters::readParamFile.";
            cout<<ss.str()<<endl;
            throw runtime_error(ss.str());
        }

        string lineBuffer;
        vector<string> keyValuePair;
        cout<<"Reading in the parameters . . ."<<endl;
        while (! paramsFile.eof() )
        {
            getline(paramsFile, lineBuffer);
            lineBuffer = StringUtils::strip(lineBuffer); //Clean off the whitespace on the ends
            if(lineBuffer[0] == '#' || lineBuffer.size() == 0){//skip comment lines
                continue;
            }

            StringUtils::split( lineBuffer, '=', keyValuePair);

            if(keyValuePair.size() != 2){
                    cout<<"keyValuePair.size() == "<<keyValuePair.size()<<endl;
                    cout<<"\nERROR: BAD PARAMETER FILE RECORD\n-------------------------------\n"
                            << lineBuffer << endl <<endl;

                            continue; //skip this record.
            }
            run_time_params_.insert( StringUtils::strip(keyValuePair[0]), StringUtils::strip( keyValuePair[1] ) );
            keyValuePair.clear();
        }
        return 1; //Return True
    }
} //end of namespace pde_tools
