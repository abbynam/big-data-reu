#ifndef _RUNTIMEPARAMETERS
#define _RUNTIMEPARAMETERS

#include <iostream>
#include <fstream>
#include <string>
#include <vector> 
#include "StringMap.h"
#include "StringUtils.h"

using namespace std;

namespace pg_tools {

/*! \brief Extracts the parameter from the run configuration file.
 * 
 * The run configuration file uses a key value format for passing 
 * parameters to the application at runtime.
 *  
 * @author Dennis Mackin
 */
    class RunTimeParameters {
    public:

        RunTimeParameters(const string &runTimeParamsFilePath);

        //RunTimeParameters(const RunTimeParameters &params);
        RunTimeParameters(const RunTimeParameters *params);

        ~RunTimeParameters(void);

        size_t readParamFile(string runTimeParamsFilePath, StringMap &runTimeParams);

        string operator[](const string &key) const { return run_time_params_[key]; };

        string operator[](const string &key) { return run_time_params_[key]; };

        /// \todo Add error checking for these return types
        double get_double(const string &key) const {
            string s(run_time_params_[key]);
            double d = strtod(s.c_str(), 0);
            return d;
        };

        float get_float(const string &key) const {
            float f = strtof(run_time_params_[key].c_str(), 0);
            return f;
        };

        vector<float> get_csv_values(const string &key) const {

            auto vals_str = StringUtils::split(run_time_params_[key], ',');
            vector<float> vals(vals_str.size(), 0.0);

            for(size_t i = 0; i < vals_str.size(); ++i){
                auto val = StringUtils::strip(vals_str[i]);
                vals[i] = strtof(val.c_str(), 0);
            }
            return vals;
        };

        int get_int(const string &key) const { return atoi(run_time_params_[key].c_str()); };

        string getParametersfile() const { return parameters_file_; };


    private:

        StringMap run_time_params_;
        string parameters_file_;

        size_t readParamFile(string runTimeParamsFilePath);

    };
} //end of namespace p_tools

#endif //__RUNTIMEPARAMETERS



