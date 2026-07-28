#include <cmath>
#include "utilities/FileUtils.h"
#include "StringUtils.h"

namespace pg_tools {

    /*! \brief Parses data files such as CSV
     * Reads data from file and puts it into a 2D array. Function is 
     * intended to with simple flat files such as CSV files.
     * 
     * @param data 2D array for holding the data
     * @param filepath The path of the file to be parsed
     * @param delimiter Simple delimiter such as a comma that separates numerical
     *              data values
     * 
     * @author Dennis Mackin
     */
    void FileUtils::ReadDataFromFile(vector< vector<float> > &data, const std::string &filePath, const std::string &delimiter){

        //make sure the data vector is empty
        data.clear();
        //Open the data file
        cout << "Reading file " << filePath.c_str() << " .\n";
        ifstream dataFile(filePath.c_str(), ios::in);

        if(! dataFile.is_open()){
            cout << "\n\nERROR: Could not open data file " << filePath << ".\n\n";
            cout << "Program aborting.\n\n";
            exit(EXIT_FAILURE);
        }else {
            cout << "Opened data file " << filePath << " for reading.\n";
        }

        //First, find out how many dimensions the data has
        string lineBuffer;
        vector<string> vecStrEvent;
//        getline(dataFile, lineBuffer);
//        StringUtils::split(lineBuffer, delimiter.c_str(), vecStrEvent);
//        auto d = vecStrEvent.size();
        vector<float> event8(8,0.0);
        vector<float> event12(12,0.0);
        vector<float> *event;
        size_t d = 0;

        while( 1 ){
            //Read the next line
            getline(dataFile, lineBuffer);
            if(dataFile.eof()) break;

            StringUtils::split(lineBuffer, delimiter.c_str(), vecStrEvent);
            d = vecStrEvent.size();
            if (d == 0){ //If 0 records are return, exit
                break;
            }else if(d == 8){
                event = &event8;
            }else if(12 == d){
                event = &event12;
            }else{
                cerr<<"WARNING BAD RECORD in "<< filePath << "\n" <<lineBuffer <<"\n\n";
                continue;
            }

            for(size_t  i=0; d > i; i++){
                (*event)[i] = atof( vecStrEvent[i].c_str() );
            }
            data.push_back(*event);


        }

        //Reset the file pointer to the beginning
        dataFile.close();

        cout<<"Read in "<< (int)data.size()<<" records from "<< filePath << " . . .\n";
    } //end of ReadDataFromFile

    ///fileCopy is makes a system call to copy a file
    size_t FileUtils::fileCopy(string sourcePath, string destinationPath, const string &copyCommand){
            string command(copyCommand);
            command += " ";
            command += sourcePath;
            command += " " + destinationPath;
            cout<<"calling system command {"<<command.c_str()<<"} . . ."<<endl;
            return system(command.c_str());
    };
    
    ///fileCopy is makes a system call to copy a file
    string FileUtils::fileCopy(const string &source_file_path, const string &output_file_path){

        std::ifstream  src(source_file_path, std::ios::binary);
        std::ofstream  dst(output_file_path, std::ios::binary);
        dst << src.rdbuf();
        
        return output_file_path;
    };
    

};
pg_tools::FileUtils futils;
