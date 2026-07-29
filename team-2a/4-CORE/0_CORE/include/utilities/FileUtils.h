#ifndef __FILEUTILS
#define __FILEUTILS

//Standard C++ Header Files
#include <vector>
#include <valarray>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
namespace pg_tools {

class FileUtils
{
public:
    FileUtils(void){};
    ~FileUtils(void){};

    static void ReadDataFromFile(vector< vector<float> > &data, const std::string &filePath, const std::string &delimiter);
    static size_t fileCopy(string sourcePath, string destinationPath, const string &copyCommand);
    static string fileCopy(const string &source_path, const string &output_folder);
};
}
#endif //__FILEUTILS
