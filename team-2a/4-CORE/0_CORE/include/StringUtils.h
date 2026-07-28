//-----------------------------------------------------------------------
// StringUtils.h
//-----------------------------------------------------------------------

#ifndef __STRINGUTILS_H_
#define __STRINGUTILS_H_

#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream> 
#include <cassert>
#include <algorithm>

using namespace std;

namespace pg_tools{
/**Class contains sting methods common to high level languages like Python and Perl.
    @author Dennis Mackin 
    @date Aug. 25, 2004
*/
class StringUtils
{

  public:
/** @param input string to be split
    @param delimiter characters delimiting the breaks in the string
    @param results string vector used to store the results
    @returns number of substrings created */
        static size_t split(const string& input, const string& delimiter, vector<string>& results);
        static size_t split(const string& input, const char delimiter, vector<string>& results);
        static vector<string> split(const string& input, const char delimiter);
        static string replace(const string& input, const string &oldText, const string &replacementText);
        static string strip(string &); ///Removes leading and trailing whitespace characters
        static string strip(const char *); ///Removes leading and trailing whitespace characters

    static size_t strtoi(string str){ return atoi(str.c_str()); };

    ///convert an integer to a string
    inline string itostring(const size_t i){char buff[10]; sprintf(buff,"%zu",i); return string(buff);}

    ///convert an double to a string
    inline string dtostring(const double d){char buff[30]; sprintf(buff,"%.8f",d); return string(buff);}
    inline static string to_upper(const string rhs){ string s = rhs; transform(s.begin(), s.end(), s.begin(), ::toupper); return s; };
    inline static string to_lower(const string rhs){ string s = rhs; transform(s.begin(), s.end(), s.begin(), ::tolower); return s; };
};

}//end of pg_tools

#endif //end of __STRINGUTILS_H_


