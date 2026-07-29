///-----------------------------------------------------------------------
/// StringUtils.cpp
///-----------------------------------------------------------------------
/// @TODO This code should be replaced with equivalent from boost.
/// @author Dennis Mackin
//-----------------------------------------------------------------------


#include "StringUtils.h"

using namespace pg_tools;

string StringUtils::replace(const string& input, const string &oldText, const string &replacementText){
    vector<string> splitResults;
    string newString, delimiter, bookEnd("BE"); 

    //bookend makes this method work for case when replacement string is at the beginning or the end
    StringUtils::split(bookEnd + input + bookEnd, oldText, splitResults);

    for(size_t  i=0; i < splitResults.size(); i++){
        newString += delimiter;
        newString += splitResults[i];
        delimiter = replacementText;
    }

    //remove the bookend
    newString.erase(newString.size() - 2);
    newString.erase(0, 2);
    return newString;
}	

vector<string> StringUtils::split(const string& input, const char delimiter){
    vector<string> results;
    StringUtils::split(input, string(1,delimiter), results);
    return results;
}

size_t StringUtils::split(const string& input, const char delimiter, vector<string>& results){
    return StringUtils::split(input, string(1,delimiter), results);
}


size_t StringUtils::split(const string& original, const string& delimiter, vector<string>& results)
{
    size_t numFound = 0;
    int currentPos=0, findPos=0;
    results.clear();
    findPos = original.find(delimiter, 0);
    while( 0 <= findPos){
        assert( currentPos < (int)original.size() );
        if( findPos - currentPos >= 1){ //make sure the token value is not empty
            numFound++;
            results.push_back( original.substr(currentPos, findPos - currentPos) );
        }

        currentPos = findPos + delimiter.size();
        findPos = original.find( delimiter, currentPos);
    }

    if(original.size() - currentPos >= 1){ //add the remainder of the original to the output
        results.push_back( original.substr(currentPos) );
        numFound++;
    }
    return numFound;
}//end of split

string StringUtils::strip(string &line)
{
    size_t l = line.size();
    if ( l == 0 ) return string("");

    size_t n = 0;
    while ( n < l && (
            (line[n] == 0)    ||
            (line[n] == ' ' ) || 
            (line[n] == '"')  ||
            (line[n] == '\'') ||
            (line[n] == '\r') ||
            (line[n] == '\n') || 
            (line[n] == '\t'))) n++;

    size_t m = l - 1;
    while ( m > 0 && (
            (line[m] == 0)    ||
            (line[m] == ' ')  || 
            (line[m] == '"')  ||
            (line[m] == '\'') ||
            (line[m] == '\n') ||
            (line[m] == '\r') ||
            (line[m] == '\t')) ) m--;

    string stripped = line.substr(n,m-n+1);
    return stripped;
}

string StringUtils::strip(const char *line)
{	
    string strLine(line);
    return StringUtils::strip(strLine);
}
