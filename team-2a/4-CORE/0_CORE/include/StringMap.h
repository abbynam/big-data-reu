#ifndef __MAP
#define __MAP

#include<vector>
#include<string>
#include<map>
#include<iostream>
#include<cstdlib>

using namespace std;
namespace pg_tools {

/*! \brief Sting map that throws error if key is not found on lookup.
 * 
 * Encapsulate std::map<string, string> and adds a check to 
 * make sure the lookup key exists in the map.
 * 
 * @author Dennis Mackin
 */
class StringMap {
    
private:
  //typedef std::pair<string,string> Pair;
  std::map<string,string> pairs;
  StringMap(StringMap const &);      // private to prevent copying
  StringMap& operator=(const StringMap&); // private to prevent copying

  
public:
  StringMap();
  const string& operator[](const string&) const;
  //inline void insert(string key, string value){ this->vec.push_back(Pair(key,value));}
  inline void insert(string key, string value){ pairs[key] = value;};
  void print_all() const;
};
} //end of namespace pg_tools

#endif //__MAP