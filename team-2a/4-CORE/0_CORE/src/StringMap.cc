#include "StringMap.h"
#include <sstream>
#include <stdexcept>

namespace pg_tools {
    
    /*! \brief Return value for key. If not found throw error.
      */     

    
    StringMap::StringMap(const StringMap&){
        
        cout<<"Calling StringMap private copy constructor ????\n"<<endl;
    }
    
    StringMap::StringMap(){
        
        cout<<"Calling StringMap default constructor\n"<<endl;
    }    
    const string& StringMap::operator[](const string& key) const {
        
        std::map<string,string>::const_iterator value = pairs.find(key);
        
        if (value != pairs.end()) return value->second;
        
        //Key not found. Throw error.
        stringstream ss;
        ss <<"\n\nKEY_LOOKUP_ERROR: key "<< key <<" not found. "<< "Exiting . . .\n"<<endl;
        throw runtime_error(ss.str());
    }

    /*! \brief Print the key/value pairs
      */    
    void StringMap::print_all() const
    {
        std::map<string,string>::const_iterator iter = pairs.begin();
        for( /* */; iter!=pairs.end(); ++iter){
            cout << iter->first << ": " << iter->second << '\n';
        }
    }
} //end of namespace pg_tools
