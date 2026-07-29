/////////////////////////////////////////////////////////////////
/// Utilities
/////////////////////////////////////////////////////////////////
/// Miscellaneous collection of functions for doing things like
/// copying stl vectors into arrays (the first function in this header)
///
/// Created 2011-12-12 by Dennis Mackin
/////////////////////////////////////////////////////////////////
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <cassert>
using namespace std;
template<typename T>void copy_vector_to_array(vector<T> stl_vector, T *array_ptr){
  for(size_t i=0; i < stl_vector.size(); ++i){
    array_ptr[i] = stl_vector[i];
  }
}
#endif //UTILITIES_H_
