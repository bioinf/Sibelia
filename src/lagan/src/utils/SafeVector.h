// SafeVector.h
// ------------
// Class for array bounds checking.

// define ENABLE_CHECKS in order to enable array bounds checking.

#ifndef SAFEVECTOR_H
#define SAFEVECTOR_H

#include <assert.h>
#include <vector>

using namespace std;

// class derived from the STL std::vector
template<class TYPE>
class SafeVector : public std::vector<TYPE>{
public:

  // miscellaneous constructors
  SafeVector () {} 
  SafeVector (size_t size) : vector<TYPE>(size) {} 
  SafeVector (size_t size, const TYPE &value) : vector<TYPE>(size, value) {} 
  SafeVector (const SafeVector &source) : vector<TYPE>(source) {}

#ifdef ENABLE_CHECKS

  // [] array bounds checking
  TYPE &operator[](size_t index){
    assert (index >= 0 && index < size());
    return std::vector<TYPE>::operator[] (index);
  }

  // [] const array bounds checking
  const TYPE &operator[] (size_t index) const {
    assert (index >= 0 && index < size());
    return std::vector<TYPE>::operator[] (index) ;
  }

#endif

};

#endif
