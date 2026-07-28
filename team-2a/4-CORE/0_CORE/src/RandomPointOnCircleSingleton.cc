
#include "RandomPointOnCircleSingleton.h"
using namespace prompt_gamma_reconstruction;

template<typename T> RandomPointOnCircleSingleton_T<T> *RandomPointOnCircleSingleton_T<T>::ptrToSelf = NULL;
template<> RandomPointOnCircleSingleton_T<float> *RandomPointOnCircleSingleton_T<float>::ptrToSelf = NULL;

