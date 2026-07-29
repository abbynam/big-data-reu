#ifndef RANDOM_Y_SINGLETON_H_
#define RANDOM_Y_SINGLETON_H_
#define _USE_MATH_DEFINES

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <ostream>
#include <sstream>
#include <vector>
#include <algorithm>

//LIBRARY INCLUDES

//Prompt Gamma Includes
#include "ComptonScatter.h"
#include "utilities/Random.h"

namespace prompt_gamma_reconstruction{

/*! \brief Returns a precalculated point on a circle
 * 
 * Produces random x,y pairs on the unit circle. The random points are 
 * generated on initialization. Users should pay attention to how this
 * class is used because the period is fairly short.
 * 
 * @author Dennis Mackin
 */
 
    
template<typename T> class RandomPointOnCircleSingleton_T { 

  public:
    inline std::pair<T,T> & getPoint(){
        if(NUMBER_OF_POINTS == current_point_) current_point_ = 0;
        return points[current_point_++];
    };

    inline void getPoint(T &x, T &y){
        x = points[current_point_].first;
        y = points[current_point_].second;
        if(NUMBER_OF_POINTS == ++current_point_) current_point_ = 0;
    };

    static RandomPointOnCircleSingleton_T<T> * Instance(){
        if( NULL == ptrToSelf){
            ptrToSelf = new RandomPointOnCircleSingleton_T<T>;
        }
        return ptrToSelf;
    };

  private:
    pg_tools::Random rand_;
    static const int NUMBER_OF_POINTS = 10000079;
    int current_point_;

    std::vector<std::pair<T,T> > points;

    static RandomPointOnCircleSingleton_T<T> *ptrToSelf;

    RandomPointOnCircleSingleton_T<T>():rand_(0),current_point_(0){ 
        for(int i=0; i<NUMBER_OF_POINTS; ++i){
            std::pair<T,T> point;
            T angle = 2.0*M_PI*rand_.Rndm();
            point.first = cos(angle);
            point.second = sin(angle);
            points.push_back(point);
        };
        std::random_shuffle(points.begin(),points.end());
    };

};//end of class

typedef RandomPointOnCircleSingleton_T<float> RandomPointOnCircleSingleton;
};//end of namespace
#endif //RANDOM_Y_SINGLETON_H_
