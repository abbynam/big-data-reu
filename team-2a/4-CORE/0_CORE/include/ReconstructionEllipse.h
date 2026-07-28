#ifndef RECONSTRUCTION_ELLIPSE_H_
#define RECONSTRUCTION_ELLIPSE_H_

//standard C++/C libraries
#include <cstdio>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
//#include "boost/shared_ptr.hpp"
#include <memory>


//local include files
#include "ComptonScatter.h"
#include "ConicSection.h"
#include "RandomSqrtSingleton.h"
#include "RandomSingleton.h"
//#include "RandomPointOnCircleSingleton.h"

using namespace pg_tools;
namespace prompt_gamma_reconstruction{


    struct AngularRange{
        float start_angle;
        float end_angle;
        float start_fraction;//starting fraction of whole
        float end_fraction; //ending fraction of whole
    };


   /*! \brief Represents the elliptical intersection of the cone
    *  and the Phantom.
    * 
    * The intersection of a cone and plane forms a 
    * 
    * @author Dennis Mackin
    */  
    class ReconstructionEllipse: public ConicSection {

    private:

		shared_ptr<Random> rand_;
        RandomSqrtSingleton *ptr_random_sqrt_;
        RandomSingleton *ptr_random_;
//        RandomPointOnCircleSingleton_T<float> *ptr_random_point_on_circle_;

    public:
        ReconstructionEllipse(const ComptonScatter & comptonScatter, shared_ptr<const PhantomVolume> phantomVolume, const size_t seed);

        bool doesConicIntersectPhantom();
        void calculateEllipseParameters(const float y, float &a, float &b, float &ellipse_x_min);

        size_t getRandomPointInPhantom0(PGVector3 &random_point);

        void getRandomPoint(const float y, const float ellipse_a, const float ellipse_b, const float ellipse_x_min,  const std::vector<AngularRange> & ranges, PGVector3 &random_point);
        void changeFramePhantomToConicSection();
        void changeFrameConicSectionToCone();
        void changeFrameConeToConicSection();
        void changeFrameConicSectionToPhantom();
        void getLineIntercepts_(float a, float b, float h, float y,
                                std::pair<float, float> point1, std::pair<float, float> point2,
                                std::vector<PGVector3> &intersection_points);
        void getPhantomEllipseIntercepts(float a, float b, float x_offset, float y, std::vector<PGVector3> &intersection_points);
        float getAngle(const float x, const float a, const float z, const float b);
        void convertPointsToAngles(const std::vector<PGVector3> &points, float ellipse_a, float ellipse_x_min, std::vector<float> &angles);
        void getAngularRanges(float y, std::vector<float> &angles, std::vector<AngularRange> &ranges);

        //void draw(float a, float b, float h, Rectangle rectanglePoints, size_t numberIntersections, PGVector3 &random_point);
        //void draw(float y);
        void print();

    };


}//end of namespace prompt_gamma_reconstruction
#endif //RECONSTRUCTION_ELLIPSE_H_
