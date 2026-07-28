#ifndef PGVECTOR3_H_
#define PGVECTOR3_H_
#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <iostream>
#include <ostream>
#include <sstream>


namespace prompt_gamma_reconstruction{

    /*! \brief 3-Vector used throughout the application
     * 
     * This reconstruction method is all about 3D transformations.
     * Having a custom 3D vector gives us the ability to make them
     * perform as we want them to perform.
     * 
     * @author Dennis Mackin
     */
    
    template<typename T> class PGVector3_T{
      public:
        T x;
        T y;
        T z;

        PGVector3_T<T>(T x=0.0, T y=0.0, T z=0.0 ): x(x), y(y), z(z){/* nothing else to do */};
        PGVector3_T<T>(const PGVector3_T<T>  &other): x(other.x), y(other.y), z(other.z){/* nothing else to do */};
        PGVector3_T<T>(PGVector3_T<T>  &other): x(other.x), y(other.y), z(other.z){/* nothing else to do */};

        inline PGVector3_T<T>  translate(PGVector3_T<T>  &in, PGVector3_T<T>  &translation_vector) {
            PGVector3_T<T> output;
            output += in;
            output += translation_vector;
            return output;
        }
        inline PGVector3_T<T>  translate(const PGVector3_T<T>  translation_vector){
            //return translate( *this, translation_vector);
            return *this + translation_vector;
        }

        inline void  set_xyz_array(T*  arr){
            arr[0] = x;
            arr[1] = y;
            arr[2] = z;
        }

        inline static PGVector3_T<T>  rotateXaxis(PGVector3_T<T>  &in, T angle) {
            PGVector3_T<T>  output;
            output.x = in.x;
            rotate(in.y, in.z, angle, output.y, output.z);
            return output;
        }

        inline static PGVector3_T<T>  rotateYaxis(PGVector3_T<T>  &in, T angle){
            PGVector3_T<T>  output;
            output.y = in.y;
            rotate(in.x, in.z, angle, output.x, output.z);
            return output;
        }

        inline void rotateYaxis(T angle){
            rotate(x, z, angle, x, z);
        }

        inline void rotateZaxis(T angle){
            rotate(x, y, angle, x, y);
        }

        inline static PGVector3_T<T>  rotateZaxis(PGVector3_T<T>  &in, T angle){
            PGVector3_T<T>  output;
            output.z = in.z;
            rotate(in.x, in.y, angle, output.x, output.y);
            return output;
        }

        inline static void rotate(T a_in, T b_in, T angle, T &a_out, T &b_out) {
            T cos_angle=0.0, sin_angle=0.0;
            sin_angle = sin(angle);
            cos_angle = cos(angle);
            a_out = a_in*cos_angle - b_in*sin_angle;
            b_out = a_in*sin_angle + b_in*cos_angle;
        }

        inline PGVector3_T<T>  operator+ (const PGVector3_T<T> & v) const
        {
            return PGVector3_T<T> ( x + v.x, y + v.y, z + v.z);
        }

        inline PGVector3_T<T>  operator- (const PGVector3_T<T> & v) const
        {
            return PGVector3_T<T>  ( x - v.x, y - v.y, z - v.z);
        }
        inline PGVector3_T<T> & operator+= (const PGVector3_T<T> & rhs)
        {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;

            return (*this);
        }

        inline PGVector3_T<T> & operator-= (const PGVector3_T<T> & rhs)
        {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;

            return (*this);
        }

        inline PGVector3_T<T> & operator*= (const PGVector3_T<T> & rhs)
        {
            x *= rhs.x;
            y *= rhs.y;
            z *= rhs.z;

            return (*this);
        }

        inline T getDistanceToPoint(const PGVector3_T<T> & rhs)
        {
            return sqrt((x-rhs.x)*(x-rhs.x)+(y-rhs.y)*(y-rhs.y)+(z-rhs.z)*(z-rhs.z));
        }      

        inline T magnitude() const
        {
            return sqrt(x*x+y*y+z*z);
        }

        inline T dotProduct(const PGVector3_T<T> & v) const{
            return x*v.x + y*v.y + z*v.z;
        }

        inline T dotProductNormalized(const PGVector3_T<T> & v) const
        {
            return  dotProduct(v)/(this->magnitude()*v.magnitude());
        }

        inline PGVector3_T<T>  normalize() const
        {
            return (*this) * (1.0f/magnitude());
        }


        inline bool operator== (const PGVector3_T<T> & v) {return (x == v.x && y == v.y && z == v.z);};

        PGVector3_T<T> operator* (const T scalar) const {return PGVector3_T<T> (scalar*x,scalar*y,scalar*z);} ;
        PGVector3_T<T> operator*= (const T scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;

            return *this;
        };

        PGVector3_T<T>  getReverse() const
        {
            PGVector3_T<T> output;
            output.x = -x;
            output.y = -y;
            output.z = -z;
            return output;
        }

        std::string print() const{
            std::stringstream ss (std::stringstream::in | std::stringstream::out);
            ss.precision(5);
            ss<<"("<<x << "," << y << "," << z << ")";
            std::cout<<ss.str()<<std::endl;
            return ss.str();
        }
    };

typedef PGVector3_T<float> PGVector3;
}


#endif //PGVECTOR3
