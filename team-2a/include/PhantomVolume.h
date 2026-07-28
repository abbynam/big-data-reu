#ifndef PHANTOM_VOLUME_H_
#define PHANTOM_VOLUME_H_
#define _USE_MATH_DEFINES

#include "RandomSingleton.h"
#include "PGVector3.h"

using namespace prompt_gamma_reconstruction;

namespace prompt_gamma_reconstruction{
    
 /*! \brief Defines the rectangular region for image reconstruction 
  * Phantom is a 3D rectangular tank. Represent the water tank
  * used in prompt gamma imaging studies.
  * 
  * @author Dennis Mackin
  */  
class PhantomVolume{
    public:
        double x_min;
        double x_max;
        size_t x_bins_;

        double y_min;
        double y_max;
        size_t y_bins_;

        double z_min;
        double z_max;
        size_t z_bins_;

        PhantomVolume(double x0, double x1, double y0, double y1, double z0, double z1):
                      x_min(x0), x_max(x1), y_min(y0), y_max(y1), z_min(z0), z_max(z1)
                      { /* DO NOTHING */};
                      
        PhantomVolume():x_min(0.0), x_max(0.0), y_min(0.0), y_max(0.0), z_min(0.0), z_max(0.0){};
        
        PhantomVolume(const PhantomVolume &lhs){
            x_min = lhs.x_min;
            x_max = lhs.x_max;
            y_min = lhs.y_min;
            y_max = lhs.y_max;
            z_min = lhs.z_min;
            z_max = lhs.z_max;
        };

        inline bool is_in_volume(const PGVector3 &point) const{
          return ( point.x < x_max && point.x > x_min
              && point.y < y_max && point.y > y_min
              && point.z < z_max && point.z > z_min);
        };

        PGVector3 get_random_point() const{
           PGVector3 p;
           RandomSingleton *_rand = prompt_gamma_reconstruction::RandomSingleton::Instance();

           p.x = _rand->getRand()*(x_max - x_min) + x_min;
           p.y = _rand->getRand()*(y_max - y_min) + y_min;
           p.z = _rand->getRand()*(z_max - z_min) + z_min;
           
           return p;
           
        };

        PGVector3 get_center_point() const{
            PGVector3 p(0.5*(x_min +x_max), 0.5*(y_min + y_max), 0.5*(z_min + z_max));
            return p;
        };

        float get_max_length() const{
            return max(x_max - x_min, max(y_max - y_min, z_max - z_min));
        }
        
        void print(){
          printf("\n---- phantom ----\n");
          printf("range x, %.3f, %.3f\n",x_min, x_max);
          printf("range y, %.3f, %.3f\n",y_min, y_max);
          printf("range z, %.3f, %.3f\n",z_min, z_max);
          printf("----------------\n");
        };
};

}
#endif //PHANTOM_VOLUME
