#ifndef KERNEL_CONIC_H_
#define KERNEL_CONIC_H_

//standard C++/C libraries

//ROOT include files

//local include files
#include "ComptonScatter.h"
#include "ConicSection.h"
#include "RandomSqrtSingleton.h"
//#include "RandomPointOnCircleSingleton.h"


namespace prompt_gamma_reconstruction{

/*! \brief Class to model conic section for use with Kernel Density Esitimation.
 *  * 
 * @author Dennis Mackin
 */
  class KernelConic {

  private:

    float apex_x_;
    float apex_y_;
    float apex_z_;
    float cos_alpha_;
    float sin_alpha_;

    float tan_phi_;
    float cos_phi_;

    float cos_theta1_;
    float sin_theta1_;

    double calculateXZrotationAngle_(double x, double z) const;

  public:

    KernelConic(const ComptonScatter &comptonScatter, shared_ptr<const PhantomVolume> phantomVolume);
    double getDistanceToPoint(PGVector3 point) const;
    void transformPointToConeAxisFrameFromPhantomFrame(PGVector3 &point) const;

    void print();
  };


}//end of namespace prompt_gamma_reconstruction
#endif //KERNEL_CONIC_H_
