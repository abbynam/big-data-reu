#include "KernelConic.h"
#include <cmath>

using namespace std;
using namespace prompt_gamma_reconstruction;


KernelConic::KernelConic(const ComptonScatter & comptonScatter, shared_ptr<const PhantomVolume> phantomVolume){

  apex_x_= comptonScatter.getConeApex().x;
  apex_y_= comptonScatter.getConeApex().y;
  apex_z_= comptonScatter.getConeApex().z;
  
  double alpha = comptonScatter.getAlpha();
  cos_alpha_ = cos(alpha);
  sin_alpha_ = sin(alpha);
  float phi_ = comptonScatter.getScatteringAngle()*M_PI/180.0;
  //tan_alpha_plus_phi_ = tan(alpha_ + phi_);
  //tan_alpha_minus_phi_ = tan(alpha_ - phi_);
  tan_phi_ = tan(phi_);
  cos_phi_ = cos(phi_);

  PGVector3 tmpVec = comptonScatter.getConeAxis();
  float theta1_ = calculateXZrotationAngle_(tmpVec.x, tmpVec.z);
  cos_theta1_ = cos(theta1_);
  sin_theta1_ = sin(theta1_);

};

void KernelConic::transformPointToConeAxisFrameFromPhantomFrame(PGVector3 &point) const{
  //translate to cone axis
  //point -= compton_scatter_.getConeApex();
  point.x -= apex_x_;
  point.y -= apex_y_;
  point.z -= apex_z_;
  
  float tmp;
  tmp = point.x*cos_theta1_ - point.z*sin_theta1_;
  point.z = point.x*sin_theta1_ + point.z*cos_theta1_;
  point.x = tmp;

  tmp = point.x*cos_alpha_ + point.y*sin_alpha_;
  point.y = - point.x*sin_alpha_ + point.y*cos_alpha_;
  point.x = tmp;

  point.y *= -1.0;
}

///////////////////////////////////////////////////////////////////////////////
double KernelConic::getDistanceToPoint(PGVector3 point) const{
///////////////////////////////////////////////////////////////////////////////
/// @param point in phantom frame for which the distance is to be calculated.
///
/// @returns the smallest distance between the point and the cone
///////////////////////////////////////////////////////////////////////////////

  transformPointToConeAxisFrameFromPhantomFrame(point);
  double distance_from_y_axis = sqrt(point.x*point.x + point.z*point.z);
  double radius = point.y*tan_phi_;
  
  double distance = abs(distance_from_y_axis - radius)*cos_phi_;

  return distance;
}

double KernelConic::calculateXZrotationAngle_(double x, double z)const {

  double rotation_angle = 0.0;
  if(0.0 == x){
    if(0.0 == z){
      rotation_angle = 0.0;
    }else{
      rotation_angle = M_PI/2.0;
    }
  }else{
    rotation_angle = fabs(atan(z/x));

    //now need to correct based on the quadrants
    if(0.0 > x){
      rotation_angle = M_PI - rotation_angle;
    }
    if(0.0 < z){ 
      rotation_angle *= -1.0;
    }
  }
  return rotation_angle;
}

void KernelConic::print(){
  
  printf("-------KernelConic---------\n");
  printf("cone apex(%.2f, %.2f, %.2f)\n", apex_x_, apex_y_, apex_z_);

  printf("--------------------------------------\n");
  std::cout<<std::endl;

};


