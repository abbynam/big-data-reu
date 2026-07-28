#include "ReconstructionEllipse.h"
#include "utilities/Random.h"

using namespace std;
using namespace prompt_gamma_reconstruction;



ReconstructionEllipse::ReconstructionEllipse(const ComptonScatter & comptonScatter, shared_ptr<const PhantomVolume> phantomVolume, const size_t seed):
                      ConicSection(comptonScatter, phantomVolume, seed){

  ptr_random_sqrt_ = RandomSqrtSingleton::Instance();
  ptr_random_ = RandomSingleton::Instance();
//  ptr_random_point_on_circle_ = RandomPointOnCircleSingleton::Instance();

  rand_ = shared_ptr<pg_tools::Random>(new pg_tools::Random());

};


void ReconstructionEllipse::getLineIntercepts_(float a, float b, float vertex_x, float y,
                                               std::pair<float, float> point1, std::pair<float, float> point2,
                                               vector<PGVector3> &intersection_x_values){
  /// Solutions take the following form after shifting x coord of center to origin ///
  /******************************************************************************

x -> -((a^2 intercept slope + Sqrt[a^2 b^2 (b^2 - intercept^2 + a^2 slope^2)])/(b^2 + a^2 slope^2))

x -> (-a^2 intercept slope + Sqrt[a^2 b^2 (b^2 - intercept^2 + a^2 slope^2)])/(b^2 + a^2 slope^2)
  ******************************************************************************/
  /// If the ellipse and the line intersect then the radical term will be real.
  /// There is no intersection if the term under the radical is negative.

  //change x coord for center of ellipse to 0
  point1.first -= a + vertex_x;
  point2.first -= a + vertex_x;

  float slope_denominator = point1.first - point2.first;
  float slope_numerator = point1.second - point2.second;
  PGVector3 point_neg, point_pos;
  //float x_neg, x_pos, z_neg, z_pos;

  if(0.0 == slope_denominator){ //then x is constant
    if( fabs(point1.first) > a ){//then the line doesn't intercept ellipse
      return;
    }
    point_pos.z = sqrt(1 - (point1.first)*(point1.first)/(a*a))*b;
    point_neg.z = -point_pos.z;
    point_neg.x = point1.first;
    point_pos.x = point1.first;
  }else if(0.0 == slope_numerator){ //then z is constant
    if( fabs(point1.second) > b ){//then the line doesn't intercept ellipse
      return;
    }
    float radical_term = sqrt(1 - point1.second*point1.second/(b*b))*a;
    point_pos.x =  radical_term;
    point_neg.x = -point_pos.x;
    point_neg.z = point1.second;
    point_pos.z = point1.second;
  }else{
    float slope = slope_numerator/slope_denominator;
    float intercept = point1.second - slope*point1.first;
    //float radical_term = a*a*b*b*(b*b -(intercept-a*slope+vertex_x*slope)*(intercept+(a+vertex_x)*slope));
    float radical_term = a*a*b*b*(b*b - intercept*intercept + a*a*slope*slope);
    if (0.0 > radical_term){//no real intercepts
      return;
    }
    float term1 = -a*a*intercept*slope;
    float denominator = b*b + a*a*slope*slope;

    point_neg.x = (term1 - sqrt(radical_term))/denominator;
    point_pos.x = (term1 + sqrt(radical_term))/denominator;
    point_pos.z = slope*point_pos.x + intercept;
    point_neg.z = slope*point_neg.x + intercept;
  }

  //Test to see that the
  if(point_neg.x >= std::min(point1.first, point2.first) && point_neg.x <= std::max(point1.first, point2.first)
    && point_neg.z >= std::min(point1.second, point2.second) && point_neg.z <= std::max(point1.second, point2.second)){
    point_neg.x += a + vertex_x;//shift x coordinate back to original origin
    intersection_x_values.push_back(point_neg);//add it to the intersections vector
  }
  if(point_pos.x >= std::min(point1.first, point2.first) && point_pos.x <= std::max(point1.first, point2.first)
    && point_pos.z >= std::min(point1.second, point2.second) && point_pos.z <= std::max(point1.second, point2.second)){
    point_pos.x += a + vertex_x;//shift x coordinate back to original origin
    intersection_x_values.push_back(point_pos);//add it to the intersections vector
  }
};

//
//void ReconstructionEllipse::draw(float y){
//  float a, b, x_min;
//  PGVector3 point;
//  calculateEllipseParameters(y, a, b, x_min);
//  std::vector<PGVector3> intersections;
//  point.y = y;
//  getPhantomEllipseIntercepts(a, b, x_min, y, intersections);
//  draw(a, b, x_min, getRotatedPhantomBoundary(), intersections.size(), point);
//
//  return;
//}




///////////////////////////////////////////////////////////////////////////////
///doesConicIntersectPhantom()
///
/// returns true if cone passes through phantom
///
/// Order of checks:
///  1. First check to see if the inside edge of the conic touches the phantom surface.
///  2. Check for any intersections at top surface
///  3. Check for any intersections at bottom surface
///////////////////////////////////////////////////////////////////////////////
bool ReconstructionEllipse::doesConicIntersectPhantom(){

  float epsilon = 1E-5; //Small offset so that we can test inside phantom, not on edge.
  PGVector3 intercept =  getConeInnerYPlaneIntercept(phantom_volume_->y_max - epsilon);
  if( phantom_volume_->is_in_volume(intercept)){
    return true;
  }
  if( phantom_volume_->is_in_volume( getConeInnerYPlaneIntercept(phantom_volume_->y_min + epsilon))){
    return true;
  }
  if( phantom_volume_->is_in_volume( getConeOuterYPlaneIntercept(phantom_volume_->y_max - epsilon))){
    return true;
  }
  if( phantom_volume_->is_in_volume( getConeOuterYPlaneIntercept(phantom_volume_->y_min + epsilon))){
    return true;
  }

  //check to see if ellipse intersects with top surface
  float ellipse_a, ellipse_b, ellipse_x_min;
  vector<PGVector3> intersection_points;

  const size_t NUMBER_Y_LAYERS = 5;
  for(size_t i=0; i<=NUMBER_Y_LAYERS; ++i){
    //check to see if ellipse intersects with bottom surface
    float y = phantom_volume_->y_min + (phantom_volume_->y_max - phantom_volume_->y_min)*((float)i/NUMBER_Y_LAYERS);
    calculateEllipseParameters(y, ellipse_a, ellipse_b, ellipse_x_min);
    getPhantomEllipseIntercepts( ellipse_a, ellipse_b, ellipse_x_min, phantom_volume_->y_min, intersection_points);
    if( intersection_points.size() > 0) return true;
  }


  //check to see if ellipse intersects at level where inner edge x = x_min
  PGVector3 axis(compton_scatter_.getConeAxis());
  PGVector3 apex(compton_scatter_.getConeApex());
  float axis_ratio = (axis.x != 0.0)? axis.y/axis.x : 0.0;
  float y = apex.y + (apex.x - phantom_volume_->x_min + epsilon)*axis_ratio /tan(alpha_ - phi_);
  if(phantom_volume_->y_max - epsilon > y && phantom_volume_->y_min + epsilon< y ){
    calculateEllipseParameters(y, ellipse_a, ellipse_b, ellipse_x_min);
    getPhantomEllipseIntercepts( ellipse_a, ellipse_b, ellipse_x_min, y, intersection_points);
    if( intersection_points.size() > 0) return true;
  }

  //check to see if ellipse intersects at level where inner edge x = x_max
  y = apex.y + (apex.x - phantom_volume_->x_max - epsilon)*axis_ratio /tan(alpha_ - phi_);
  if(phantom_volume_->y_max - epsilon > y && phantom_volume_->y_min + epsilon < y ){
    calculateEllipseParameters(y, ellipse_a, ellipse_b, ellipse_x_min);
    getPhantomEllipseIntercepts( ellipse_a, ellipse_b, ellipse_x_min, y, intersection_points);
  }
  if( intersection_points.size() > 0) return true;

  //no cone/phantom intersection
  return false;
}




///////////////////////////////////////////////////////////////////////////////
void ReconstructionEllipse::calculateEllipseParameters(const float y,
                                                       float &a,
                                                       float &b,
                                                       float &ellipse_x_min){
///////////////////////////////////////////////////////////////////////////////
/// @param y is the phantom y coordinate
/// @param a is half the length of the major axis
/// @param b is half the length of the minor axis
/// @param ellipse_x_min is minimal x value in the reference frame with
/// the x axis is colinear with the ellipse major axis
///////////////////////////////////////////////////////////////////////////////
    float delta_y = getComptonScatter().getConeApex().y - y;
//     assert(delta_y > 0.0);

    //float ellipse_x_max = delta_y*tan_alpha_plus_phi_;
    //a = (ellipse_x_max - ellipse_x_min)/2.0; //half length of major axis
    //ellipse_x_min = delta_y*tan_alpha_minus_phi_;
    //float ellipse_x_max = delta_y*tan_alpha_plus_phi_;
    PGVector3 inner_intercept = getConeInnerYPlaneIntercept(y);
    PGVector3 inner_intercept_ellipse_frame = inner_intercept;
    inner_intercept_ellipse_frame -= getComptonScatter().getConeApex();
    inner_intercept_ellipse_frame = PGVector3::rotateYaxis(inner_intercept_ellipse_frame, theta1_);

    PGVector3 outer_intercept = getConeOuterYPlaneIntercept(y);
    PGVector3 outer_intercept_ellipse_frame = outer_intercept;
    outer_intercept_ellipse_frame -= getComptonScatter().getConeApex();
    outer_intercept_ellipse_frame = PGVector3::rotateYaxis(outer_intercept_ellipse_frame, theta1_);
    ellipse_x_min = inner_intercept_ellipse_frame.x;
    a = abs((outer_intercept_ellipse_frame.x - inner_intercept_ellipse_frame.x)/2.0); //half length of major axis
    b = delta_y*tan_phi_;

    return;
}
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
void ReconstructionEllipse::getAngularRanges(float y, std::vector<float> &angles, std::vector<AngularRange> &ranges){
///////////////////////////////////////////////////////////////////////////////
/// @param angles is a vector of angles for phantom/ellipse intercepts
/// @param ranges is the output of angular
/// @param b is half the length of the minor axis
/// @param ellipse_x_min is minimal x value in the reference frame with
/// the x axis is colinear with the ellipse major axis
///////////////////////////////////////////////////////////////////////////////

  //if the outer most point on the ellipse is in the phantom then the first
  // intersection is an exit from the phantom. Thus the range from the
  // maximum to the first intersetion is inside the phantom.
  AngularRange range;
  float total_range = 0.0;
  if( this->phantom_volume_->is_in_volume(getConeOuterYPlaneIntercept(y))){
    //range.opening_angle = *itr;
    angles.insert(angles.begin(), 0.0);
    angles.push_back(2.0*M_PI);
    //range.start_angle = 0.0;
    //range.end_angle = *itr;
    //ranges.push_back(range);
    //total_range += *itr;
    //++itr;
  }

  std::vector<float>::iterator itr = angles.begin();
  for(;itr != angles.end(); ++itr){
    float opening_angle = *(itr + 1) - *itr;
    total_range += opening_angle;
    range.start_angle = *itr;
    itr++;
    range.end_angle = *itr;
    ranges.push_back(range);
  }
  float total_fraction = 0.0;
  std::vector<AngularRange>::iterator itr_ranges = ranges.begin();
  for(;itr_ranges != ranges.end(); ++itr_ranges){
    itr_ranges->start_fraction = total_fraction;
    total_fraction += (itr_ranges->end_angle - itr_ranges->start_angle)/total_range;
    itr_ranges->end_fraction = total_fraction;
  }

  if(total_fraction < 0.9999 || total_fraction > 1.00001){
    print();
    draw(y);
    cout<<"Total Fraction is only: "<<total_fraction<<".\n";
    stringstream ss;
    ss << "Error in ReconstructionEllipse::getAngularRanges\n"
       << "Sum of range fractions = " << total_fraction
       << ".\n Should be 1.0."<<endl;
    string error_message(ss.str());
    throw runtime_error(error_message) ;
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
void ReconstructionEllipse::getRandomPoint(const float y, const float ellipse_a, const float ellipse_b, const float ellipse_x_min,
                                           const std::vector<AngularRange> & ranges, PGVector3 &random_point){
///////////////////////////////////////////////////////////////////////////////

  float random_fraction = this->rand_->Rndm();
//   float random_fraction = (float)RandomSingleton::Instance()->getRand();
  float random_angle = 0.0;
  random_point.y = y;

  std::vector<AngularRange>::const_iterator itr =ranges.begin();
  for( ; itr != ranges.end(); ++itr){
    if( random_fraction < itr->end_fraction){
      random_angle = (itr->end_angle - itr->start_angle)*(random_fraction - itr->start_fraction) + itr->start_angle;
      random_point.x = ellipse_a * cos(random_angle) + ellipse_x_min + ellipse_a;
      random_point.z = ellipse_b * sin(random_angle);
      transformPointFromXZprimeToPhantom(random_point);
      if(! phantom_volume_->is_in_volume(random_point)){
        cout<<"WARNING: random point is not involume."<<endl;
        print();
        draw(random_point.y);
        string error_message("ERROR: ReconstructionEllipse::getRandomPoint\n");
        error_message += "failed to generate random point in the phantom.";
        cout<<error_message<<endl;
        throw runtime_error(error_message) ;

      }

      return;
    }
  }

}
///////////////////////////////////////////////////////////////////////////////
void ReconstructionEllipse::getPhantomEllipseIntercepts(float a, float b, float x_offset,
                                                        float y, vector<PGVector3> &intersection_points) {
    intersection_points.clear();

    getLineIntercepts_(a, b, x_offset, y,
                        rotated_phantom_boundary_.x_min_z_min, rotated_phantom_boundary_.x_min_z_max,
                        intersection_points);
    getLineIntercepts_(a, b, x_offset, y,
                        rotated_phantom_boundary_.x_min_z_max, rotated_phantom_boundary_.x_max_z_max,
                        intersection_points);
    getLineIntercepts_(a, b, x_offset, y,
                        rotated_phantom_boundary_.x_max_z_max, rotated_phantom_boundary_.x_max_z_min,
                        intersection_points);
    getLineIntercepts_(a, b, x_offset, y,
                        rotated_phantom_boundary_.x_max_z_min, rotated_phantom_boundary_.x_min_z_min,
                        intersection_points);
    return;
}


//
//int ReconstructionEllipse::getRandomPointInPhantom0(PGVector3 &random_point){
//
//  const int max_num_tries = 5000;
//  float ellipse_a=0.0, ellipse_b=0.0, ellipse_x_min=0.0;
//  for(int number_tries=0; number_tries < max_num_tries; ++number_tries){
//    random_point.y = rand_->Uniform( phantom_volume_->y_min, phantom_volume_->y_max);
//
//    calculateEllipseParameters(random_point.y, ellipse_a, ellipse_b, ellipse_x_min);
//
//    vector<PGVector3> intersections;
//    getPhantomEllipseIntercepts(ellipse_a, ellipse_b, ellipse_x_min, random_point.y, intersections);
//    PGVector3 outer_intercept = getConeOuterYPlaneIntercept(random_point.y);
//    if(intersections.size()==0 && ! phantom_volume_->is_in_volume(outer_intercept)){
//      continue;
//    }
//
//    vector<float> angles;
//    convertPointsToAngles(intersections, ellipse_a, ellipse_x_min, angles );
//    std::vector<AngularRange> ranges;
//    getAngularRanges(random_point.y, angles, ranges);
//    getRandomPoint(random_point.y, ellipse_a, ellipse_b, ellipse_x_min, ranges, random_point);
//
//    if(phantom_volume_->is_in_volume(random_point)) {
//      if(false == isLikelyOriginSet_) setLikelyOrigin(random_point);
//      return number_tries;
//    }
//  }
//  //return max_num_tries;
//  print();
//  draw(random_point.y);
//  string error_message("ERROR: ReconstructionEllipse::getRandomPointInPhantom\n");
//  error_message += "failed to generate random point in the phantom.";
//  cout<<error_message<<endl;
//  throw runtime_error(error_message) ;
//}



///////////////////////////////////////////////////////////////////////////////
/// getAngle
///////////////////////////////////////////////////////////////////////////////
/// calculates the angle for a point on an ellipse specified by x,z,a,b.
/// IMPORTATNT-- this method assumes that the ellipse is centered on (0,0)
/// Ellipse equation is x*x/(a*a) + (z*z)/(b*b) = 1
/// @returns angle for point on ellipse
///////////////////////////////////////////////////////////////////////////////
float ReconstructionEllipse::getAngle(const float x, const float a, const float z, const float b){

  float theta = acos(x/a);
  if( z < 0.0) theta = 2*M_PI - theta;
  return theta;
}


///////////////////////////////////////////////////////////////////////////////
/// convertPointsToAngles
///////////////////////////////////////////////////////////////////////////////
/// Creates a list of angles to represent a list of points on an ellipse.
/// @param points is the semi-major axis
/// @param ellipse_a is the semi-major axis
/// @param ellipse_x_min is the offset in the x direction of minimal edge of the ellipse
///
/// Ellipse equation is x*x/(a*a) + (z*z)/(b*b) = 1
/// @returns vector of angles via parameter angles
///////////////////////////////////////////////////////////////////////////////
void ReconstructionEllipse::convertPointsToAngles(const vector<PGVector3> &points,
      float ellipse_a,
      float ellipse_x_min,
      vector<float> &angles){

  vector<PGVector3>::const_iterator itr = points.begin();
  for( ; itr != points.end(); ++itr){
    float centered_x = itr->x - ellipse_a - ellipse_x_min;
    float angle = getAngle(centered_x, ellipse_a, itr->z, 0.0);
    angles.push_back(angle);
  }
  std::sort(angles.begin(),angles.end());
  return;
}
///////////////////////////////////////////////////////////////////////////////



void ReconstructionEllipse::print(){

  printf("-------ReconstructionEllipse---------\n");
  printf("cone apex(%.2f, %.2f, %.2f)\n",
              compton_scatter_.getConeApex().x, compton_scatter_.getConeApex().y, compton_scatter_.getConeApex().z);
  printf("cone axis(%.2f, %.2f, %.2f)\n", compton_scatter_.getConeAxis().x, compton_scatter_.getConeAxis().y, compton_scatter_.getConeAxis().z);


  PGVector3 top_intercept_inner = this->getConeInnerYPlaneIntercept(this->phantom_volume_->y_max);
  printf("\ny top inner intercept(%.2f, %.2f, %.2f)\n", top_intercept_inner.x, top_intercept_inner.y, top_intercept_inner.z);
  PGVector3 top_intercept_axis = this->getConeAxisYPlaneIntercept(this->phantom_volume_->y_max);
  printf("y top axis intercept (%.2f, %.2f, %.2f)\n", top_intercept_axis.x, top_intercept_axis.y, top_intercept_axis.z);
  PGVector3 top_intercept_outer = this->getConeOuterYPlaneIntercept(this->phantom_volume_->y_max);
  printf("y top outer intercept(%.2f, %.2f, %.2f)\n\n", top_intercept_outer.x, top_intercept_outer.y, top_intercept_outer.z);

  PGVector3 middle_intercept_inner = this->getConeInnerYPlaneIntercept(0.0);
  printf("\ny mid inner intercept(%.2f, %.2f, %.2f)\n", middle_intercept_inner.x, middle_intercept_inner.y, middle_intercept_inner.z);
  PGVector3 middle_intercept_axis = this->getConeAxisYPlaneIntercept(0.0);
  printf("y mid axis intercept (%.2f, %.2f, %.2f)\n", middle_intercept_axis.x, middle_intercept_axis.y, middle_intercept_axis.z);
  PGVector3 middle_intercept_outer = this->getConeOuterYPlaneIntercept(0.0);
  printf("y mid outer intercept(%.2f, %.2f, %.2f)\n\n", middle_intercept_outer.x, middle_intercept_outer.y, middle_intercept_outer.z);

  PGVector3 bottom_intercept = this->getConeInnerYPlaneIntercept(this->phantom_volume_->y_min);
  printf("y bottom inner intercept(%.2f, %.2f, %.2f)\n", bottom_intercept.x, bottom_intercept.y, bottom_intercept.z);
  bottom_intercept = this->getConeAxisYPlaneIntercept(this->phantom_volume_->y_min);
  printf("y bottom axis intercept(%.2f, %.2f, %.2f)\n", bottom_intercept.x, bottom_intercept.y, bottom_intercept.z);
  bottom_intercept = this->getConeOuterYPlaneIntercept(this->phantom_volume_->y_min);
  printf("y bottom outer intercept(%.2f, %.2f, %.2f)\n", bottom_intercept.x, bottom_intercept.y, bottom_intercept.z);

  printf("gamma likely origin(%.2f, %.2f, %.2f)\n", likely_origin_.x, likely_origin_.y,likely_origin_.z);
  printf("phi (%.2f)\n", phi_*180.0/M_PI);
  printf("alpha (%.2f)\n", alpha_*180.0/M_PI);
  printf("theta1 (%.2f)\n", theta1_*180.0/M_PI);
  printf("theta2 (%.2f)\n", theta2_*180.0/M_PI);
  printf("--------------------------------------\n");
  std::cout<<std::endl;

};
