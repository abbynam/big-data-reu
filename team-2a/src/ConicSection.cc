#include "ConicSection.h"
#include "PhantomVolume.h"
//#include "RandomPointOnCircleSingleton.h"

using namespace std;
//using namespace prompt_gamma_reconstruction;

namespace prompt_gamma_reconstruction{

ConicSection::ConicSection(ComptonScatter comptonScatter, shared_ptr<const PhantomVolume> phantomVolume, const size_t seed)
    : alpha_(0), cos_alpha_(0), sin_alpha_(0), phi_(0), theta1_(0), theta2_(0), tan_phi_(0), cos_phi_(0),
        cos_theta1_(0), sin_theta1_(0), weight_(0),
        y_max_cone_axis_frame_(0), y_min_cone_axis_frame_(0), is_y_range_set_(0), compton_scatter_(comptonScatter),
        phantom_volume_(phantomVolume), rand_(seed), inverse_square_param_(1.0){

    ptr_random_sqrt_ = RandomSqrtSingleton::Instance();
    ptr_random_ = RandomSingleton::Instance();

    isInitialOriginSet_ = false;
    isLikelyOriginSet_ = false;
    isPhantomParallelToConeFrame_ = false;

    alpha_ = comptonScatter.getAlpha();
    cos_alpha_ = cos(alpha_);
    sin_alpha_ = sin(alpha_);

    phi_ = comptonScatter.getScatteringAngle();
    tan_phi_ = tan(phi_);
    cos_phi_ = cos(phi_);

    theta1_ = calculateXZrotationAngle(compton_scatter_.getConeAxis().x, compton_scatter_.getConeAxis().z);
    cos_theta1_ = cos(theta1_);
    sin_theta1_ = sin(theta1_);

    theta2_ = alpha_ - M_PI/2.0;


    setWeight(1.0);

    mc_truth = shared_ptr<MC_Truth> (new MC_Truth);
};

    void ConicSection::setYrange_(){
        auto y_range = this->getPhantomVolumeYRange();

        y_min_cone_axis_frame_ = y_range.first;
        y_max_cone_axis_frame_ = y_range.second;

        is_y_range_set_ = true;

    };

    ///Find shortest distance from cone and line defined by p1 and p2
    pair<float, PGVector3> ConicSection::getDistanceToLine(const PGVector3 p1, const PGVector3 p2) const{

        pg_tools::Random rand;
        const float EPSILON = 0.5;
        auto p = p1 + (p2 - p1) * rand.Rndm();
        auto v1 = (p2 - p1).normalize();
        auto pca = getPointOfClosestApproach(p);

        auto v2 = (pca - p);
        auto cos_theta = v1.dotProduct(v2.normalize());
        auto length_projected_onto_v1 = cos_theta * v2.magnitude();
        auto p_new = p + v1 * length_projected_onto_v1;
        while(p_new.getDistanceToPoint(p) > EPSILON){
            p = p_new;
            pca = getPointOfClosestApproach(p);
            v2 = pca - p;
            cos_theta = v1.dotProduct(v2.normalize());
            length_projected_onto_v1 = cos_theta * v2.magnitude();
            p_new = p + v1 * length_projected_onto_v1;
        }
        auto dca = pca.getDistanceToPoint(p_new);
        return make_pair(dca, pca);
    }


///////////////////////////////////////////////////////////////////////////////
float ConicSection::getDistanceToPoint(PGVector3 point) const{
///////////////////////////////////////////////////////////////////////////////
/// @param point in phantom frame for which the distance is to be calculated.
///
/// @returns the smallest distance between the point and the cone (distance of closest approach)
///////////////////////////////////////////////////////////////////////////////

    transformPointToConeAxisFrameFromPhantomFrame(point);
    float distance_from_y_axis = sqrt(point.x*point.x + point.z*point.z);
    float radius = point.y*tan_phi_;
    float distance = abs(distance_from_y_axis - radius)*cos_phi_;
    if(distance != distance){
        stringstream ss;
        ss<<"ERROR: Distance to point is nan for point("<<point.x<<","<<point.y<<","<<point.z<<").\n";
        cout<<ss.str()<<endl;
        throw runtime_error(ss.str());
    }
    return distance;
}


///////////////////////////////////////////////////////////////////////////////
float ConicSection::getDistanceToPoint(const PGVector3 &scatter1, const PGVector3 &scatter2, const float scatter_angle, const PGVector3 &point){
///////////////////////////////////////////////////////////////////////////////
/// Uses the scatter positions and the scatter angle to calculate the DCA.
///
///
/// @returns the smallest distance between the point and the cone (distance of closest approach)
///////////////////////////////////////////////////////////////////////////////

    PGVector3 v2 = scatter2 - scatter1;
    PGVector3 v1 = point - scatter1;

    float cos_theta = v1.dotProductNormalized(v2);
    float theta = M_PI - acos(cos_theta);
    float y = v1.magnitude() * cos_theta;
    float radius = y * tan(scatter_angle);
    float distance_from_cone_axis = y * tan(theta);
    float dca = abs((distance_from_cone_axis - radius) * cos(scatter_angle));

    assert(dca >= 0.0);

    return dca;
}


vector<float> ConicSection::getDistanceToPoints(vector<PGVector3> &points) const{
///////////////////////////////////////////////////////////////////////////////
/// @param point in phantom frame for which the distance is to be calculated.
///
/// @returns the smallest distance between the point and the cone (distance of closest approach)
///////////////////////////////////////////////////////////////////////////////

  transformPointsToConeAxisFrameFromPhantomFrame(points);

  auto numberOfPoints = points.size();
  vector<float> distances(numberOfPoints);
  vector<float> distances_from_y_axis(numberOfPoints);
  vector<float> radi(numberOfPoints);

  for(auto i=0u; i<numberOfPoints; ++i){
    distances_from_y_axis[i] = sqrt(points[i].x*points[i].x + points[i].z*points[i].z);
    radi[i] = points[i].y*tan_phi_;
    distances[i] = abs(distances_from_y_axis[i] - radi[i])*cos_phi_;
  }
  return distances;
}

///////////////////////////////////////////////////////////////////////////////
PGVector3 ConicSection::getPointOfClosestApproach(PGVector3 point) const{
///////////////////////////////////////////////////////////////////////////////
/// @param point in phantom frame for which the distance is to be calculated.
///
/// @returns the smallest distance between the point and the cone
///////////////////////////////////////////////////////////////////////////////


  transformPointToConeAxisFrameFromPhantomFrame(point);
  PGVector3 pca; //point of closest approach

  float distance_from_y_axis = sqrt(point.x*point.x + point.z*point.z);
  float radius = point.y*tan_phi_;
  float dca = abs(distance_from_y_axis - radius)*cos_phi_;
  float sin_phi = tan_phi_*cos_phi_;
  float delta_y = dca * sin_phi;

  //if distance_from_y_axis is greater than radius, then the point is outside cone
  // if point is outside cone than delta_y is positive, otherwise point is in cone
  // and delta y is negative
  pca.y = (distance_from_y_axis > radius)? point.y + delta_y : point.y - delta_y;
  float radius_pca =  pca.y*tan_phi_;

  float theta = M_PI * 0.5f;
  if(1E-10 < point.z*point.z) theta = atan(point.x/point.z);
  pca.x = sin(theta) * radius_pca;
  pca.z = cos(theta) * radius_pca;


  //Ensure that pca and the point are the same sign in the cone frame
  if(point.x * pca.x < 0.0) pca.x *= -1.0;
  if(point.z * pca.z < 0.0) pca.z *= -1.0;

  transformPointToPhantomFrameFromConeAxisFrame(pca);

  return pca;
}


///Calculates the angle needed to rotate the
/// y=apex.y axis in order to get the cone axis to lie
/// along the new X axis so that the conic section
/// major axes will be in the x direction.
float ConicSection::calculateXZrotationAngle(float x, float z){
        x += 1.0E-10; //Fix case where cone axis is parallel to y axis
        float vx = x/sqrt(x*x + z*z);
        float phi = acos(vx);

        //If Z is positive, reverse the rotation. Avoid if branch.
        phi *= (z <= 0) - (z > 0);

        return phi;
}

void ConicSection::transformPointToConeAxisFrameFromPhantomFrame_orig(PGVector3 &point) const{
  //translate to cone axis
  point -= compton_scatter_.getConeApex();
  point.rotateYaxis(theta1_);
  point.rotateZaxis(- alpha_);
  point.y *= -1.0;
}

void ConicSection::transformPointsToConeAxisFrameFromPhantomFrame(vector<PGVector3> &points) const{

    auto numberOfPoints = points.size();
    auto apex = compton_scatter_.getConeApex();

    vector<float> tmp(numberOfPoints);

    float cos_theta1 = cos_theta1_;
    float sin_theta1 = sin_theta1_;
    float cos_alpha = cos_alpha_;
    float sin_alpha = sin_alpha_;


    size_t i = 0u;
    //SLOWS-->#pragma omp parallel for
    for(i = 0u; i<numberOfPoints; ++i){
        points[i].x -= apex.x;
        points[i].y -= apex.y;
        points[i].z -= apex.z;

        tmp[i] = points[i].x*cos_theta1 - points[i].z*sin_theta1;
        points[i].z = points[i].x*sin_theta1 + points[i].z*cos_theta1;
        points[i].x = tmp[i];

        tmp[i] = points[i].x*cos_alpha + points[i].y*sin_alpha;
        points[i].y = - points[i].x*sin_alpha + points[i].y*cos_alpha;
        points[i].x = tmp[i];

        points[i].y *= -1.0;
    }
}


void ConicSection::transformPointToConeAxisFrameFromPhantomFrame(PGVector3 &point) const{
    //translate to cone axis
    point -= compton_scatter_.getConeApex();

    //rotate around y to make z=0 and x = combination of phantom x&z directions
    auto tmp = point.x*cos_theta1_ - point.z*sin_theta1_;
    point.z = point.x*sin_theta1_ + point.z*cos_theta1_;
    point.x = tmp;

    //rotate around z axis so that Y is now aligned with cone axis
    tmp = point.x*cos_alpha_ - point.y*sin_alpha_;
    point.y = point.x*sin_alpha_ + point.y*cos_alpha_;
    point.x = tmp;

    if(point.x != point.x){
        stringstream ss;
        ss<<"ERROR: Transformed point is nan for point("<<point.x<<","<<point.y<<","<<point.z<<").\n";
        cout<<ss.str()<<endl;
        throw runtime_error(ss.str());
    }


}


void ConicSection::transformPointToPhantomFrameFromConeAxisFrame(PGVector3 &point) const{

    //rotate around z axis
    auto tmp = point.x*cos_alpha_ + point.y*sin_alpha_;
    point.y = -point.x*sin_alpha_ + point.y*cos_alpha_;
    point.x = tmp;

    //rotate around y axis
    tmp = point.x*cos_theta1_ + point.z*sin_theta1_;
    point.z = -point.x*sin_theta1_ + point.z*cos_theta1_;
    point.x = tmp;

    point += compton_scatter_.getConeApex();
}


PGVector3 ConicSection::getConeAxisYPlaneIntercept(float plane_y){

  PGVector3 output;
  PGVector3 axis = getComptonScatter().getConeAxis();
  PGVector3 apex = getComptonScatter().getConeApex();
  float delta_y = apex.y - plane_y;
  output.x = apex.x - delta_y *( axis.x/axis.y);
  output.z = apex.z - delta_y *( axis.z/axis.y);
  output.y = plane_y;

  return output;
};


pair<float, float> ConicSection::getPhantomVolumeYRange(){
    auto apex = this->getComptonScatter().getConeApex();

    auto phantom_center = phantom_volume_->get_center_point();
    transformPointToConeAxisFrameFromPhantomFrame(phantom_center);
    float y_start = phantom_center.y;

    float y_step_size = 10.0;
    float search_stop_near = y_step_size;
    float search_stop_far = apex.magnitude() + 2.0 * phantom_volume_->get_max_length();


    float y_range_min = search_stop_far;
    float y_range_max = search_stop_near;

    const size_t TEST_ANGLES = 60;
    PGVector3 test_point;

    float angle = 0;
    bool values_set = false;
    bool layer_in_volume = false;

    vector<float> testpoints_x(TEST_ANGLES, 0.0);
    vector<float> testpoints_z(TEST_ANGLES, 0.0);
    for(size_t i = 0; i < TEST_ANGLES; ++i) {
        angle = (1.0 - 2.0)*(i % 2) * i*2.0*M_PI/TEST_ANGLES;
        testpoints_x[i] = tan_phi_*cos(angle);
        testpoints_z[i] = tan_phi_*sin(angle);
    }

    float y = y_start;
    auto search_layer = [&](float y){
        for(size_t i=0; i<TEST_ANGLES + 1; ++i){

            test_point.x = y*testpoints_x[i];
            test_point.y = y;
            test_point.z = y*testpoints_z[i];
            transformPointToPhantomFrameFromConeAxisFrame(test_point);

            if(phantom_volume_->is_in_volume(test_point)){
                y_range_min = min(y_range_min, y);
                y_range_max = max(y_range_max, y);
                values_set = true;
                layer_in_volume = true;
                break;
            }
        }
    };

    //search from middle away from cone appex
    while(y < search_stop_far){
        layer_in_volume = false;
        search_layer(y);
        if(values_set && !layer_in_volume) break; //We've found the volume, and now we are out of the farside of phantom
        y += y_step_size;
    }

    //search from middle toward cone apex
    y = y_start - y_step_size;
    while(y > search_stop_near){
        layer_in_volume = false;
        search_layer(y);
        if(values_set && !layer_in_volume) break; //We've found the volume, and now we are out of the top of phantom
        y -= y_step_size;
    }

    return make_pair(y_range_min - 0.5*y_step_size, y_range_max + 0.5*y_step_size);
};


///////////////////////////////////////////////////////////////////////////////
float ConicSection::solveConeEquation(float x, float y, float phi){
///////////////////////////////////////////////////////////////////////////////
/// @param x the independent coordinate orthogonal to the cone axis
/// @param y the coordinate in the direction of the cone axis
/// @param phi angle between cone axis and cone surface
/// @returns z the dependent coordinate orthogonal to the cone axis
///////////////////////////////////////////////////////////////////////////////
    //now solve cone equation to get zprimeprime
    float r = y*tan(phi);
    if ( 0 > r*r - x*x ){
        stringstream error_message;
        error_message <<"Error in ConicSection::solveConeEquation\n"
                    <<"Complex root indicates that there is no solution to to the cone equation"
                    "of the form ("<<x<<","<<y<<",z).\n";
        throw runtime_error(error_message.str()) ;
    }
    return sqrt(r*r - x*x);
}


long ConicSection::getRandomPointInPhantom(PGVector3 &random_point, size_t max_num_tries){

    if(! is_y_range_set_) setYrange_();

    for(size_t iNumTries=0; iNumTries<max_num_tries; ++iNumTries){

        // DSM 20151230 This is the original Y randomization code for SOE

        //DSM 2016-04-06 This version evenly distributes the randoms on the cone
//        random_point.y = sqrt(
//                rand_.Rndm() *
//                        (y_max_cone_axis_frame_*y_max_cone_axis_frame_ - y_min_cone_axis_frame_*y_min_cone_axis_frame_)
//                + y_min_cone_axis_frame_*y_min_cone_axis_frame_
//        );
//
//        random_point.y =  pow(rand_.Rndm(),0.5) * (y_max - y_min) + y_min

//        random_point.y = rand_.Rndm() *
//                (y_max_cone_axis_frame_*y_max_cone_axis_frame_ - y_min_cone_axis_frame_*y_min_cone_axis_frame_)
//                + y_min_cone_axis_frame_;
        //DSM - this version includes an inverse square correction
//        random_point.y = sqrt(
//                pow(rand_.Rndm(),inverse_square_param_) *
//                (y_max_cone_axis_frame_*y_max_cone_axis_frame_ - y_min_cone_axis_frame_*y_min_cone_axis_frame_)
//                + y_min_cone_axis_frame_*y_min_cone_axis_frame_
//        );
        //DSM the inverse square correction is now part of the point algorithm
        random_point.y = sqrt(rand_.Rndm()*
                (y_max_cone_axis_frame_*y_max_cone_axis_frame_ - y_min_cone_axis_frame_*y_min_cone_axis_frame_)
                + y_min_cone_axis_frame_*y_min_cone_axis_frame_
        );

        float random_angle = 2.0*M_PI*rand_.Rndm();
        float r = random_point.y *this->tan_phi_;
        random_point.x = r*cos(random_angle);
        random_point.z = r*sin(random_angle);


        transformPointToPhantomFrameFromConeAxisFrame(random_point);

        if(phantom_volume_->is_in_volume(random_point)){
            if(!isLikelyOriginSet_) setLikelyOrigin(random_point); //set only if 1st point found
            return ++iNumTries;
        }
    }
    return -1;
}

//
//    int ConicSection::getRandomPointInPhantom(PGVector3 &random_point, int max_num_tries){
//        PGVector3 tmp;
//
//        for( int iNumTries=0; iNumTries<max_num_tries; ++iNumTries){
//            tmp.x = rand_.Rndm()*(this->phantom_volume_->x_max - this->phantom_volume_->x_min) + this->phantom_volume_->x_min;
//            tmp.y = rand_.Rndm()*(this->phantom_volume_->y_max - this->phantom_volume_->y_min) + this->phantom_volume_->y_min;
//            tmp.z = rand_.Rndm()*(this->phantom_volume_->z_max - this->phantom_volume_->z_min) + this->phantom_volume_->z_min;
//
//            random_point = this->getPointOfClosestApproach(tmp);
//
//            if(phantom_volume_->is_in_volume(random_point)){
//                if(!isLikelyOriginSet_) setLikelyOrigin(random_point); //set only if 1st point found
//                return ++iNumTries;
//            }
//        }
//
//        return -1;
//    }


size_t ConicSection::getRandomStepInPhantom(const PGVector3 &current_point, PGVector3 &new_point, float step_size, size_t max_num_tries){

  PGVector3 tmpPoint;
  step_size *= 1.0/sqrt(3.0); ///Used to make 3D steps approximate step size in length

  for(size_t iNumTries=0; iNumTries<max_num_tries; ++iNumTries){
    tmpPoint = current_point;
    tmpPoint.x += rand_.Gaus(0.0, step_size);
    tmpPoint.y += rand_.Gaus(0.0, step_size);
    tmpPoint.z += rand_.Gaus(0.0, step_size);
    new_point = getPointOfClosestApproach(tmpPoint);


    if(phantom_volume_->is_in_volume(new_point)){
        if(!isLikelyOriginSet_) setLikelyOrigin(new_point);
        return ++iNumTries;
    }
  }
  return -1;
}



void ConicSection::setMCTruth(float *scattering_angle, PGVector3 *position, float *energy_deposition, float *incident_energy, float initial_energy, float *origin){

  mc_truth->scattering_angle = vector<float>(scattering_angle, scattering_angle+2);
  mc_truth->position = vector<PGVector3>(position, position+3);
  mc_truth->energy_deposition = vector<float>(energy_deposition, energy_deposition+3);
  mc_truth->incident_energy = vector<float>(incident_energy, incident_energy+3);
  mc_truth->initial_energy = initial_energy;
  mc_truth->origin_true = vector<float>(origin, origin+3);
};


void ConicSection::setLikelyOrigin(const PGVector3 point){

    likely_origin_.x = point.x;
    likely_origin_.y = point.y;
    likely_origin_.z = point.z;

    isLikelyOriginSet_ = true;

    if( false == isInitialOriginSet_){
        setInitialOrigin(point);
        isInitialOriginSet_ = true;
    }
};

 void ConicSection::setInitialOrigin(const PGVector3 point){

    initial_origin_.x = point.x;
    initial_origin_.y = point.y;
    initial_origin_.z = point.z;
};

void ConicSection::print(){
    printf("-------ReconstrunctionEllipse---------\n");
    printf("cone apex(%.2f, %.2f, %.2f)\n",
           compton_scatter_.getConeApex().x, compton_scatter_.getConeApex().y, compton_scatter_.getConeApex().z);
    printf("cone axis(%.2f, %.2f, %.2f)\n",
           compton_scatter_.getConeAxis().x, compton_scatter_.getConeAxis().y, compton_scatter_.getConeAxis().z
    );
    printf("phi (%.2f)\n", phi_);
    printf("alpha (%.2f)\n", alpha_);
    printf("theta1 (%.2f)\n", theta1_);
    printf("theta2 (%.2f)\n", theta2_);
    printf("--------------------------------------\n");
};

}
