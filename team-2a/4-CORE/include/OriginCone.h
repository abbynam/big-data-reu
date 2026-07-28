#ifndef ORIGIN_CONE_H_
#define ORIGIN_CONE_H_

#include <vector>
#include "PGVector3.h"
#include "OriginCone.h"
#include "ConicSection.h"

using namespace std;
namespace prompt_gamma_reconstruction{

/*! \brief Class of minimum data needed to describe origin cone.
 *
 *  * 
 * @author Dennis Mackin
 * @data 2018-07-17
 */
    class OriginCone {

    public:
        PGVector3 apex;
        PGVector3 axis;
        float cos_angle;

        OriginCone(){/*cout << "constructing OriginCone (default)"<<endl;*/};
        OriginCone(const PGVector3 &ap, const PGVector3 &ax, const float &cos_ag):apex(ap), axis(ax), cos_angle(cos_ag){};
        static vector<OriginCone> build_origin_cone_array(const vector<ConicSection> &cs){
            vector<OriginCone> oc(cs.size());

            for(size_t i = 0; i < cs.size(); ++i){
                oc[i].apex = cs[i].getScatterInfo()->getConeApex();
                oc[i].axis = cs[i].getScatterInfo()->getConeAxis().normalize();
                oc[i].cos_angle = cos(cs[i].getScatterInfo()->getConeOpeningAngle());

                PGVector3 test_point(0.0, 0.0, 0.0);
                float dca_new = oc[i].get_DCA(test_point);
                float dca_old = cs[i].getDistanceToPoint(test_point);
                if( (dca_new - dca_old)*(dca_new - dca_old) > 1.0E-0){
                    cout<< "ERROR: DCA does not agree: " << dca_new << ", " << cs[i].getDistanceToPoint(test_point) << endl;
                    cout << "Theta1:" << cs[i].getScatterInfo()->getTheta1Degrees() <<", angle:"<< oc[i].cos_angle*180/M_PI <<endl<<endl;
                }
            }

            return oc;
        };

        float get_DCA(const PGVector3 &p) const{

            PGVector3 v1 = p - apex;

            float cos_theta = axis.normalize().dotProduct(v1.normalize());
            cos_theta += (cos_theta < 0.0f)*1.1E-4f - (cos_theta > 0.0f)*1.1E-4f; //floating posize_t error can make abs(cos_theta)>1
            float v1_mag =  v1.magnitude();
            float y = v1_mag * cos_theta;

            float distance_from_cone_axis = sqrt(v1_mag*v1_mag - y*y);

            float radius = y*sqrt(1.0f/(cos_angle*cos_angle) - 1.0f);
            float dca = abs((distance_from_cone_axis - radius) * cos_angle);

            if(!(0.0f <= dca)){
                printf("angle=%.2f, apex=(%.1f, %.1f, %.1f), axis=(%.3f, %.3f, %.3f)\n",
                       cos_angle,       apex.x,apex.y,apex.z,    axis.x, axis.y, axis.z);
                printf("point=(%.1f, %.1f, %.1f)\n", p.x, p.y, p.z);
                printf("dfca=%.1f, radius=%.1f, cos_theta=%.4f, y=%.6f, v1_mag=%.6f, DIFF=%0.6f\n", distance_from_cone_axis, radius, cos_theta, y, v1_mag, v1_mag*v1_mag - y*y);
                cout << endl;
                throw runtime_error("ERROR: DCA is not greater than or equal to 0.");
            }

          return dca;
        }

    PGVector3 get_PCA(const PGVector3 &p) const{

        PGVector3 v1 = p - apex;
        float cos_theta = axis.normalize().dotProduct(v1.normalize());
        cos_theta += (cos_theta < 0.0f)*1E-4 - (cos_theta > 0.0f)*1E-4; //floating posize_t error can make abs(cos_theta)>1
        float v1_mag =  v1.magnitude();
        float y = v1_mag * cos_theta;
        float radius = y*sqrt(1.0f/(cos_angle*cos_angle) - 1.0f);
        float distance_from_cone_axis = sqrt(v1_mag*v1_mag - y*y);
        float dca = abs((distance_from_cone_axis - radius) * cos_angle);

        if(!(0.0f <= dca)){
            printf("angle=%.2f, apex=(%.1f, %.1f, %.1f), axis=(%.3f, %.3f, %.3f)\n",
                   cos_angle,       apex.x,apex.y,apex.z,    axis.x, axis.y, axis.z);
            printf("point=(%.1f, %.1f, %.1f)\n", p.x, p.y, p.z);
            cout << endl;
            throw runtime_error("ERROR in OriginCone::get_PCA: DCA is not greater than or equal to 0.");
        }

        PGVector3 p0 = apex + axis * (v1_mag/cos_theta);
        PGVector3 vector_to_cone_axis = (p0 - p).normalize();
        PGVector3 pca = p + vector_to_cone_axis * dca;

        return pca;
    }
  };
}//end of namespace prompt_gamma_reconstruction
#endif //ORIGIN_CONE_H_
