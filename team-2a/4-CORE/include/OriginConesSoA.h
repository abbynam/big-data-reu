#ifndef ORIGIN_CONES_SOA_H_
#define ORIGIN_CONES_SOA_H_

#include <vector>
#include "PGVector3.h"
#include "OriginCone.h"

namespace prompt_gamma_reconstruction{

/*! \brief Structure of array for cone information
 *     SoA format is preferred for SIMD (same instruction, multiple data) and GPU
 *     code.
 *
 *  *
 * @author Dennis Mackin
 * @data 2018-10-01
 */
    class OriginConesSoA {

    public:
        std::vector<float> apex_x;
        std::vector<float> apex_y;
        std::vector<float> apex_z;

        std::vector<float> axis_x;
        std::vector<float> axis_y;
        std::vector<float> axis_z;

        std::vector<float> cos_angle;

    void push_back(const OriginCone &c){
        apex_x.push_back(c.apex.x);
        apex_y.push_back(c.apex.y);
        apex_z.push_back(c.apex.z);

        axis_x.push_back(c.axis.x);
        axis_y.push_back(c.axis.y);
        axis_z.push_back(c.axis.z);
        cos_angle.push_back(c.cos_angle);
    };

    void add_cones(const vector<OriginCone> &c){
        for(size_t i=0; i<c.size(); ++i) push_back(c[i]);
    }

    vector<float> get_arranged_memory() const{
        const size_t num_float_vectors = 7;
        const size_t N = apex_x.size();
        vector<float> m(num_float_vectors * N);
        for(size_t i=0; i < N; ++i){
            m[i]     = apex_x[i];
            m[i + 1*N] = apex_y[i];
            m[i + 2*N] = apex_z[i];
            m[i + 3*N] = axis_x[i];
            m[i + 4*N] = axis_y[i];
            m[i + 5*N] = axis_z[i];
            m[i + 6*N] = cos_angle[i];
        }
        return m;
    };

    size_t getLength() const{
        return apex_x.size();
    }

  };
}//end of namespace prompt_gamma_reconstruction
#endif //ORIGIN_CONES_SOA_H_
