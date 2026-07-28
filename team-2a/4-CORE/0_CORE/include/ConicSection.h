#ifndef CONIC_SECTION_H_
#define CONIC_SECTION_H_

//#include <assert.h>
#include <cassert>
#include <stdexcept>

#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>

#include <memory>


//prompt gamma includes
#include "ComptonScatter.h"
#include "PhantomVolume.h"
#include "RandomSqrtSingleton.h"
#include "RandomSingleton.h"
//#include "RandomPointOnCircleSingleton.h"
#include "PGVector3.h"
#include "TripleScatter.h"
#include "utilities/Random.h"


namespace prompt_gamma_reconstruction {

///Structure used to store true MC values so that they can be compared to the values smeared with detector effects
/*! \brief Monte Carlo calculated gamma emission positions.
 * 
 * The original gamma emission positions calculated by the 
 * Monte Carlo are stored so that they can be compared to the position
 * smeared to simulate detector effects and the positions
 * determined by the reconstruction algorithm.
 * 
 * @author Dennis Mackin
 */
    class MC_Truth {

    public:
        std::vector<float> scattering_angle;
        std::vector<PGVector3> position;
        std::vector<float> energy_deposition;
        std::vector<float> incident_energy;
        float initial_energy;
        std::vector<float> origin_true;


        MC_Truth() : scattering_angle(2), position(3), energy_deposition(2),
                     incident_energy(3), initial_energy(0.0), origin_true(3) {/* DO NOTHING */};

        void print() const {
            printf("----- MC Truth -----\n");
            printf("scatter angles: %.3f, %.3f\n", scattering_angle[0], scattering_angle[1]);
            printf("scatter positions: %s, %s, %s\n", position[0].print().c_str(), position[1].print().c_str(),
                   position[2].print().c_str());
            printf("energy deposition: %.3f, %.3f, %.3f\n", energy_deposition[0], energy_deposition[1],
                   energy_deposition[2]);
            printf("incident energy: %.3f, %.3f, %.3f\n", incident_energy[0], incident_energy[1], incident_energy[2]);
            printf("initial energy: %.3f\n", initial_energy);
            printf("origin: %.3f, %.3f, %.3f\n", origin_true[0], origin_true[1], origin_true[2]);
            printf("---- END MC Truth ---\n\n");
        };
    };


/*! \brief Base class of the reconstruction elements ellipse and parabola.
 * 
 * The intersection of cones and planes form either ellipses or 
 * parabolas. This base class is for both reconstruction ellipses 
 * and reconstruction parabolas. ConicSection objects store the
 * cone specific data. This class also defines the affine transformations
 * to and from this cone's reference frame and the phantom reference frame.
 * 
 * @author Dennis Mackin
 */
    class ConicSection {

    private:
        RandomSqrtSingleton *ptr_random_sqrt_;
        RandomSingleton *ptr_random_;

    protected:
        //values can be calculated once and stored
        float alpha_; ///angle between cone axis and line from apex to xz plane
        float cos_alpha_;
        float sin_alpha_;

        float phi_; ///scattering angle and angle between cone axis and cone surface
        float theta1_; ///rotation angle making x axis the axis of symmetry or major axis
        float theta2_; ///rotation angle for rotation around z' axis making y'' the cone axis
        float tan_phi_;
        float cos_phi_;

        float cos_theta1_;
        float sin_theta1_;

        float weight_;

        float y_max_cone_axis_frame_;
        float y_min_cone_axis_frame_;
        bool is_y_range_set_;

        float inverse_square_param_;
        pg_tools::Random rand_;


        ///structure mc_truth stores a pointer to the true MC values so that
        /// comparisons can be made between the smeared and unsmeared values
        //shared_ptr<MC_Truth> mc_truth;
        //shared_ptr<Scatter> scatter_;
        std::shared_ptr<MC_Truth> mc_truth;
        std::shared_ptr<Scatter> scatter_;


        bool isLikelyOriginSet_;
        PGVector3 likely_origin_; ///best guess as to origin of the related gamma

        bool isInitialOriginSet_;
        PGVector3 initial_origin_; ///initial guess as to origin of the related gamma

        ComptonScatter compton_scatter_;
        shared_ptr<const PhantomVolume> phantom_volume_;
        Rectangle rotated_phantom_boundary_;

        void setPhantomVolume_(shared_ptr<const PhantomVolume> phantomVolume) { phantom_volume_ = phantomVolume; };

        void getXZprimeBoundaries();
        bool isPhantomParallelToConeFrame_;
        void set_y_min_cone_axis_frame(float y) { y_min_cone_axis_frame_ = y; };
        void set_y_max_cone_axis_frame(float y) { y_max_cone_axis_frame_ = y; };
        void setYrange_();


    public:

        //accessors
        inline float getAlpha() const { return alpha_; };

        inline float getPhi() const { return phi_; };

        inline float getTheta1() const { return theta1_; };

        inline float getTheta2() const { return theta2_; };

        inline float getWeight() const { return weight_; };

        inline void setWeight(const float weight) { weight_ = weight; };

        inline float get_y_max_cone_axis_frame() { return y_max_cone_axis_frame_; };

        inline float get_y_min_cone_axis_frame() { return y_min_cone_axis_frame_; };

        void setMCTruth(float *scattering_angle, PGVector3 *position, float *energy_deposition, float *incident_energy,
                        float initial_energy, float *origin);


        //Creation methods
        ConicSection(ComptonScatter comptonScatter, shared_ptr<const PhantomVolume> phantomVolume, const size_t seed);

        virtual ~ConicSection() {};

        void setInitialOrigin(const PGVector3 point);

        virtual bool doesConicIntersectPhantom() { return true; };

        virtual size_t getNumberOfPhantomIntercepts(float y) { return 0; };

        inline void setScatter(shared_ptr<Scatter> s) { scatter_ = s; };

        inline void setInverseSquareParam(float p) { inverse_square_param_ = p; };

        inline const ComptonScatter &getComptonScatter() const { return compton_scatter_; };

        inline const MC_Truth &getMCTruth() { return *mc_truth; };

        inline const shared_ptr<Scatter> &getScatterInfo() const { return scatter_; };


        /// returns a random posize_t that is on the cone surface and in the
        /// phantom volume.
        long getRandomPointInPhantom(PGVector3 &randPoint, size_t num_random_tries = 1000);

        /// alternative random walk approach to finding the next point.
        size_t getRandomStepInPhantom(const PGVector3 &current_point, PGVector3 &new_point, float step_size,
                                   size_t max_num_tries); ///Find new posize_t on cone based on old point


        //Transformation methods
        void transformPointToConeAxisFrameFromPhantomFrame(PGVector3 &point) const;

        void transformPointsToConeAxisFrameFromPhantomFrame(vector<PGVector3> &points) const;

        void transformPointToConeAxisFrameFromPhantomFrame_orig(PGVector3 &point) const;

        void transformPointToPhantomFrameFromConeAxisFrame(PGVector3 &point) const;

        PGVector3 getConeAxisYPlaneIntercept(float y);

        pair<float, float> getPhantomVolumeYRange();

        static float calculateXZrotationAngle(float x, float z);

        //DCA and PCA are two crucial concepts and methods
        // for the reconstruction of Compton Camera images.
        float getDistanceToPoint(const PGVector3 point) const;
        static float getDistanceToPoint(const PGVector3 &scatter1, const PGVector3 &scatter2, const float angle, const PGVector3 &point);

        //Find shortest distance from cone and line z= ax + by +
        pair<float, PGVector3> getDistanceToLine(const PGVector3 p1, const PGVector3 p2) const;

        vector<float> getDistanceToPoints(vector<PGVector3> &points) const;


        PGVector3 getPointOfClosestApproach(const PGVector3 point) const;

        //Algorithm methods
        inline const Rectangle getRotatedPhantomBoundary() { return rotated_phantom_boundary_; };

        static float solveConeEquation(float x, float y, float phi);

        bool isPointOnCone(PGVector3 point);

        const inline PGVector3 &getLikelyOrigin() const { return likely_origin_; };

        inline PGVector3 &getInitialOrigin() { return initial_origin_; };

        void setLikelyOrigin(const PGVector3 point);

        //Output methods
        virtual void print();

        virtual void draw(float y) {/*do nothing here*/};
    };

}

#endif //CONIC_SECTION_H_
