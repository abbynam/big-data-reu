#ifndef DENSITY_MATRIX_H_
#define DENSITY_MATRIX_H_
#define _USE_MATH_DEFINES

//standard C++ includes
#include <cmath>
#include <cstdio>
#include <valarray>
#include <vector>
#include <stdexcept>

//PromptGamma includes
#include "ComptonScatter.h"
#include "PGVector3.h"
#include "DensityEstimator.h"
#include "utilities/Random.h"


using namespace std;
namespace prompt_gamma_reconstruction {

    class DensityMatrix : public DensityEstimator {

    public:
        DensityMatrix(
                double x_min, double x_max, size_t x_bins,
                double y_min, double y_max, size_t y_bins,
                double z_min, double z_max, size_t z_bins,
                bool poisson_binning = 0) :
                x0_min_(x_min), x0_max_(x_max), x0_bins_(x_bins),
                x_min_(x_min), x_max_(x_max), x_bins_(x_bins),
                y0_min_(y_min), y0_max_(y_max), y0_bins_(y_bins),
                y_min_(y_min), y_max_(y_max), y_bins_(y_bins),
                z0_min_(z_min), z0_max_(z_max), z0_bins_(z_bins),
                z_min_(z_min), z_max_(z_max), z_bins_(z_bins),
                num_bins_(x_bins * y_bins * z_bins), counts_vector_(-1.0, 1), poisson_binning_(poisson_binning),
                bad_bin_(-1), rand_(4321) {
            if (0 == x0_bins_) x0_bins_ = 1;
            x_bin_size_reciprocal_ = static_cast<double>(x0_bins_) / (x0_max_ - x0_min_);

            if (0 == y_bins_) y_bins_ = 1;
            y_bin_size_reciprocal_ = static_cast<double>(y0_bins_) / (y0_max_ - y0_min_);

            if (0 == z0_bins_) z0_bins_ = 1;
            z_bin_size_reciprocal_ = static_cast<double>(z0_bins_) / (z0_max_ - z0_min_);

            setup_bins_();
        };

        DensityMatrix() :
                x_min_(-1.0), x_max_(1.0), x_bins_(1), x_bin_size_reciprocal_(1.0),
                y_min_(-1.0), y_max_(1.0), y_bins_(1), y_bin_size_reciprocal_(1.0),
                z_min_(-1.0), z_max_(1.0), z_bins_(1), z_bin_size_reciprocal_(1.0),
                num_bins_(1), counts_vector_(1), poisson_binning_(0), bad_bin_(-1), rand_(4321) {
            cout << "DensityMatrix::DensityMatrix(empty constructor):" << counts_vector_.size() << endl;
        };


        DensityMatrix(const DensityMatrix &other) : x_min_(-1.0), x_max_(1.0), x_bins_(1), x_bin_size_reciprocal_(1.0),
                                                    y_min_(-1.0), y_max_(1.0), y_bins_(1), y_bin_size_reciprocal_(1.0),
                                                    z_min_(-1.0), z_max_(1.0), z_bins_(1), z_bin_size_reciprocal_(1.0),
                                                    num_bins_(1), counts_vector_(1), poisson_binning_(0), bad_bin_(-1),
                                                    rand_(4321) {

            make_copy(other);
        }


        ~DensityMatrix() {
            if (-1 != bad_bin_) {
                cerr << "WARNING: the value for DensityMatrix.bad_bin = " << bad_bin_ << ". It should be equal to -1.\n"
                << "This error indicates that an attempt was made to add density outside the valid region." << endl;
            }
        }

        // functions accepts event position vector and returns cooresponding density value
        inline float getDensity(const PGVector3 &pos) const {
            size_t bin_number = static_cast<size_t>(getBinNumber_(pos));

            if (counts_vector_.size() <= bin_number) {
                stringstream ss;
                ss << "ERROR in DensityMatrix::getDensity(" << pos.x << "," << pos.y << "," << pos.z <<
                ") returned bin " << bin_number << "\n" << counts_vector_.size() - 1 << " is maximum bin.\n";
                ss << "x range(" << this->x_min_ << ", " << x_max_ <<")\n";
                ss << "y range(" << this->y_min_ << ", " << y_max_ <<")\n";
                ss << "z range(" << this->z_min_ << ", " << z_max_ <<")\n";
                ss << "bins(" << x_bins_ <<", " << y_bins_ << ", " << z_bins_ <<")" <<endl;
                ss << "bin_size_reciprocals (" << x_bin_size_reciprocal_ <<", " << y_bin_size_reciprocal_ << ", " << z_bin_size_reciprocal_ <<")" <<endl;
                cout << ss.str() << endl;
                throw runtime_error(ss.str());
            }
            return counts_vector_[bin_number];
        };

        vector<float> getDensities(const vector<PGVector3> &positions) const {
            vector<float> densities(positions.size());
            for(size_t i = 0; i < densities.size(); ++i){
                densities[i] = getDensity(positions[i]);
            }
            return densities;
        };


        inline void fill(const PGVector3 &pos, float weight) {
        #pragma omp atomic
            operator[](pos) += weight;
        };


        inline float &operator[](const PGVector3 &pos) {
            int bin_number = getBinNumber_(pos);
            if (bin_number < 0) {
                stringstream ss;
                ss << "ERROR in DensityMatrix::operator[]: point(" << pos.x << "," << pos.y << "," << pos.z <<
                ") returned bin "
                << bin_number << ".\n";
                cout << ss.str() << endl;
                //print();
                throw runtime_error(ss.str());
            }
            return counts_vector_[bin_number];
        };


        inline void updateMatrix(const PGVector3 &oldPos, const PGVector3 &newPos) {
            updateMatrix(oldPos, newPos, 1.0);
        };

        inline void updateMatrix(const PGVector3 &oldPos, const PGVector3 &newPos, const float &weight) {

            auto bin = getBinNumber_(oldPos);
            #pragma omp atomic
            counts_vector_[bin] -= weight;

            bin = getBinNumber_(newPos);
            #pragma omp atomic
            counts_vector_[bin] += weight;
        };

        DensityMatrix &operator=(float rhs) {
            counts_vector_ = rhs;
            return *this;
        }


        DensityMatrix &operator=(const DensityMatrix &rhs) {
            make_copy(rhs);
            return *this;
        }


        DensityMatrix &operator+=(const DensityMatrix &rhs) {
            if (counts_vector_.size() != rhs.counts_vector_.size()) {
                stringstream ss;
                ss << "ERROR in DensityMatrix &operator+= :\n";
                ss << "matrix size(" << counts_vector_.size() << ") != right hand side(" << rhs.counts_vector_.size() <<
                ").\n\n";
                cout << ss.str() << endl;
                throw runtime_error(ss.str());
            }
            counts_vector_ += rhs.counts_vector_;
            return *this;
        }


        inline bool is_in_volume(const PGVector3 &point) const {
            return (
                    point.x <= x_max_ && point.x >= x_min_
                    && point.y <= y_max_ && point.y >= y_min_
                    && point.z <= z_max_ && point.z >= z_min_);
        };


        void clear() {
            counts_vector_ = 0.0;
            setup_bins_();
        };


        void print() const {

            //check for assignment to self
            cout << "--- DENSITY MATRIX ------\n";
            cout << "x min:      " << x_min_ << "\n";
            cout << "x max:      " << x_max_ << "\n";
            cout << "x bins:     " << x_bins_ << "\n";
            cout << "x 1/bin size: " << x_bin_size_reciprocal_ << "\n\n";

            cout << "y min:      " << y_min_ << "\n";
            cout << "y max:      " << y_max_ << "\n";
            cout << "y bins:     " << y_bins_ << "\n";
            cout << "y 1/bin size: " << y_bin_size_reciprocal_ << "\n\n";

            cout << "z min:      " << z_min_ << "\n";
            cout << "z max:      " << z_max_ << "\n";
            cout << "z bins:     " << z_bins_ << "\n";
            cout << "z 1/bin size: " << z_bin_size_reciprocal_ << "\n\n";
        }


        void dump_bins() const {
            printf("\nDensity bins\n--------------------\n");
            for (size_t i = 0; i < counts_vector_.size(); ++i) {
                printf("%lu, %.4f\n", (unsigned long) i, counts_vector_[i]);
            }
            printf("--------------------\n\n");
        };

        inline PGVector3 getBinCenter(size_t bin) const {

            float x_bin_size = (x_max_ - x_min_) / float(x_bins_);
            float y_bin_size = (y_max_ - y_min_) / float(y_bins_);
            float z_bin_size = (z_max_ - z_min_) / float(z_bins_);

            size_t z_bin_num = static_cast<int>(floor(static_cast<float>(bin) / static_cast<float>(x_bins_ * y_bins_)));
            size_t y_bin_num = static_cast<int>(floor(
                    static_cast<float>(bin % (z_bins_ * y_bins_)) / static_cast<float>(x_bins_)));
            size_t x_bin_num = bin % x_bins_;
            PGVector3 center;
            center.x = x_min_ + (static_cast<float>(x_bin_num) + 0.5) * x_bin_size;
            center.y = y_min_ + (static_cast<float>(y_bin_num) + 0.5) * y_bin_size;
            center.z = z_min_ + (static_cast<float>(z_bin_num) + 0.5) * z_bin_size;

            return center;
        };

        vector<PGVector3> getBinCenters() const {
            vector<PGVector3> bin_centers(counts_vector_.size());
            for(size_t i=0; i<counts_vector_.size(); ++i){
                bin_centers[i] = getBinCenter(i);
            }
            return bin_centers;
        };


        string get3DDose(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax,  size_t nz, float zmin, float zmax) const {

            //Do not allow ranges outside of the reconstruction volume
            xmin = std::max(xmin, x_min_);
            xmax = std::min(xmax, x_max_);
            ymin = std::max(ymin, y_min_);
            ymax = std::min(ymax, y_max_);
            zmin = std::max(zmin, z_min_);
            zmax = std::min(zmax, z_max_);

            stringstream ss;
            ss.precision(4);

            ss << nx << " " << ny << " " << nz << endl;

            auto binedges_lambda = [&](float vmin, float vmax, size_t bins) {
                for (size_t i = 0; i <= bins; ++i) {
                    ss << vmin + float(i) * (vmax - vmin) / float(bins) << ",";
                };
                ss << endl;
            };
            binedges_lambda(xmin, xmax, nx);
            binedges_lambda(ymin, ymax, ny);
            binedges_lambda(zmin, zmax, nz);

            float step_size_x = (xmax - xmin)/float(nx);
            float step_size_y = (ymax - ymin)/float(ny);
            float step_size_z = (zmax - zmin)/float(nz);

            PGVector3 p;

            for (size_t i = 0; i < nz; ++i) {
                p.z = zmin + (float(i) + 0.5) * step_size_z;
                for (size_t j = 0; j < ny; ++j) {
                    p.y = ymin + (float(j) + 0.5) * step_size_y;
                    for (size_t k = 0; k < nx; ++k) {
                        p.x = xmin + (float(k) + 0.5) * step_size_x;
                        ss << getDensity(p) << ",";
                    }
                }
                ss << "\n";
            }
            return ss.str();
        };

        string get3DDose(size_t nx, size_t ny, size_t nz) const {

            return get3DDose(nx, x_min_, x_max_, ny, y_min_, y_max_, nz, z_min_, z_max_);
        };


        string get3DDose() const {
            return get3DDose(x_bins_, y_bins_, z_bins_);
        };

    private:
        float x0_min_;
        float x0_max_;
        size_t x0_bins_;
        float x_min_;
        float x_max_;
        size_t x_bins_;

        float x_bin_size_reciprocal_;

        float y0_min_;
        float y0_max_;
        size_t y0_bins_;
        float y_min_;
        float y_max_;
        size_t y_bins_;
        float y_bin_size_reciprocal_;

        float z0_min_;
        float z0_max_;
        size_t z0_bins_;
        float z_min_;
        float z_max_;
        size_t z_bins_;
        float z_bin_size_reciprocal_;

        size_t num_bins_;
        //we use a 1D valarray to represent the 3D dimensions
        // so that we can get a contiguous block of memory
        valarray<float> counts_vector_;
        bool poisson_binning_;
        ///bad_bin_ is returned when the bin_number is incorrect
        /// It should always be zero.
        long bad_bin_;

        pg_tools::Random rand_;

        void make_copy(const DensityMatrix &model) {

            //check for assignment to self
            if (this == &model) return;

            x_min_ = model.x_min_;
            x_max_ = model.x_max_;
            x_bins_ = model.x_bins_;
            if (0 == x_bins_) x_bins_ = 1;
            x_bin_size_reciprocal_ = static_cast<float>(x_bins_) / (x_max_ - x_min_);

            y_min_ = model.y_min_;
            y_max_ = model.y_max_;
            y_bins_ = model.y_bins_;
            if (0 == y_bins_) y_bins_ = 1;
            y_bin_size_reciprocal_ = static_cast<float>(y_bins_) / (y_max_ - y_min_);

            z_min_ = model.z_min_;
            z_max_ = model.z_max_;
            z_bins_ = model.z_bins_;
            if (0 == z_bins_) z_bins_ = 1;
            z_bin_size_reciprocal_ = static_cast<float>(z_bins_) / (z_max_ - z_min_);

            num_bins_ = model.num_bins_;
            bad_bin_ = model.bad_bin_;

            //need to makes sure vectors are same size before
            // copying.
            counts_vector_.resize(model.counts_vector_.size());
            counts_vector_ = model.counts_vector_;
        }

        std::shared_ptr<DensityEstimator> clone() const{

            std::shared_ptr<DensityMatrix> p(new DensityMatrix(*this));
            return p;
        }

        inline size_t getBinNumber_(const PGVector3 &pos) const {

            size_t bin_number = static_cast<int>((pos.x - x_min_) * x_bin_size_reciprocal_);
            bin_number += x_bins_ * static_cast<int>((pos.y - y_min_) * y_bin_size_reciprocal_);
            bin_number += y_bins_ * x_bins_ * std::min(static_cast<size_t>((pos.z - z_min_) * z_bin_size_reciprocal_), z_bins_ - 1);

            return bin_number;
        };

        void setup_bins_() {

            if (poisson_binning_) {
                x_bins_ = x0_bins_ + 1;
                float bin_size = (x_max_ - x_min_) / (x_bins_);
                float offset = rand_.Rndm() * bin_size;
                x_min_ = x0_min_ - offset;
                x_max_ = x_min_ + x_bins_ * bin_size;

                y_bins_ = y0_bins_ + 1;
                bin_size = (y0_max_ - y0_min_) / (y0_bins_);
                offset = rand_.Rndm() * bin_size;
                y_min_ = y0_min_ - offset;
                y_max_ = y_min_ + y_bins_ * bin_size;

                z_bins_ = z0_bins_ + 1;
                bin_size = (z0_max_ - z0_min_) / (z0_bins_);
                offset = rand_.Rndm() * bin_size;
                z_min_ = z0_min_ - offset;
                z_max_ = z_min_ + z_bins_ * bin_size;
            }
            num_bins_ = x_bins_ * y_bins_ * z_bins_;
            counts_vector_.resize(num_bins_);
            //counts_vector_ = 1.0;
            counts_vector_ = 0;
        }
    };

}
#endif //DENSITY_MATRIX
