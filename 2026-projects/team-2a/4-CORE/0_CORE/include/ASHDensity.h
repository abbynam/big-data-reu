#ifndef ASH_DENSITY_H_
#define ASH_DENSITY_H_
#define _USE_MATH_DEFINES

//standard C++ includes
#include <cmath>
#include <cstdio>
#include <valarray>
#include <vector>
#include <stdexcept>
#include <omp.h>

//PromptGamma includes
#include "ComptonScatter.h"
#include "PGVector3.h"
#include "DensityEstimator.h"
#include "utilities/Random.h"


using namespace std;
namespace prompt_gamma_reconstruction {

    class ASHDensity : public DensityEstimator {

    public:
        ASHDensity(size_t num_shifts,
                float x_min, float x_max, size_t x_bins,
                float y_min, float y_max, size_t y_bins,
                float z_min, float z_max, size_t z_bins
        ): num_histograms_(num_shifts),
           x_min_(x_min), x_max_(x_max), x_bins_(x_bins),
           y_min_(y_min), y_max_(y_max), y_bins_(y_bins),
           z_min_(z_min), z_max_(z_max), z_bins_(z_bins)
        {
            x_bin_size_reciprocal_ = static_cast<float>(x_bins_) / (x_max_ - x_min_);
            y_bin_size_reciprocal_ = static_cast<float>(y_bins_) / (y_max_ - y_min_);
            z_bin_size_reciprocal_ = static_cast<float>(z_bins_) / (z_max_ - z_min_);
            num_histograms_reciprocal_ = 1.0/static_cast<float>(num_histograms_);

            x_bins_++;
            y_bins_++;
            z_bins_++;

            counts_vector_.resize(x_bins_*y_bins_*z_bins_*num_histograms_);
            build_shifted_volumes_();
        };

        ASHDensity() :
                x_min_(-1.0), x_max_(1.0), x_bins_(1),
                y_min_(-1.0), y_max_(1.0), y_bins_(1),
                z_min_(-1.0), z_max_(1.0), z_bins_(1), num_histograms_(1)
        {
            cout << "ASHDensity::ASHDensity(empty constructor):" << counts_vector_.size() << endl;
            build_shifted_volumes_();
        };

        ~ASHDensity() {
            cout << "destroying ASHDensity . . ." << endl;
        }

        // functions accepts event position vector and returns cooresponding density value
        inline float getDensity(const PGVector3 &pos0) const {

            size_t bin = 0;
            float density = 0.0;
            PGVector3 pos = pos0;

            //jitter values
//            pos.x += -(pos.x <= x_min_)*1E-6 + (pos.x >= x_max_)*1E-6;
//            pos.y += -(pos.y <= y_min_)*1E-6 + (pos.y >= y_max_)*1E-6;
//            pos.z += -(pos.z <= z_min_)*1E-6 + (pos.z >= z_max_)*1E-6;

            for(size_t k=0; k < num_histograms_; ++k){

                bin = static_cast<size_t>((pos.x - x_mins_[k]) * x_bin_size_reciprocal_);
                bin += x_bins_ * static_cast<size_t>((pos.y - y_mins_[k]) * y_bin_size_reciprocal_);
                bin += y_bins_ * x_bins_ * static_cast<size_t>((pos.z - z_mins_[k]) * z_bin_size_reciprocal_);
                bin += k*z_bins_*y_bins_*x_bins_;

                if (counts_vector_.size() <= bin) {
                    stringstream ss;
                    ss << "ERROR in ASHDensity::getDensity(" << pos.x << "," << pos.y << "," << pos.z <<
                    ") returned bin " << bin << "\n" << counts_vector_.size() - 1 << " is maximum bin.\n";
                    ss << "x range(" << x_min_ << ", " << x_max_ <<")\n";
                    ss << "y range(" << y_min_ << ", " << y_max_ <<")\n";
                    ss << "z range(" << z_min_ << ", " << z_max_ <<")\n";
                    ss << "bins(" << x_bins_ <<", " << y_bins_ << ", " << z_bins_ <<")" <<endl;
                    ss << "bin_size_reciprocals (" << x_bin_size_reciprocal_ <<", " << y_bin_size_reciprocal_ << ", " << z_bin_size_reciprocal_ <<")" <<endl;
                    cout << ss.str() << endl;
                    throw runtime_error(ss.str());
                }

                density += this->counts_vector_[bin];
            }
            density *= num_histograms_reciprocal_;

            return density;
        };


        inline void fill(const PGVector3 &pos, float weight) {

            size_t bin = 0;
//            float density = 0.0;
            for(size_t k=0; k < num_histograms_; ++k) {
                bin = static_cast<size_t>((pos.x - x_mins_[k]) * x_bin_size_reciprocal_);
                bin += x_bins_ * static_cast<size_t>((pos.y - y_mins_[k]) * y_bin_size_reciprocal_);
                bin += y_bins_ * x_bins_ * static_cast<size_t>((pos.z - z_mins_[k]) * z_bin_size_reciprocal_);
                bin += k * z_bins_ * y_bins_ * x_bins_;

                //DSM 96% of the octane runtime is in this one line
                #pragma omp atomic update
                counts_vector_[bin] += weight;
            }
        };


        inline void updateMatrix(const PGVector3 &oldPos, const PGVector3 &newPos) {
            updateMatrix(oldPos, newPos, 1.0);
        };


        inline void updateMatrix(const PGVector3 &oldPos, const PGVector3 &newPos, const float &weight) {
            fill(oldPos, -1.0*weight);
            fill(newPos, weight);
        };


        /*! \brief Shallow copy of the density matrix -- not a deep copy of the ASH object  */
        ASHDensity &operator=(float rhs) {
            counts_vector_ = rhs;
            return *this;
        }

        std::shared_ptr<DensityEstimator> clone() const{

            cout << "Cloning ASHDensity . . ." << endl;
            std::shared_ptr<ASHDensity> p(new ASHDensity());
            cout << "Copying counts_vector . . ." << endl;
            p->counts_vector_ = this->counts_vector_;
            cout << "New counts_vector_ length: "<< p->counts_vector_.size() << " . . ." << endl;
            p->x_mins_ = this->x_mins_;
            p->x_min_ = this->x_min_;
            p->x_max_ = this->x_max_;
            p->x_bins_ = this->x_bins_;

            p->x_bin_size_reciprocal_ = this->x_bin_size_reciprocal_;

            p->y_mins_ = this->y_mins_;
            p->y_min_ = this->y_min_;
            p->y_max_ = this->y_max_;
            p->y_bins_ = this->y_bins_;
            p->y_bin_size_reciprocal_ = this->y_bin_size_reciprocal_;

            p->z_mins_ = this->z_mins_;
            p->z_min_ = this->z_min_;
            p->z_max_ = this->z_max_;
            p->z_bins_ = this->z_bins_;
            p->z_bin_size_reciprocal_ = this->z_bin_size_reciprocal_;

            p->num_histograms_ = this->num_histograms_;
            p->num_histograms_reciprocal_ = this->num_histograms_reciprocal_;

            return p;
        };

        /*! \brief Shallow copy of the density matrix -- not a deep copy of the ASH object  */
        ASHDensity &operator=(const ASHDensity &rhs) {
            cout << "Copying ASHDensity . . ." << endl;
            counts_vector_ = rhs.counts_vector_;
            x_mins_ = rhs.x_mins_;
            x_min_ = rhs.x_min_;
            x_max_ = rhs.x_max_;
            x_bins_ = rhs.x_bins_;

            x_bin_size_reciprocal_ = rhs.x_bin_size_reciprocal_;

            y_mins_ = rhs.y_mins_;
            y_min_ = rhs.y_min_;
            y_max_ = rhs.y_max_;
            y_bins_ = rhs.y_bins_;
            y_bin_size_reciprocal_ = rhs.y_bin_size_reciprocal_;

            z_mins_ = rhs.z_mins_;
            z_min_ = rhs.z_min_;
            z_max_ = rhs.z_max_;
            z_bins_ = rhs.z_bins_;
            z_bin_size_reciprocal_ = rhs.z_bin_size_reciprocal_;

            num_histograms_ = rhs.num_histograms_;
            num_histograms_reciprocal_ = rhs.num_histograms_reciprocal_;

            return *this;
        }


        ASHDensity &operator+=(const ASHDensity &rhs) {
            if (counts_vector_.size() != rhs.counts_vector_.size()) {
                stringstream ss;
                ss << "ERROR in ASHDensity &operator+= :\n";
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

        vector<float> getDensities(const vector<PGVector3> &points, vector<float> &densities) const{
            size_t i = 0;
            size_t num_points = points.size();
            densities.resize(num_points);

            #pragma omp parallel for
            for(i = 0; i < num_points; ++i){
                densities[i] = getDensity(points[i]);
            }

            return densities;
        }

        vector<float> getDensities(const vector<PGVector3> &points) const{
            vector<float> densities(points.size());
            for(size_t i = 0; i < points.size(); ++i){
                densities[i] = getDensity(points[i]);
            }
            return densities;
        }

        valarray<float> getBins() const{ return counts_vector_; };
        void setBins(valarray<float> bins) {
            if (counts_vector_.size() != bins.size()) {
                stringstream ss;
                ss << "ERROR: ASHDensity::setBins failed.\n";
                ss << "matrix size(" << counts_vector_.size() << ") != right hand side(" << bins.size() << ").\n\n";
                cout << ss.str() << endl;
                throw runtime_error(ss.str());
            }
            counts_vector_ = bins;
        };


        void clear() {
            counts_vector_ = 0.0;
        };


        void print() const{

            cout<<"----- ASH DENSITY MATRIX ------\n";
            cout<<"----- "<< num_histograms_ <<" histograms ------\n";
            cout<<"x min:      "<< x_min_ <<"\n";
            cout<<"x max:      "<< x_max_ <<"\n";
            cout<<"x bins:     "<< x_bins_ <<"\n\n";

            cout<<"y min:      "<< y_min_ <<"\n";
            cout<<"y max:      "<< y_max_ <<"\n";
            cout<<"y bins:     "<< y_bins_ <<"\n\n";

            cout<<"z min:      "<< z_min_ <<"\n";
            cout<<"z max:      "<< z_max_ <<"\n";
            cout<<"z bins:     "<< z_bins_ <<"\n\n";
        }


        /*! \brief Finds the center of bin
         */
        inline PGVector3 getBinCenter(size_t bin) const{

            float x_bin_size = (x_max_ - x_min_)/float(x_bins_);
            float y_bin_size = (y_max_ - y_min_)/float(y_bins_);
            float z_bin_size = (z_max_ - z_min_)/float(z_bins_);

            size_t z_bin_num = static_cast<int>(floor(static_cast<float>(bin)/ static_cast<float>(x_bins_ * y_bins_)));
            size_t y_bin_num = static_cast<int>(floor(static_cast<float>(bin % (z_bins_*y_bins_))/ static_cast<float>(x_bins_)));
            size_t x_bin_num = bin % x_bins_;
            PGVector3 center;
            center.x = x_min_ + (static_cast<float>(x_bin_num) + 0.5) * x_bin_size;
            center.y = y_min_ + (static_cast<float>(y_bin_num) + 0.5) * y_bin_size;
            center.z = z_min_ + (static_cast<float>(z_bin_num) + 0.5) * z_bin_size;

            return center;
        };


        string get3DDose(size_t nx, float xmin, float xmax, size_t ny, float ymin, float ymax,  size_t nz, float zmin, float zmax) const{

            cout << "ASHDensity::get3DDose: Producing the dose output string . . ." << endl;
            cout << "Output range: (" << xmin <<"," << xmax << "," << nx << ") (" << ymin << "," << ymax << "," << ny <<") (" << zmin << "," << zmax << ","  << nz << ")" << endl;
            stringstream ss;
            ss.precision(7);

            //Do not allow values outside of the reconstruction volume
            xmin = std::max(xmin, x_min_);
            xmax = std::min(xmax, x_max_);
            ymin = std::max(ymin, y_min_);
            ymax = std::min(ymax, y_max_);
            zmin = std::max(zmin, z_min_);
            zmax = std::min(zmax, z_max_);

            ss<<nx<<" "<<ny<<" "<<nz<<endl;

            auto binedges_lambda = [&](float vmin, float vmax, size_t bins){
                for(size_t i = 0; i <= bins; ++i){
                    ss << vmin + float(i)*(vmax - vmin)/float(bins) << " ";
                };
                ss<<endl;
            };
            binedges_lambda(xmin, xmax, nx);
            binedges_lambda(ymin, ymax, ny);
            binedges_lambda(zmin, zmax, nz);

            PGVector3 p;
            vector<PGVector3> points(ny * nx);
            //vector<float> densities(ny * nx);
            for(size_t i=0; i< nz; ++i){
                p.z = zmin + (float(i) + 0.5)*(zmax - zmin)/float(nz);

                for(size_t j=0; j<ny; ++j){
                    p.y = ymin + (float(j) + 0.5)*(ymax - ymin)/float(ny);
                    for(size_t k=0; k<nx; ++k){
                        auto step_size = (xmax - xmin)/float(nx);
                        p.x = xmin + (float(k) + 0.5)*step_size;
                        points[j*nx + k].x = p.x;
                        points[j*nx + k].y = p.y;
                        points[j*nx + k].z = p.z;
                    }
                }
                auto densities = getDensities(points);
//                auto densities = getDensities(points);
                for_each(densities.begin(), densities.end(), [&](float density){ss << density << ",";});
                if(i%10 == 0) cout << "ASHDensity::get3DDose(): Calculated slice " << i << " . . ." << endl;
                ss<<"\n";
            }
            return ss.str();
        } ;

        string get3DDose() const {
            auto scalar = log(M_E + float(num_histograms_));
            return get3DDose(min(int(x_bins_*scalar), 300), min(int(y_bins_*scalar),300), min(int(z_bins_*scalar),300));
        }


        string get3DDose(size_t nx, size_t ny, size_t nz) const{
            return get3DDose(x_min_, x_max_, nx, y_min_, y_max_, ny, z_min_, z_max_, nz);
        } ;

    private:
        size_t num_histograms_;
        valarray<float> x_mins_;
        float x_min_;
        float x_max_;
        size_t x_bins_;

        float x_bin_size_reciprocal_;

        valarray<float> y_mins_;
        float y_min_;
        float y_max_;
        size_t y_bins_;
        float y_bin_size_reciprocal_;

        valarray<float> z_mins_;
        float z_min_;
        float z_max_;
        size_t z_bins_;
        float z_bin_size_reciprocal_;

        float num_histograms_reciprocal_;

        // counts_vector will be x_bins*y_bins*zbins**num_histograms in length
        valarray<float> counts_vector_;
        //DSM  too slow ==>  vector<omp_lock_t> bin_locks_;


        /*! \brief Builds the shifted histograms for averaging
           *
           * Methods determines the positions of the bin edges for each
           * histogram. The shifts are determined individually for each
           * x,y,z dimension and then randomized. This method requires
           * N histograms rather than N^3 if shifts were were applied one
           * dimension at a time.
           */
        void build_shifted_volumes_(){
            float x_bin_size = 1.0/x_bin_size_reciprocal_;
            float y_bin_size = 1.0/y_bin_size_reciprocal_;
            float z_bin_size = 1.0/z_bin_size_reciprocal_;
            x_mins_.resize(num_histograms_);
            y_mins_.resize(num_histograms_);
            z_mins_.resize(num_histograms_);


            for(size_t k=0; k < num_histograms_; ++k){
                x_mins_[k] = x_min_ - x_bin_size + k*(x_bin_size/static_cast<float>(num_histograms_));
                y_mins_[k] = y_min_ - y_bin_size + k*(y_bin_size/static_cast<float>(num_histograms_));
                z_mins_[k] = z_min_ - z_bin_size + k*(z_bin_size/static_cast<float>(num_histograms_));
            }

            std::random_device rd;
            std::mt19937 g(rd());
            g.seed(123);
            std::shuffle(&x_mins_[0], &x_mins_[num_histograms_], g);
            std::shuffle(&y_mins_[0], &y_mins_[num_histograms_], g);
            std::shuffle(&z_mins_[0], &z_mins_[num_histograms_], g);
        };
    };

}
#endif //ASH_DENSITY
