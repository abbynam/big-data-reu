#ifndef DCA_KERNEL_DENSITY_ESTIMATOR_H_
#define DCA_KERNEL_DENSITY_ESTIMATOR_H_
#define _USE_MATH_DEFINES

//standard C++ includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <valarray>
#include <vector>
#include <stdexcept>

//ROOT includes
#include "TH3F.h"

//PromptGamma includes
#include "ComptonScatter.h"
#include "PGVector3.h"
#include "DensityEstimator.h"

using namespace std;
namespace prompt_gamma_reconstruction{
    
/*! \brief Estimate density using kernel density estimation rather than
 * histograms.
 * 
 * Gaussian kernels are used to estimate the density. This method produces a 
 * smoother results and is preferend for small number of events. However, 
 * it is order M*N*N where M is the number of iteration and M is the 
 * number of scatter cones -- it can be quite slow.
 * 
 * 
 * @author Dennis Mackin
 */
class DCAKernelDensityEstimator: public DensityEstimator{

  public:
    DCAKernelDensityEstimator(vector< ConicSection> *conic_sections,
                              const float bandwidth, const int num_threads,
                              double x_min, double x_max, int x_bins,
                              double y_min, double y_max, int y_bins,
                              double z_min, double z_max, int z_bins):
                              conic_sections_ptr_(conic_sections), bandwidth_(bandwidth), //num_threads_(num_threads),
                              x_min_(x_min), x_max_(x_max), x_bins_(x_bins),
                              y_min_(y_min), y_max_(y_max), y_bins_(y_bins),
                              z_min_(z_min), z_max_(z_max), z_bins_(z_bins)
    {


    };

    inline  void fill(const PGVector3 &pos, double weight) {
        /* DO NOTHING */
        return;
    };


    inline float operator[](const PGVector3 &pos) {
        return getDensity(pos);
    };

    inline void updateMatrix(const PGVector3 &oldPos, const PGVector3 &newPos, const float &weight) {
        cout<<"WARNING! DCAKernelDensityEstimator::updateMatrix is undefined . . ."<<endl;
        return;
    };

    inline void updateMatrix(const PGVector3 &oldPos, const PGVector3 &newPos) {
        cout<<"WARNING! DCAKernelDensityEstimator::updateMatrix is undefined . . ."<<endl;
        return;
    };


    //DCAKernelDensityEstimator &operator=(double rhs){
    DCAKernelDensityEstimator &operator=(double rhs){
      cout<<"WARNING!  operator=(double) undefined for DCAKernelDensityEstimator . . ."<<endl;
      return *this;
    }

    DCAKernelDensityEstimator &operator=(const DCAKernelDensityEstimator &rhs){
        make_copy(rhs);
        return *this;
    }

    inline bool is_in_volume(const PGVector3 &point) const{
        cout<<"WARNING! DCAKernelDensityEstimator::is_in_volume is undefined . . ."<<endl;
        return true;
    };

    PGVector3 getBinCenter(int bin) const{
        cout<<"WARNING! DCAKernelDensityEstimator::getBinCenter is undefined . . ."<<endl;
        return PGVector3();
    };


    void clear(){
        /* Nothing to clear for this estimator.*/
        return;
    };

   void print() const{
        cout<<"---- DCA KERNEL DENSITY ESTIMATOR ----\n";
    }


    // functions accepts event position vector and returns corresponding density value
    float getDensity(const PGVector3 &pos, const float bandwidth) const {

        const float sqrt_2PI = sqrt(2.0*M_PI);
        size_t num_cones = conic_sections_ptr_->size();
        float density = 0.0;
        float term = 0.0;
        float denominator = 1.0/(2.0*bandwidth*bandwidth);
        for(auto iter = conic_sections_ptr_->begin(); iter != conic_sections_ptr_->end(); ++iter){
            
            float dca = iter->getDistanceToPoint(pos);
            if(dca < 3.0 * bandwidth){
                term = exp(-dca*dca*denominator);
                density += term;
            }
        }

        float normalization_constant = 1.0/(sqrt_2PI*bandwidth* num_cones);
        return density*normalization_constant;
    };

    float getDensity(const PGVector3 &pos) const { return getDensity(pos, bandwidth_); };
    
    
    vector<float> getDensities(const vector<PGVector3> &points) const{ 
       return getDensities(points, bandwidth_);
    };
    
    vector<float> getDensities(const vector<PGVector3> &points, const float bandwidth) const { 
       
       auto numberOfCones = conic_sections_ptr_->size();
       auto numberOfPoints = points.size();
       
       const float sqrt_2PI = sqrt(2.0*M_PI);
       float normalization = 1.0/(sqrt_2PI*bandwidth* numberOfCones);
       float sigma3Inverse = 1.0/(3.0f*bandwidth);
       float denominator = -0.5/(bandwidth*bandwidth);
       
       vector<float> dcas(points.size());
       vector<float> densities(points.size(), 0.0);
       vector<PGVector3> nonconstPoints = points;
       
       size_t j = 0;
       for(j = 0; j < numberOfCones; ++j){
          auto distances = (*conic_sections_ptr_)[j].getDistanceToPoints(nonconstPoints);
          for(auto i=0u; i < numberOfPoints; ++i) {
             if(1.0f > distances[i]*sigma3Inverse) densities[i] += normalization* exp(dcas[i]*dcas[i] * denominator);
          }
       }
       return densities;
    };

    
    TH3F *getRootHist() const{
        TH3F *rootHist = new TH3F("densityHist","Density Estimation Matrix",
                                  x_bins_, x_min_, x_max_,
                                  y_bins_, y_min_, y_max_,
                                  z_bins_, z_min_, z_max_);
        PGVector3 p;
        //  #pragma omp parallel for num_threads(num_threads_), private(p,j,k)
        // ROOT is not thread safe. Need to make a thread safe hist and then copy to root.
        std::vector<PGVector3> points(x_bins_*y_bins_*z_bins_);
        
         for(int i=1; i<=z_bins_; ++i){
            p.z = z_min_ + (i - 0.5)*(z_max_ - z_min_)/float(z_bins_);
            for(int j=1; j<=y_bins_; ++j){
               p.y = y_min_ + (j - 0.5)*(y_max_ - y_min_)/float(y_bins_);
               for(int k=1; k<=x_bins_; ++k){
                    p.x = x_min_ + (k - 0.5)*(x_max_ - x_min_)/float(x_bins_);
                    points.push_back(p);
                    float density = getDensity(p);
                    int bin = rootHist->GetBin(i,j,k);
                    rootHist->SetBinContent(bin, density);
               }
           }
          cout<<"filled plane "<< i;
         }
        auto densities = getDensities(points);
         
       return rootHist;
    };


    TH2F getRootHist_xz(double y_depth=0.0) const{
        
        char buffer[100];
        sprintf(buffer, "density_hist_xz_y%.2f", y_depth);
        TH2F hist(buffer,"Density Estimation Matrix",
                                  z_bins_, z_min_, z_max_,
                                 x_bins_, x_min_, x_max_);
        PGVector3 p;
        vector<PGVector3> points(x_bins_*z_bins_);
        p.y = y_depth;
        double z_step = (z_max_ - z_min_)/float(z_bins_);
        double x_step = (x_max_ - x_min_)/float(x_bins_);
        for(int i = 1; i <= z_bins_; ++i){
            p.z = z_min_ + (i - 0.5)*z_step;
            for(int j=1; j<= x_bins_; ++j){
                p.x = x_min_ + (j - 0.5)*x_step;
                points.push_back(p);
            }
        }
        
        auto densities = getDensities(points);
        for(int i = 1; i <= z_bins_; ++i){
            for(int j=1; j<= x_bins_; ++j){
               auto bin = j + (i - 1)*x_bins_;
               hist.SetBinContent(bin, densities[bin]);
            }
        }        

        return hist;
    };
    
   TH2F *getRootHist_yz(double x_depth=0.0) const{
      TH2F * hist = new TH2F("densityHist_yz","Density Estimation Matrix",
                                  y_bins_, y_min_, y_max_,
                                  z_bins_, z_min_, z_max_);
        PGVector3 p;
        p.x = x_depth;
        for(int i = 0; i <= z_bins_; ++i){
            p.z = z_min_ + (i - 0.5)*(z_max_ - z_min_)/float(z_bins_);
            for(int j=0; j<= y_bins_; ++j){
               p.y = y_min_ + (j - 0.5)*(y_max_ - y_min_)/float(y_bins_);
               hist->SetBinContent(hist->GetBin(j,i), getDensity(p));
            }
        }            
        return hist;
   };
    
    
    TH2F *getRootHist_xy(double z_depth=0.0){
        
        TH2F * hist = new TH2F("densityHist_xy","Density Estimation Matrix",
                                    x_bins_, x_min_, x_max_,
                                    y_bins_, y_min_, y_max_);
        PGVector3 p;
        p.z = z_depth;
        for(int i = 0; i <= y_bins_; ++i){
            p.y = y_min_ + (i - 0.5)*(y_max_ - y_min_)/float(y_bins_);
            for(int j=0; j<= x_bins_; ++j){
                p.x = x_min_ + (j - 0.5)*(x_max_ - x_min_)/float(x_bins_);
                hist->SetBinContent(hist->GetBin(j,i), getDensity(p));
            }
        }      
        return hist;
    };
    
    
    
   TH1F *getRootHist_z(const int bins, const float range_min, const float range_max, const float bandwidth){
    
        char buffer[301];
        sprintf(buffer,"DCA_z_%.0f", bandwidth*10000);
        TH1F *hist = new TH1F(buffer,"Cental Axis DCA density", bins, range_min, range_max);

        PGVector3 p(0.0,0.0,0.0);
        vector<float> bin_vals(bins + 1,0.0);

        //#pragma omp parallel for num_threads(num_threads_), private(p, density)
        for(int iBin = 1; iBin <= bins; ++iBin){
            p.z = hist->GetBinCenter(iBin);
            hist->SetBinContent(iBin, getDensity(p, bandwidth));
        }
        
        return hist;
    };    

    TH1F *getRootHist_z(const int bins, const float range_min, const float range_max){
        return getRootHist_z(bins, range_min, range_max, bandwidth_);
    };

  private:

    vector<ConicSection> *conic_sections_ptr_;
    float bandwidth_;
    //int num_threads_;

    shared_ptr<DensityEstimator> density_estimator_ptr_;

    double x_min_;
    double x_max_;
    int x_bins_;

    double y_min_;
    double y_max_;
    int y_bins_;

    double z_min_;
    double z_max_;
    int z_bins_;
    

    void make_copy(const DCAKernelDensityEstimator &model){

        //check for assignment to self
        if(this == &model) return;

        bandwidth_ = model.bandwidth_;
        conic_sections_ptr_ = model.conic_sections_ptr_;

        x_min_ = model.x_min_;
        x_max_ = model.x_max_;
        x_bins_ = model.x_bins_;
        if(0 == x_bins_) x_bins_ = 1;

        y_min_ = model.y_min_;
        y_max_ = model.y_max_;
        y_bins_ = model.y_bins_;
        if(0 == y_bins_) y_bins_ = 1;

        z_min_ = model.z_min_;
        z_max_ = model.z_max_;
        z_bins_ = model.z_bins_;
        if(0 == z_bins_) z_bins_ = 1;

        cout<<"DCAKernelDensityEstimator::make_copy:"<<endl;
    }
};

}
#endif //DCA_KERNEL_DENSITY_ESTIMATOR_H_
