//
// Created by dsmackin on 1/9/16.
//

#ifndef TRUST11_DENSITYESTIMATORFACTORY_H
#define TRUST11_DENSITYESTIMATORFACTORY_H


namespace prompt_gamma_reconstruction{
    class DensityEstimatorFactory{

    public:
        static shared_ptr<DensityEstimator> createDensityEstimator(const pg_tools::RunTimeParameters &params, shared_ptr<const PhantomVolume> pv_ptr){
            size_t density_estimator_type = params.get_int("DENSITY_ESTIMATOR_TYPE");

            if(2 == density_estimator_type){//using averaged shifted histograms
                cout<<"Using Averaged Shifted Histograms for density estimation . . ."<<endl;
                size_t number_of_shifts = params.get_int("NUMBER_OF_SHIFTS");
                return shared_ptr<DensityEstimator>(
                        new ASHDensity(number_of_shifts,
                                       pv_ptr->x_min, pv_ptr->x_max, pv_ptr->x_bins_,
                                       pv_ptr->y_min, pv_ptr->y_max, pv_ptr->y_bins_,
                                       pv_ptr->z_min, pv_ptr->z_max, pv_ptr->z_bins_ ));
            }else{
                cout<<"Using a standard 3D histogram for density estimation . . ."<<endl;
                return  shared_ptr<DensityMatrix>(
                        new DensityMatrix(
                                pv_ptr->x_min, pv_ptr->x_max, pv_ptr->x_bins_,
                                pv_ptr->y_min, pv_ptr->y_max, pv_ptr->y_bins_,
                                pv_ptr->z_min, pv_ptr->z_max, pv_ptr->z_bins_ ) );
            }
        };
        static shared_ptr<DensityMatrix> createDensityMatrix(shared_ptr<const PhantomVolume> pv_ptr, const size_t nx, const size_t ny, const size_t nz){

            return  shared_ptr<DensityMatrix>( new DensityMatrix(
                            pv_ptr->x_min, pv_ptr->x_max, nx,
                            pv_ptr->y_min, pv_ptr->y_max, ny,
                            pv_ptr->z_min, pv_ptr->z_max, nz)
            );
        };
    };
}

#endif //TRUST11_DENSITYESTIMATORFACTORY_H
