a#include <stdio.h>
#include <math.h>
#include <cstring>
#include "cuda_fp16.h"
#include "cuda_functions.h"
#include "OriginConesSoA.h"
#include <cuda_runtime.h>

#include "../include/PGVector3.h"
#include "../include/OriginConesSoA.h"

using namespace std;
using namespace prompt_gamma_reconstruction;


 __device__ float calc_dca(float v1x, float v1y, float v1z,
                             float axis_x, float axis_y, float axis_z, float cos_angle){



    float v1mag = __fsqrt_rd(v1x*v1x + v1y*v1y + v1z*v1z);
    //float axis_mag = __fsqrt_rd(axis_x*axis_x + axis_y*axis_y + axis_z*axis_z);

    //float cos_theta = (v1x*axis_x + v1y*axis_y + v1z*axis_z)/(axis_mag*v1mag);
    float cos_theta = (v1x*axis_x + v1y*axis_y + v1z*axis_z)/(1.0*v1mag);

    cos_theta += (cos_theta < 0.0f)*1.0E-4f - (cos_theta > 0.0f)*1.0E-4f; //floating point error can make abs(cos_theta)>1

    //float y = v1mag * cos_theta;
    float dca = fabs((__fsqrt_rd(v1mag*v1mag - v1mag * cos_theta*v1mag * cos_theta) - v1mag * cos_theta*__fsqrt_rd(1.0f/(cos_angle*cos_angle) - 1.0f)) * cos_angle);
//     float dca = fabs((__fsqrt_rd(v1mag*v1mag - y*y) - y*__fsqrt_rd(1.0f/(cos_angle*cos_angle) - 1.0f)) * cos_angle);

     //Including this if make the code about 33% faster. I have no explanation.
    if(!(0.0f <= dca)) {
        printf("angle=%.2f, apex=(%.1f, %.1f, %.1f), axis=(%.3f, %.3f, %.3f), p=(%.3f, %.3f, %.3f)\n",
               cos_angle, v1x, v1y, v1z, axis_x, axis_y, axis_z, v1x, v1y, v1z);
//        printf("costheta %.3f, y %.2f, v1mag %.2f, axis_mag %.2f\n", cos_theta, y, v1mag, axis_mag);

        printf("is this still fast?\n");
    }

     return dca;
}


__device__ void calc_center(unsigned int bin, float* roi, float *p){

    if(bin > roi[2] * roi[5] * roi[8]) {
        printf("BAD BIN %d, %.3f, %.3f, %.3f\n", bin, roi[2], roi[5], roi[8]);
        bin = roi[2] * roi[5] * roi[8];
    }
    //assert(bin < roi[2] * roi[5] * roi[8]);

    unsigned int zbin = bin/(roi[2]*roi[5]);
    unsigned int ybin = (bin - zbin*(roi[2]*roi[5]))/roi[2];
    unsigned int xbin = bin - zbin*(roi[2]*roi[5]) - ybin* roi[2];

    float x_step = (roi[1] - roi[0])/roi[2];
    float y_step = (roi[4] - roi[3])/roi[5];
    float z_step = (roi[7] - roi[6])/roi[8];

    p[0] = roi[0] + (xbin + 0.5)*x_step;
    p[1] = roi[3] + (ybin + 0.5)*y_step;
    p[2] = roi[6] + (zbin + 0.5)*z_step;

    return;
}

__global__
void getDCA(float *matrix, float *volume_grid_info, int num_bins, float *cone_data, int num_cones, float bandwidth){
/*
    matrix: 3D array of voxel intensities
    volume_grid_info: array of 9 floats with range and number of bins for x,y,z dimensions
    num_bins: total number of voxels
    cone_data: 7 arrays, one after another (cone apex x,y,z, cone axis x,y,z, cosine of cone opening angle), 7*num_cones floats
    num_cones: number of cones used in the reconstructions
    bandwidth: smoothing factor for kernel density smoothing

*/
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    float tmp = 0.0f;
    float v1mag = 0.0f;

    float cos_theta = 0.0f;
    float p[3];

    //calc_center(12, volume_grid_info, p);
    if(i < num_bins) {
        calc_center(i, volume_grid_info, p);
        matrix[i] = 0.0f;

        //with unroll 46 sec.without 39 -- don't use unroll
        for(size_t j=0; j< num_cones; ++j){

            //calculate x coordinate
            //tmp = volume_grid_info[0 * num_bins + i] - cone_data[num_cones * 0 + j]; //15 sec to here
            tmp = p[0] - cone_data[num_cones * 0 + j]; //15 sec to here
            v1mag = tmp * tmp;
            cos_theta = tmp * cone_data[num_cones * 3 + j];

            //calculate y coordinate
            //tmp = volume_grid_info[1 * num_bins + i] - cone_data[num_cones * 1 + j];
            tmp = p[1] - cone_data[num_cones * 1 + j];
            v1mag += tmp * tmp;
            cos_theta += tmp * cone_data[num_cones * 4 + j]; //11 sec to here

            //calculate z coordinate
            //tmp = volume_grid_info[2 * num_bins + i] - cone_data[num_cones * 2 + j];
            tmp = p[2] - cone_data[num_cones * 2 + j];
            v1mag += tmp * tmp;
            cos_theta += tmp * cone_data[num_cones * 5 + j];

            v1mag = __fsqrt_rd(v1mag); //18 sec to here
            cos_theta /= v1mag;
            cos_theta += (cos_theta < 0.0f)*1.0E-4f - (cos_theta > 0.0f)*1.0E-4f; //floating point error can make abs(cos_theta)>1
            tmp = cone_data[num_cones * 6 + j];   //20 sec to here

            tmp = fabs((__fsqrt_rd(v1mag * v1mag - v1mag * cos_theta * v1mag * cos_theta) - v1mag * cos_theta * __fsqrt_rd(1.0f / (tmp * tmp) - 1.0f)) * tmp); //46 sec to here

            tmp = tmp / (bandwidth + 1.0E-5f);
            tmp = (tmp < 1.0f) * 0.75f * (1.0f - tmp * tmp);

            //DSM - this if with print statement somehow reduces the number of registers and speeds up the code
            if(!(0.0f <= tmp)) {
//                printf("angle=%.2f, apex=(%.1f, %.1f, %.1f), axis=(%.3f, %.3f, %.3f), axis=(%.3f, %.3f, %.3f)\n",
//                       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, cos_theta, v1mag, x);
                printf("%.1f, %.1f,%.1f,%.1f,%.1f, %.1f, %.1f, %.1f, %.1f, %0.1f\n", 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, cos_theta, v1mag,
                       tmp, 0.0f);
//                printf("1234567890123456789012345678901234567890 %200.1f\n",x);
//                printf("is this still fast?\n");
//                printf("angle=%.2f, apex=(%.1f, %.1f, %.1f), axis=(%.3f, %.3f, %.3f), axis=(%.3f, %.3f, %.3f)\n",0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, cos_theta, v1mag, x);
            }
          //  voxel_total += x;
            matrix[i] += tmp;
        }
        //if(i%100000 == 0) printf("%d %.3f\n", i, matrix[i]);
        //matrix[i] = voxel_total;
    }
};



void populate_density_matrix_cuda(vector<float> &density_matrix, const float *volume_grid_info, const OriginConesSoA &cone_soa, const float bandwidth){
    //Get single vector<float> of the apex, direction, and opening angle of the cones for copying to GPU
    vector<float> soa = cone_soa.get_arranged_memory();

    size_t N = cone_soa.getLength();
    size_t M = density_matrix.size();
    size_t N_data_items = 7; //3 apex coordinates, 3 axis coordinates, 1 cone angle
    size_t volume_grid_info_size = 9 * sizeof(float);

    cout <<"soa: "<<soa[0 + 3] << " " << soa[N + 3] <<" "<< soa[2*N + 3] <<" "<<soa[3*N + 3] <<" "<< soa[N_data_items*N - 1 + 3] <<endl;

    float *d_volume_grid_info = 0;

    float *d_density_matrix;
    assert(soa.size() % 7 == 0);

    float *d_cone_data;
    unsigned cone_data_size = N * 7 * sizeof(float);

    cudaMalloc((float **)&d_cone_data, cone_data_size);
    cout << "CUDA copy soa: " << cudaGetErrorString(cudaMemcpy(d_cone_data, &soa[0], cone_data_size, cudaMemcpyHostToDevice)) << endl;

    cudaMalloc((float **)&d_volume_grid_info, volume_grid_info_size);
    cudaMemcpy(d_volume_grid_info, &volume_grid_info[0], volume_grid_info_size, cudaMemcpyHostToDevice);

    cout << "CUDA allocate memory: " << cudaGetErrorString(cudaMalloc((float **)&d_density_matrix, M * sizeof(float))) << endl;
    cout << "CUDA copy to gpu result: " << cudaGetErrorString(cudaMemcpy(d_density_matrix, &density_matrix[0], M * sizeof(float), cudaMemcpyHostToDevice)) << endl;

    size_t threads = 256; //32 and 256 give similar results
    getDCA<<<(M + threads - 1)/threads, threads>>>(d_density_matrix, d_volume_grid_info, M, d_cone_data, N, bandwidth);

    cout << "CUDA copy back result: " << cudaGetErrorString(cudaMemcpy(&density_matrix[0], d_density_matrix,
                                                                       M * sizeof(float), cudaMemcpyDeviceToHost)) << endl;

    cudaFree(d_cone_data);
    cudaFree(d_volume_grid_info);
};


float getDensity(const PGVector3 p, const vector<OriginCone> &cones, float bandwidth_inv){
    float density = 0.0;
    float u = 0.0;

    //estimate density using Epanechinikov kernel
    for(size_t i=0; i < cones.size(); ++i){
        u = cones[i].get_DCA(p)*bandwidth_inv;
        density += (u < 1.0f)*0.75f*(1.0f - u*u);
    }
    return density;
};

