cmake_minimum_required(VERSION 3.5)
project(coreGPU LANGUAGES CXX CUDA)


#set(BASEDIR "/Users/edraeger/Development/PromptGamma/CORE")
set(BASEDIR ".")
include_directories(${BASEDIR}/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -O3 -Wall -Wno-reorder -fno-omit-frame-pointer -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -O3 -fno-omit-frame-pointer")
#-Wsign-compare
#-Wunused-but-set-variable
#-Wunused-variable
#-Wsequence-point
#-Wformat=
#-Wall

set(SOURCE_FILES
        include/AbstractImageAlgorithmBuilder.h
        include/AbstractResultsFormatter.h
        include/ASHDensity.h
        include/ComptonScatter.h
        include/Co60Algorithm.h
        include/Co60ResultsFormatter.h
        include/ConicSection.h
        include/CSVEventsLoader.h
        include/DCALineAlgorithm.h
        include/DCALineAlgorithmBuilder.h
        include/DCALineEventsLoader.h
        include/DensityEstimator.h
        include/DensityEstimatorFactory.h
        include/DensityMatrix.h
        include/DetectorEffectsEventsLoader.h
        include/DoubleScatter.h
        include/EventFilters.h
        include/EventDataAlgorithm.h
        include/EventDataAlgorithmBuilder.h
        include/EventsLoaderFactory.h
        include/ImageAlgorithm.h
        #    include/ImageAlgorithmFactory.h
        include/KernelConic.h
        include/OctaneAlgorithm.h
        include/OctaneEMAlgorithm.h
        include/KEMAlgorithm.h
        include/OctaneAlgorithmBuilder.h
        include/KEMAlgorithmBuilder.h
        #    include/PGSimulatedEventsLoader.h
        include/PGVector3.h
        include/PhantomVolume.h
        include/PhantomVolumeBuilder.h
        include/PointAlgorithm.h
        include/PointAlgorithmBuilder.h
        #    include/RandomPointOnCircleSingleton.h
        include/RandomSingleton.h
        include/RandomSqrtSingleton.h
        include/ReconstructionEllipse.h
        include/ReconstructionParabola.h
        include/ReconstructorBuilder.h
        include/ReconstructorTemplate.h
        include/ResultsFormatterFactory.h
        include/RunTimeParameters.h
        include/Scatter.h
        include/SimpleResultsFormatter.h
        #    include/SimulatedCCEventsLoader.h
        include/SOEAlgorithm.h
        include/SOEAlgorithmBuilder.h
        include/StringMap.h
        include/StringUtils.h
        include/TripleScatter.h
        include/SBPAlgorithmBuilder.h
        include/SBPAlgorithm.h
        include/VectorAlgorithm.h
        src/Co60Algorithm.cc
        src/Co60ResultsFormatter.cc
        src/ComptonScatter.cc
        src/ConicSection.cc
        src/CSVEventsLoader.cc
        src/DCALineAlgorithm.cc
#        src/DCALineAlgorithmBuilder.cc
        src/DCALineEventsLoader.cc
        src/DoubleScatter.cc
        src/EventFilters.cc
        src/EventDataAlgorithm.cc
        src/EventsLoaderFactory.cc
        src/FileUtils.cc
        src/ImageAlgorithm.cc
        src/KernelConic.cc
        src/OctaneAlgorithm.cc
        src/OctaneEMAlgorithm.cc
        src/KEMAlgorithm.cc
        #    src/PGSimulatedEventsLoader.cc
        src/PointAlgorithm.cc
        src/RandomSingleton.cc
        src/RandomSqrtSingleton.cc
        #    src/ReconstructionEllipse.cc
        #    src/ReconstructionParabola.cc
        src/ReconstructorTemplate.cc
        src/ResultsFormatterFactory.cc

        src/RunTimeParameters.cc
        src/SimpleResultsFormatter.cc
        #    src/SimulatedCCEventsLoader.cc
        src/SOEAlgorithm.cc
        src/StringMap.cc
        src/StringUtils.cc
        src/TripleScatter.cc
        src/VectorAlgorithm.cc
        src/SBPAlgorithm.cc
        cuda/cuda_functions.cu

        include/OriginConesSoA.h

        )

add_library(core_objects STATIC ${SOURCE_FILES})
target_compile_features(core_objects PUBLIC cxx_std_11 )
set_target_properties(core_objects PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable(coreGPU src/runRecon.cu)
set_property(TARGET coreGPU PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(coreGPU PRIVATE core_objects)