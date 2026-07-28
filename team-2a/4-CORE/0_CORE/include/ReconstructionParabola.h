#ifndef RECONSTRUCTION_PARABOLA_H_
#define RECONSTRUCTION_PARABOLA_H_

//standard C++/C libraries
#include <string>
#include <vector>

//local include files
#include "ComptonScatter.h"
#include "ConicSection.h"


namespace prompt_gamma_reconstruction{


  class ReconstructionParabola: public ConicSection {

  private:

  public:
    ReconstructionParabola(const ComptonScatter & comptonScatter, shared_ptr<const PhantomVolume> phantomVolume, const size_t seed);
    ~ReconstructionParabola(){ /*std::cout<<"deleting ReconstructionParabola . . ."<<std::endl; */};

    bool doesConicIntersectPhantom();
    inline size_t getRandomPointInPhantom0(PGVector3 &randPoint, size_t num_random_tries=0);

    void changeFramePhantomToConicSection();
    void changeFrameConicSectionToCone();
    void changeFrameConeToConicSection();
    void changeFrameConicSectionToPhantom();

    void getLineIntercepts_(double a, double h, double y,
                            std::pair<double, double> &point1, std::pair<double, double> &point2,
                            std::vector<PGVector3> &intersection_x_values) const;
    void draw_parabola(double a, double h, Rectangle rectanglePoints, size_t numberIntersections, double y);
    
    void testTransformations();
    void print();
  };

}//end of namespace prompt_gamma_reconstruction
#endif //RECONSTRUCTION_PARABOLA_H_
