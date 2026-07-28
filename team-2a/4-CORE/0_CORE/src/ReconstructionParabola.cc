#include "ReconstructionParabola.h"


using namespace std;
using namespace prompt_gamma_reconstruction;



ReconstructionParabola::ReconstructionParabola(const ComptonScatter & comptonScatter, shared_ptr<const PhantomVolume> phantomVolume, const size_t seed):
                      ConicSection(comptonScatter, phantomVolume, seed){


  //testTransformations();
};


void ReconstructionParabola::getLineIntercepts_(double a, double vertex_x, double y,
                            std::pair<double, double> &point1, std::pair<double, double> &point2,
                            std::vector<PGVector3> &intersection_x_values) const{
    /// Solutions take the following form ///
    /******************************************************************************
    x -> h + 2 a slope + 2 Sqrt[a (intercept + slope (h + a slope))]

    x -> h + 2 a slope - 2 Sqrt[a (intercept + slope (h + a slope))]
    ******************************************************************************/
    /// If the parabola and the line intersect then the radical term will be real.
    /// There is no intersection if the term under the radical is negative.
    double slope_denominator = point1.first - point2.first;
    double slope_numerator = point1.second - point2.second;
    double x_neg, x_pos, z_neg, z_pos;

    if(0.0 == slope_denominator){ //then x is constant
        if(0.0 > point1.first - vertex_x){//no real intercepts
          return;
        }
        z_pos = sqrt(4.0*a*(point1.first - vertex_x));
        z_neg = -z_pos;
        x_neg = point1.first;
        x_pos = point1.first;
    }else if(0.0 == slope_numerator){ //then z is constant,
                      //exactly one intercept but not necessarily in phantom
    x_pos = (0.25/a)*point1.second*point1.second + vertex_x;
    z_pos = point1.second;
    x_neg = point1.second + 1.0; //move 2nd x value off the line
    z_neg = point1.second;
    }else{// 0 or exactly two intersections
        double slope = slope_numerator/slope_denominator;
        double intercept = point1.second - slope*point1.first;
        double under_radical_term = a*(intercept + slope*(vertex_x + a*slope));
        if( 0.0 > under_radical_term){//no real intercepts
            return;
        }
        x_pos = vertex_x +2.0*a*slope + 2.0*sqrt(under_radical_term);
        x_neg = vertex_x +2.0*a*slope - 2.0*sqrt(under_radical_term);

        //We figure out what z values to pick from the slope
        z_neg = slope*x_pos + intercept;
        z_pos = slope*x_pos + intercept;
        if( 0.0 > slope){//z_neg is paired with x_pos
            std::swap(z_neg,z_pos);
        }
    }

    //Test to see that the
    if(x_neg >= std::min(point1.first, point2.first) && x_neg <= std::max(point1.first, point2.first)
        && z_neg >= std::min(point1.second, point2.second) && z_neg <= std::max(point1.second, point2.second)){
        intersection_x_values.push_back(PGVector3(x_neg, y, z_neg));
    }
    if(x_pos >= std::min(point1.first, point2.first) && x_pos <= std::max(point1.first, point2.first)
        && z_pos >= std::min(point1.second, point2.second) && z_pos <= std::max(point1.second, point2.second)){
        intersection_x_values.push_back(PGVector3(x_pos, y, z_pos));
    }
};
/*

void ReconstructionParabola::draw_parabola(double a, double vertex_x, Rectangle rectanglePoints, size_t number_intersections, double y){

    static size_t number_drawn = 0;
    const size_t MAX_NUM_IMAGES = 100;
    ++number_drawn;
    if( MAX_NUM_IMAGES < number_drawn) return;

    char buffer[5000];
    TCanvas canvas("dc","Parabola and Phantom", 800, 600);
    TPad textPad("pad1","Data elements",0.02,0.02,0.24,0.98);
    TPad shapePad("pad2","Shapes",0.26,0.02,0.98,0.98);

    textPad.SetFillColor(kBlue - 10);
    shapePad.SetFillColor(kWhite);
    textPad.Draw();
    shapePad.Draw();

    shapePad.cd();

    TLine line;

    //find x min and max
    vector<double> tmpVec;
    tmpVec.push_back(rectanglePoints.x_max_z_max.first);
    tmpVec.push_back(rectanglePoints.x_min_z_max.first);
    tmpVec.push_back(rectanglePoints.x_max_z_min.first);
    tmpVec.push_back(rectanglePoints.x_min_z_min.first);
    tmpVec.push_back( vertex_x );

    //std::sort(tmpVec.begin(), tmpVec.end());
    //double x_range_min = tmpVec[0] - 0.1*abs(tmpVec[0]);
    //double x_range_max = 1.1*(tmpVec[tmpVec.size() - 1]);

    //find z min and max
    //tmpVec.clear();
    tmpVec.push_back(rectanglePoints.x_max_z_max.second);
    tmpVec.push_back(rectanglePoints.x_min_z_max.second);
    tmpVec.push_back(rectanglePoints.x_max_z_min.second);
    tmpVec.push_back(rectanglePoints.x_min_z_min.second);
    //tmpVec.push_back(b);
    //tmpVec.push_back(-b);
    std::sort(tmpVec.begin(), tmpVec.end());
    double z_range_min = tmpVec[0] - 0.1*abs(tmpVec[0]);
    double z_range_max =  1.1*tmpVec[tmpVec.size() - 1];

    sprintf(buffer,"sqrt((4.0*%f)*(x - %.30f))", a, vertex_x);
    TF1 parabola_positive("parabola",buffer,vertex_x*0.95,z_range_max);
    sprintf(buffer,"-sqrt((4.0*%f)*(x - %.30f))", a, vertex_x);
    //TF1 parabola_negative("parabola",buffer,vertex_x,z_range_max);
    TF1 parabola_negative("parabola",buffer,vertex_x*0.95,z_range_max);
    parabola_positive.SetLineWidth(2);
    parabola_positive.SetLineColor(kBlack);
    parabola_negative.SetLineWidth(2);
    parabola_negative.SetLineColor(kBlack);

    sprintf(buffer, "Parabola and Phantom: %d intersections; x(mm); z (mm)", number_intersections);
    TH1F axisHist("myHist",buffer, 1, z_range_min, z_range_max);
    axisHist.GetXaxis()->CenterTitle();
    axisHist.GetYaxis()->CenterTitle();
    axisHist.SetMinimum(z_range_min);
    axisHist.SetMaximum(z_range_max);
    double range_length = z_range_max - z_range_min;
    axisHist.SetStats(0);
    axisHist.Fill( (z_range_min+z_range_max)/2.0);
    axisHist.SetLineWidth(0);
    axisHist.SetLineColor(kWhite);
    axisHist.Draw();
    parabola_negative.Draw("SAME");
    parabola_positive.Draw("SAME");



    line.SetLineWidth(3);
    line.SetLineColor(kGreen);
    line.DrawLine(rectanglePoints.x_min_z_max.first, rectanglePoints.x_min_z_max.second,
                  rectanglePoints.x_max_z_max.first, rectanglePoints.x_max_z_max.second);
    line.SetLineColor(kRed);
    line.DrawLine(rectanglePoints.x_max_z_max.first, rectanglePoints.x_max_z_max.second,
                  rectanglePoints.x_max_z_min.first, rectanglePoints.x_max_z_min.second);
    line.SetLineColor(kGreen+3);
    line.DrawLine(rectanglePoints.x_max_z_min.first, rectanglePoints.x_max_z_min.second,
                  rectanglePoints.x_min_z_min.first, rectanglePoints.x_min_z_min.second);
    line.SetLineColor(kBlue+3);
    line.DrawLine(rectanglePoints.x_min_z_min.first, rectanglePoints.x_min_z_min.second,
                  rectanglePoints.x_min_z_max.first, rectanglePoints.x_min_z_max.second);
    line.SetLineWidth(1);
    line.DrawLine((rectanglePoints.x_min_z_max.first + rectanglePoints.x_max_z_max.first)/2.0,
                  (rectanglePoints.x_min_z_max.second + rectanglePoints.x_max_z_max.second)/2.0,
                  (rectanglePoints.x_min_z_min.first + rectanglePoints.x_max_z_min.first)/2.0,
                  (rectanglePoints.x_min_z_min.second + rectanglePoints.x_max_z_min.second)/2.0);
    line.DrawLine((rectanglePoints.x_min_z_max.first + rectanglePoints.x_min_z_min.first)/2.0,
                  (rectanglePoints.x_min_z_max.second + rectanglePoints.x_min_z_min.second)/2.0,
                  (rectanglePoints.x_max_z_max.first + rectanglePoints.x_max_z_min.first)/2.0,
                  (rectanglePoints.x_max_z_max.second + rectanglePoints.x_max_z_min.second)/2.0);



    textPad.cd();
    TLatex latex;
    latex.SetTextSize(latex.GetTextSize()*1.8);
    double text_x_position = z_range_min + range_length * 0.1;
    text_x_position = 0.05;
    //double text_z_position = z_range_min + range_length * 0.9;
    double z_increment = 0.06;
    double z_position = 0.9;

    sprintf(buffer,"Intersections: %d", number_intersections);
    latex.DrawLatex(text_x_position, z_position,buffer);
    z_position -= z_increment;

    sprintf(buffer,"y: %.2f", y);
    latex.DrawLatex(text_x_position, z_position,buffer);
    z_position -= z_increment;

    sprintf(buffer,"#phi: %.2f", phi_*180/M_PI);
    latex.DrawLatex(text_x_position, z_position,buffer);
    z_position -= z_increment;

    sprintf(buffer,"#alpha: %.2f", alpha_*180/M_PI);
    latex.DrawLatex(text_x_position,z_position,buffer);
    z_position -= z_increment;

    sprintf(buffer,"#theta: %.2f", theta1_*180/M_PI);
    latex.DrawLatex(text_x_position,z_position,buffer);
    z_position -= z_increment;

    sprintf(buffer,"tan(#alpha - #phi): %.2f", tan(alpha_ - phi_) );
    latex.DrawLatex(text_x_position,z_position,buffer);
    z_position -= z_increment;

    sprintf(buffer,"a: %.2f", a );
    latex.DrawLatex(text_x_position,z_position,buffer);
    z_position -= z_increment;

    sprintf(buffer,"vertex_x: %.2f", vertex_x );
    latex.DrawLatex(text_x_position,z_position,buffer);
    z_position -= z_increment;

    sprintf(buffer,"apex: (%.1f,%.1f,%.1f)",
      this->getComptonScatter().getConeApex().x,
      this->getComptonScatter().getConeApex().y,
      this->getComptonScatter().getConeApex().z);
    latex.DrawLatex(text_x_position, z_position,buffer);
    z_position -= z_increment;

    sprintf(buffer,"axis: (%.1f,%.1f,%.1f)",
      this->getComptonScatter().getConeAxis().x,
      this->getComptonScatter().getConeAxis().y,
      this->getComptonScatter().getConeAxis().z);
    latex.DrawLatex(text_x_position, z_position,buffer);
    z_position -= z_increment;


    sprintf(buffer,"parabola_%10d_%.1f_%.1f.png", int(phi_*180.0/M_PI), a, vertex_x);
    canvas.SaveAs(buffer);
}
*/

bool ReconstructionParabola::doesConicIntersectPhantom(){

  return true;

}


size_t  ReconstructionParabola::getRandomPointInPhantom0(PGVector3 &randPoint, size_t num_random_tries){

    PGVector3 new_point;
    PGVector3 apex = this->getComptonScatter().getConeApex();
    //PGVector3 axis = this->getComptonScatter().getConeAxis();

    randPoint.y = 0.0;
    PGVector3 min_intercept = getConeInnerYPlaneIntercept(randPoint.y);
    //PGVector3 axis_intercept = getConeAxisYPlaneIntercept(randPoint.y);
    PGVector3 min_intercept_ellipse_frame = min_intercept.translate(apex * -1.0);
    min_intercept_ellipse_frame = PGVector3::rotateYaxis(min_intercept_ellipse_frame, theta1_);

    // parabola_x_min is in the prime frame

    //double ellipse_x_min =  delta_y*( (alpha_ - phi_ > 0)? tan(alpha_ - phi_):tan(phi_ - alpha_ ) );
    double parabola_x_min = min_intercept_ellipse_frame.x;


    //need to rotate x,y around the z axis
    PGVector3 point_on_cone = min_intercept_ellipse_frame;
    //3 is arbitrary and has to be on the parabola because of our choice of the prime frame
    double arbitrary_shift = 3.0;
    point_on_cone.x += arbitrary_shift;

    //rotate to axis frame
    point_on_cone = PGVector3::rotateZaxis(point_on_cone, -alpha_);

    point_on_cone.z = ConicSection::solveConeEquation(point_on_cone.x, point_on_cone.y, phi_);

    //rotate back to the parabola frame
    point_on_cone = PGVector3::rotateZaxis(point_on_cone, alpha_);

    double a = 0.25*point_on_cone.z*point_on_cone.z/point_on_cone.x;

    vector<PGVector3> intersection_points;

    getLineIntercepts_(a, parabola_x_min, randPoint.y,
                      rotated_phantom_boundary_.x_min_z_min, rotated_phantom_boundary_.x_min_z_max,
                      intersection_points);
    getLineIntercepts_(a, parabola_x_min, randPoint.y,
                      rotated_phantom_boundary_.x_min_z_max, rotated_phantom_boundary_.x_max_z_max,
                      intersection_points);
    getLineIntercepts_(a, parabola_x_min, randPoint.y,
                      rotated_phantom_boundary_.x_max_z_max, rotated_phantom_boundary_.x_max_z_min,
                      intersection_points);
    getLineIntercepts_(a, parabola_x_min, randPoint.y,
                      rotated_phantom_boundary_.x_max_z_min, rotated_phantom_boundary_.x_min_z_min,
                      intersection_points);

    draw_parabola(a, parabola_x_min, rotated_phantom_boundary_, intersection_points.size(), randPoint.y);

    randPoint.x = 10;
    randPoint.y = 10;
    randPoint.z = 10;

    return 0;
};

///////////////////////////////////////////////////////////////////////////////
/// a member function which tests the transformations of points from
/// the phantom to the cone axis frame.
///
/// @see transformPointToConeAxisFrame(const PGVector3 point)
/// @return void
///////////////////////////////////////////////////////////////////////////////
void ReconstructionParabola::testTransformations(){
    //Find a test point on the cone.
    const PGVector3 cone_axis = compton_scatter_.getConeAxis();
    const PGVector3 cone_apex = compton_scatter_.getConeApex();
    PGVector3 test_point, test_result, test_result_apex, test_result_origin;


    // transform point ellipse min, the point on the ellipse
    // major axis which is closest to the origin in the phantom
    // frame, to the axis frame and then test to show that it is
    // still on the cone.
    double xprime = compton_scatter_.getConeApex().y * tan(alpha_ - phi_);
    test_point.y = 0.0;

    if(0.0 == cone_axis.x){
      test_point.x = 0;
      test_point.z = cone_axis.z;
    }else{
      //test_point.x = xprime*cos(atan(cone_axis.z/cone_axis.x));
      //test_point.z = -xprime*sin(atan(cone_axis.z/cone_axis.x));
      //angDiff = atan(cone_axis.z/cone_axis.x) - theta1_;
      test_point.x = xprime*cos(theta1_) + compton_scatter_.getConeApex().x;
      test_point.z = xprime*sin(theta1_) + compton_scatter_.getConeApex().z;
    }
    transformPointToConeAxisFrameFromPhantomFrame(test_point);

    double distance_off_cone = sqrt(test_result.x*test_result.x + test_result.z*test_result.z) - sqrt(tan(phi_)*tan(phi_)*test_result.y*test_result.y);

    cout<<"distance off cone (<0>): "<< distance_off_cone << endl;

    // test apex transformed to cone frame.
    // transformed point should be (0,0,0) and should
    // have distance 0.
    test_result_apex = cone_apex;
    transformPointToConeAxisFrameFromPhantomFrame(test_result_apex);
    double distance_off_cone_apex = test_result_apex.x*test_result_apex.x
                                  + test_result_apex.z*test_result_apex.z
                                  - tan(phi_)*tan(phi_)*test_result_apex.y*test_result_apex.y;

    cout<<"distance off cone apex(0): "<< distance_off_cone_apex << endl;


    //test a point on the cone axis. It should return x,z values ~0 but
    // distance from cone should be pretty large.
    transformPointToConeAxisFrameFromPhantomFrame(test_result_origin);
    double distance_off_cone_origin = test_result_origin.x*test_result_origin.x
                                  + test_result_origin.z*test_result_origin.z
                                  - tan(phi_)*tan(phi_)*test_result_origin.y*test_result_origin.y;

    cout<<"distance off cone origin: "<< distance_off_cone_origin << endl;

    //test a point on the cone axis. It should return x,z values ~0 but
    // distance from cone should be pretty large.
    double x = (-cone_apex.y)*cone_axis.x/cone_axis.y + cone_apex.x;
    double z = (-cone_apex.y)*cone_axis.z/cone_axis.y + cone_apex.z;
    PGVector3 axis_y_plane_intercept(x,0.0,z);
    transformPointToConeAxisFrameFromPhantomFrame(axis_y_plane_intercept);
}


void ReconstructionParabola::print(){
    printf("-------ReconstructionParabola---------\n");
    printf("cone apex(%.2f, %.2f, %.2f)\n", compton_scatter_.getConeApex().x, compton_scatter_.getConeApex().y, compton_scatter_.getConeApex().z);
    printf("cone axis(%.2f, %.2f, %.2f)\n", compton_scatter_.getConeAxis().x, compton_scatter_.getConeAxis().y, compton_scatter_.getConeAxis().z);
    printf("gamma origin(%.2f, %.2f, %.2f)\n", likely_origin_.x, likely_origin_.y,likely_origin_.z);
    printf("phi (%.2f)\n", phi_);
    printf("alpha (%.2f)\n", alpha_);
    printf("theta1 (%.2f)\n", theta1_);
    printf("theta2 (%.2f)\n", theta2_);
    printf("--------------------------------------\n");
};
