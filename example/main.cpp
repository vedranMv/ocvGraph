/**
 *  Examples of using OCVGraph class for plotting data
 *
 *  First example plots data read by e.g. distance radar sensor into a
 *  polar coordinate system. Length of ray represents a distance to potential
 *  obstacle while angle is angle of sensor at the particular measurement.
 *  Half-circles give a way of relating ray distance to real-life measurements.
 *  Second example shows how to clear graph on run-time, set new center and
 *  plot a parabola
 */

#include <iostream>
#include <cstdlib>
#include <cstdint>

#include "../ocvgraph.h"

//  Initialize graph with given image size
OCVGraph plot(200, 400);
double scale = 2.5;

using namespace std;


int main(int argc, char **argv)
{
    srand(time(NULL));

    /*
     *  Example 1: Plot random lines in polar coordinates, starting from origin
     */
     // Set center to be in the middle of bottom edge
    plot.SetCenter(200, 200);
    //  Plot circles in the origin of coordinate system
    plot.Circle(30*scale);
    plot.Text("30cm", cv::Point2i(30*scale-20,2));
    plot.Circle(60*scale);
    plot.Text("60cm", cv::Point2i(60*scale-20,2));
    plot.Circle(80*scale);
    plot.Text("80cm", cv::Point2i(80*scale-20,2));

    //  Add some random polar lines to the plot from origin
    for (uint16_t i = 10; i <170; i++)
        plot.LinePolar((double)i, (double)(rand() % 80) * scale);

    //  Plot graph
    cv::imshow("Plot", plot.GetMatImg());
    cv::waitKey();

    /*
     *  Example 2: Plot parabola using PolyN function
     */
     // Clear old plot
    plot.Clear();
    //  Shift origin 200px from left edge and 100 pixels from top edge
    plot.SetCenter(200, 100);
    //  Plot circle with center at origin
    plot.Circle(20);

    //  Plot parabola with equation y=2x^2+x-1
    vector<double> cof;
    cof.push_back(-1);   //  Coefficient with x^0
    cof.push_back(1);   //  Coefficient with x^1
    cof.push_back(0.3);   //  Coefficient with x^2
    plot.PolyN(cof);

    cv::imshow("Plot", plot.GetMatImg());
    cv::waitKey();

    return 0;
}
