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

    /***************************************************************************
     **  Example 1: Plot random lines in polar coordinates, starting from origin
     **************************************************************************/
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

    /***************************************************************************
     ****************************  Example 2: Plot parabola using PolyN function
     **************************************************************************/
     // Clear old plot and make origin of graph in center of the image
    plot.Clear(true);
    //  Plot data axes with divisions every 10 units in X & Y direction
    plot.AddAxes(10, 10);

    //  Plot parabola with equation y=0.1x^2-20*x+940 on interval [+80, +120]
    vector<double> cof;
    cof.push_back(940);     //  Coefficient with x^0
    cof.push_back(-20);     //  Coefficient with x^1
    cof.push_back(0.1);     //  Coefficient with x^2
    plot.PolyN(cof, 80, 120);

    //  Plot graph
    cv::imshow("Plot", plot.GetMatImg());
    cv::waitKey();

    /***************************************************************************
     ****************************  Example 3: Append couple of lines to the plot
     **************************************************************************/
    //  Plot equation y=-2x
    cof.clear();
    cof.push_back(0);       //  Coefficient with x^0
    cof.push_back(-0.5);    //  Coefficient with x^1
    plot.PolyN(cof);

    //  Plot equation y=0.3x+10, on interval from -70 to 50
    cof.clear();
    cof.push_back(10);     //  Coefficient with x^0
    cof.push_back(0.3);    //  Coefficient with x^1
    plot.PolyN(cof, -70, 50);

    //  Plot equation y=1x+12, on interval from -50 to 50
    cof.clear();
    cof.push_back(12);      //  Coefficient with x^0
    cof.push_back(1);       //  Coefficient with x^1
    plot.PolyN(cof, -50, 50);

    plot.AddAxes(10, 10);
    cv::imshow("Plot", plot.GetMatImg());
    cv::waitKey();

    return 0;
}
