/**
 *  Examples of using OCVGraph class for plotting data
 *
 *  First example plots data read by e.g. distance radar sensor into a
 *  polar coordinate system. Length of ray represents a distance to potential
 *  obstacle while angle is angle of sensor at the particular measurement.
 *  Half-circles give a way of relating ray distance to real-life measurements.
 *  Second example shows how to clear graph on run-time, set new center and
 *  plot a parabola.
 *  Third example appends three lines of different color to the graph from
 *  example 2, and draws a legend showing which color graph represents
 *  which function.
 */

#include <iostream>
#include <cstdlib>
#include <cstdint>

#include "../ocvgraph.h"

//  Initialize graph with given image size
OCVGraph plot(400, 600);

using namespace std;


int main(int argc, char **argv)
{
    srand(time(NULL));

    /***************************************************************************
     **  Example 1: Plot random lines in polar coordinates, starting from origin
     **************************************************************************/
     // Set center to be in the middle of bottom edge
    plot.SetCenter(300, 400);
    //  Set X range of graph to be between +/-90 and Y automatic
    plot.SetXRangeKeepAspectR(-90, 90);
    //  Plot circles in the origin of coordinate system
    plot.Circle(30);
    plot.Circle(60);
    plot.Circle(80);

    //  Add some random polar lines to the plot from origin
    for (uint16_t i = 10; i <170; i++)
        plot.LinePolar((double)i, (double)(rand() % 80));

    //  Plot axes with division every 10 steps, every 2nd division line larger
    plot.AddAxes(10, 10, 2, 2);
    //  Plot graph
    cv::imshow("Plot", plot.GetMatImg());
    cv::waitKey();

    /***************************************************************************
     ****************************  Example 2: Plot parabola using PolyN function
     **************************************************************************/
     // Clear old plot and make origin of graph in center of the image
    plot.Clear(true);
    //  Move center towards top
    plot.SetCenter(300, 70);
    //  Define range of X axis, keep Y proportional, move center in X if needed
    plot.SetXRangeKeepAspectR(-100, 200, true);
    //  Plot axes with division every 10 steps, every 4th division line larger
    plot.AddAxes(10, 10, 4, 4);

    //  Plot parabola with equation y=0.1x^2-20x+940 on interval [+80, +120]
    vector<double> cof;
    cof.push_back(940);     //  Coefficient with x^0
    cof.push_back(-20);     //  Coefficient with x^1
    cof.push_back(0.1);     //  Coefficient with x^2
    plot.PolyN(cof, 80, 120, COLOR_BLACK);

    plot.AddToLegend(0, COLOR_BLACK, "f(x)=0.1x^2-20x+940");

    //  Plot graph
    cv::imshow("Plot", plot.GetMatImg());
    cv::waitKey();

    /***************************************************************************
     ***************  Example 3: Append couple of lines to the plot & add legend
     **************************************************************************/
    //  Plot equation y=-0.5x
    cof.clear();
    cof.push_back(0);       //  Coefficient with x^0
    cof.push_back(-0.5);    //  Coefficient with x^1
    plot.PolyN(cof,0, 0, COLOR_YELLOW);

    //  Plot equation y=0.3x+10, on interval from -70 to 50
    cof.clear();
    cof.push_back(10);     //  Coefficient with x^0
    cof.push_back(0.3);    //  Coefficient with x^1
    plot.PolyN(cof, -70, 50, COLOR_PURPLE);

    //  Plot equation y=1x+12, on interval from -50 to 50
    cof.clear();
    cof.push_back(12);      //  Coefficient with x^0
    cof.push_back(1);       //  Coefficient with x^1
    plot.PolyN(cof, -50, 50, COLOR_BROWN);

    //  Add legend entries and plot legend
    plot.AddToLegend(1, COLOR_YELLOW, "f(x)=-0.5x");
    plot.AddToLegend(2, COLOR_PURPLE, "f(x)=0.3x+10");
    plot.AddToLegend(3, COLOR_BROWN, "f(x)=1x+12");
    plot.AppendLegend();

    cv::imshow("Plot", plot.GetMatImg());
    cv::waitKey();

    return 0;
}
