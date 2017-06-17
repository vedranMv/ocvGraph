/**
 *  Date created: 1.6.2017
 *  Author: Vedran Mikov
 *
 *  OpenCV-based simple graph plotter. Allows user to set center of the graph at
 *  any pixel in the image and then plot simple curves in reference to that point
 *  as origin of coordinate system.
 *
 *  Current API tries to follow OpenCV argument naming and order as much as possible
 *
 *  @version 1.2.0
 *  V1.0.0 - 1.6.2017
 *  +Custom origin of coordinate system at any point of image. Support for
 *  plotting circles, lines (in Cartesian and polar coordinate system) and
 *  writing text.
 *  V1.1.0
 *  +Save plot to a given path
 *  V1.2.0 - 17.6.2017
 *  +Added clear function to reset the graph
 *  +Added option to plot N-degree polynomials in Cartesian, within given range
 *  +Show data axes, variable number of divisions for X and Y axes
 *
 *  TODO:
 *  +Add ability to plot arbitrary functions
 *  +Add support for zoom in/out(scaling)
 */

#ifndef OCVGRAPH_H
#define OCVGRAPH_H

#include <opencv2/opencv.hpp>
#include <string>


class OCVGraph
{
public:
    OCVGraph(OCVGraph &graph);
    OCVGraph(cv::Mat &graphplot);
    OCVGraph(int height, int width, const cv::Scalar &s =cv::Scalar::all(255));
    ~OCVGraph();

    void            SetCenter(int uc, int vc);
    void            SetCenter(cv::Point2i center);
    cv::Point2i&    GetCenter();

    void            SetScale(double scale);
    double          GetScale();

    void LineCartesian(cv::Point2i p1, cv::Point2i p2 = cv::Point2i(0,0),
                       cv::Scalar color = cv::Scalar::all(0),
                       int thickness=1, int lineType=8, int shift=0);

    void LinePolar(double angle, double rad, cv::Point2i p2 = cv::Point2i(0,0),
                   cv::Scalar color = cv::Scalar::all(0),
                   int thickness=1, int lineType=8, int shift=0);

    void Circle(double rad, cv::Point2i p1 = cv::Point2i(0,0),
                cv::Scalar color = cv::Scalar::all(0), int thickness=1,
                int lineType=8, int shift=0);

    void PolyN(std::vector<double>&coefs, double xmin=0, double xmax=0);

    void Text(std::string txt, cv::Point2i p1, int fontFace=2, double scale=0.3,
              cv::Scalar color=cv::Scalar::all(128), int thickness=1,
              int lineType=8, bool blo=false);

    void AddAxes(double xticks = 1.0, double yticks = 1.0);


    void        Export(std::string path);
    cv::Mat&    GetMatImg();
    void        Clear(bool resetCenter = false);

    OCVGraph& operator=(OCVGraph& arg);
    OCVGraph& operator=(OCVGraph arg);


private:
    cv::Point2d _UVtoXY(cv::Point2i imageCoord);
    cv::Point2i _XYtoUV(cv::Point2d graphCoord);

private:
    cv::Mat         _graphHolder;
    cv::Point2i     _center;
    cv::Scalar      _background;
};

#endif // OCVGRAPH_H
