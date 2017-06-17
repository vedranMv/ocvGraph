#include "ocvgraph.h"
#include <vector>
#include <cmath>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

///-----------------------------------------------------------------------------
///                      Class constructor & destructor                 [PUBLIC]
///-----------------------------------------------------------------------------

OCVGraph::OCVGraph(OCVGraph &graph)
    : _graphHolder(graph._graphHolder), _center(graph._center)
{
    _background = cv::Scalar::all(255);
}

OCVGraph::OCVGraph(cv::Mat &graphplot)
    : _graphHolder(graphplot), _center(graphplot.cols/2, graphplot.rows/2)
{
    _background = cv::Scalar::all(255);
}

OCVGraph::OCVGraph(int height, int width, const cv::Scalar& s)
    : _graphHolder(height+1, width+1, CV_8UC3, s), _center(width/2, height/2)
{
    _background = s;
}

OCVGraph::~OCVGraph()
{}

/**
 * @brief Set center of a graph in image coordinate (which pixel corresponds to 0,0)
 * @param uc Pixel's x-coordinate (along width of the image)
 * @param vc Pixel's y-coordinate (along height of the image)
 */
void OCVGraph::SetCenter(int uc, int vc)
{
    _center = cv::Point2i(uc, vc);
}
void OCVGraph::SetCenter(cv::Point2i center)
{
    _center = center;
}

/**
 * @brief Get center of a graph(0,0)
 * @return Center point(0,0) of graph in image coordinates
 */
cv::Point2i& OCVGraph::GetCenter()
{
    return _center;
}

/**
 * @brief Draw a line in Cartesian coordinate system between points p1 and p2
 * @param p1 Starting point of the line in graph coordinates(!!)
 * @param p2 Ending point of the line in graph coordinates(!!)
 */
void OCVGraph::LineCartesian(cv::Point2i p1, cv::Point2i p2, cv::Scalar color,
                             int thickness, int lineType, int shift)
{
    cv::line(_graphHolder, _XYtoUV(p1), _XYtoUV(p2), color, thickness,
             lineType, shift);
}

/**
 * @brief Draw a line from p1 in polar coordinate system given angle and radius
 * @param angle Angle in degrees
 * @param rad
 * @param p1 Starting point of the line in graph coordinates(!!)
 * @param color
 */
void OCVGraph::LinePolar(double angle, double rad, cv::Point2i p1, cv::Scalar color,
                         int thickness, int lineType, int shift)
{
    //  Convert angle to radians
    angle = angle * 3.14159625 / 180;
    //  Calculate ending point in Cartesian coordinates
    cv::Point2i p2(rad * cos(angle) + p1.x, rad * sin(angle) + p1.y);

    //  Plot line in Cartesian coordinate
    LineCartesian(p1, p2, color, thickness, lineType, shift);
}

/**
 * @brief Draw a circle with given center, radius and (optionally) color
 * @param rad Radius of circle
 * @param color
 * @param p1 Center of circle in graph coordinates
 */
void OCVGraph::Circle(double rad, cv::Point2i p1, cv::Scalar color,
                      int thickness, int lineType, int shift)
{
    cv::circle(_graphHolder, _XYtoUV(p1), rad, color, thickness,
               lineType, shift);
}

void OCVGraph::PolyN(std::vector<double>&coefs, double xmin, double xmax)
{
    cv::Point2i lastP;
    bool firstEntry = false;
    //  Function boundaries in image coordinates
    int umin = 0, umax = _graphHolder.cols;

    //  Update user-provided boundaries if given
    if (xmin != xmax)
    {
        umin = _XYtoUV(cv::Point2i(xmin, 0)).x;
        umax = _XYtoUV(cv::Point2i(xmax, 0)).x;
    }

    //  Compute function value only for X within image coordinates
    for (int iu = umin; iu < umax; iu++)
    {
        //  Function value saved here
        double y = 0;
        //  Translate image coordinate U into graph coordinate 'x'
        cv::Point2d tmp = _UVtoXY(cv::Point2i(iu, 0));

        //  Loop through all coefficients and calculate 'y'
        for (int i = (coefs.size()-1); i >= 0; i--)
            y += coefs[i] * pow(tmp.x, (double)(i));

        //  If within image boundaries turn pixel black to construct curve
        cv::Point2i imgC = _XYtoUV(cv::Point2d(tmp.x, y));

        if ((imgC.y >= 0) && (imgC.y < _graphHolder.rows) && firstEntry)
        {
            cv::line(_graphHolder, imgC, lastP, cv::Scalar::all(0));
            //cv::Vec3b &pixel = _graphHolder.at<cv::Vec3b>(imgC.y, imgC.x);
            //pixel[0] = pixel[1] = pixel[2] = 0;
        }
        else if ((imgC.y >= 0) && (imgC.y < _graphHolder.rows) && !firstEntry)
        {
            firstEntry = true;
        }
        else
            continue;

        lastP = imgC;
    }
}

/**
 * @brief Add text to a specific point in graph
 * @param txt
 * @param p1
 */
void OCVGraph::Text(std::string txt, cv::Point2i p1)
{
    cv::putText(_graphHolder, txt, _XYtoUV(p1), 2, 0.3, cv::Scalar::all(128));
}

/**
 *  Export current graph image into an image file at given location
 *  @param path Path to image file
 */
void OCVGraph::Export(std::string path)
{
    cv::imwrite(path, _graphHolder);
}

/**
 * @brief Get image of a graph currently stored in this object
 * @return Mat type image of current graph
 */
cv::Mat& OCVGraph::GetMatImg()
{
    return _graphHolder;
}

void OCVGraph::Clear(bool resetCenter)
{
    cv::Mat tmp (_graphHolder.rows, _graphHolder.cols, CV_8UC3, _background);

    tmp.copyTo(_graphHolder);

    if (resetCenter)
        _center = cv::Point2i(_graphHolder.cols/2, _graphHolder.rows/2);
}

/**
 * @brief Convert U,V pixel coordinates into a corresponding X,Y graph coordinates
 * @param imageCoord Coordinates of pixel in image coordinates
 * @return cv::Point2i X,Y graph coordinates
 */
cv::Point2d OCVGraph::_UVtoXY(cv::Point2i imageCoord)
{
    return cv::Point2i(imageCoord.x-_center.x, _center.y-imageCoord.y);
}

/**
 * @brief Convert graph coordinates to image coordinates
 * @param graphCoord
 * @return Image coordinates corresponding to graph coordinates
 */
cv::Point2i OCVGraph::_XYtoUV(cv::Point2d graphCoord)
{
    return cv::Point2i(graphCoord.x+_center.x, _center.y-graphCoord.y);
}


OCVGraph& OCVGraph::operator=(OCVGraph& arg)
{
    _graphHolder = arg._graphHolder;
    _center = arg._center;
    return *this;
}

OCVGraph& OCVGraph::operator=(OCVGraph arg)
{
    _graphHolder = arg._graphHolder;
    _center = arg._center;
    return *this;
}
