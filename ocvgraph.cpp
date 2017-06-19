#include "ocvgraph.h"
#include <cmath>


#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

std::string dtos(double arg)
{
    std::string retVal;

    if (arg < 0.0)
        retVal += '-';

    arg = std::abs(arg);
    retVal = std::to_string((int)arg);

    //  Have decimals, extract max two
    if (arg-(trunc(arg)) != 0.0)
    {
        retVal += '.';
        retVal += std::to_string((int)((arg-trunc(arg))*100.0));
    }

    return retVal;
}

inline bool isModZero(const double x, const double y)
{
    double decp, intp;

    decp = modf (x/y , &intp);

    return (decp == 0) || (decp > 0.99999);
}

///-----------------------------------------------------------------------------
///                      Class constructor & destructor                 [PUBLIC]
///-----------------------------------------------------------------------------

OCVGraph::OCVGraph(OCVGraph &graph)
    : _graphHolder(graph._graphHolder), _center(graph._center)
{
    _background = cv::Scalar::all(255);
    _xlim[0] = -graph._graphHolder.cols/2;
    _xlim[1] =  graph._graphHolder.cols/2;
    _ylim[0] = -graph._graphHolder.rows/2;
    _ylim[1] =  graph._graphHolder.rows/2;
}

OCVGraph::OCVGraph(cv::Mat &graphplot)
    : _graphHolder(graphplot), _center(graphplot.cols/2, graphplot.rows/2)
{
    _background = cv::Scalar::all(255);
    _xlim[0] = -graphplot.cols/2;
    _xlim[1] =  graphplot.cols/2;
    _ylim[0] = -graphplot.rows/2;
    _ylim[1] =  graphplot.rows/2;
}

OCVGraph::OCVGraph(int height, int width, const cv::Scalar& s)
    : _graphHolder(height, width, CV_8UC3, s), _center(width/2, height/2)
{
    _background = s;
    _xlim[0] = -width/2;
    _xlim[1] =  width/2;
    _ylim[0] = -height/2;
    _ylim[1] =  height/2;
}

OCVGraph::~OCVGraph()
{}

///-----------------------------------------------------------------------------
///                      Manipulation of center of graph                [PUBLIC]
///-----------------------------------------------------------------------------

/**
 * @brief Set center of a graph in image coordinate (which pixel corresponds to 0,0)
 * @note Resets axis limits
 * @param uc Pixel's x-coordinate (along width of the image)
 * @param vc Pixel's y-coordinate (along height of the image)
 */
void OCVGraph::SetCenter(int uc, int vc)
{
    _center = cv::Point2i(uc, vc);

    //  Adjust limits
    _xlim[0] = -uc;
    _xlim[1] =  _graphHolder.cols-uc;
    _ylim[0] = -(_graphHolder.rows-vc);
    _ylim[1] = vc;
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

///-----------------------------------------------------------------------------
///                      Manipulation of scale of graph                 [PUBLIC]
///-----------------------------------------------------------------------------

void OCVGraph::SetXRange(/*double xmin, */double xmax)
{
    //_xlim[0] = xmin;
    _xlim[1] = xmax;
}
void OCVGraph::SetYRange(/*double ymin, */double ymax)
{
    //_ylim[0] = ymin;
    _ylim[1] = ymax;
}

///-----------------------------------------------------------------------------
///                      Plotting functions                             [PUBLIC]
///-----------------------------------------------------------------------------

/**
 * @brief Draw a line in Cartesian coordinate system between points p1 and p2
 * @param p1 Starting point of the line in graph coordinates(!!)
 * @param p2 Ending point of the line in graph coordinates(!!)
 */
void OCVGraph::LineCartesian(cv::Point2d p1, cv::Point2d p2, cv::Scalar color,
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
void OCVGraph::LinePolar(double angle, double rad, cv::Point2d p1, cv::Scalar color,
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
void OCVGraph::Circle(double rad, cv::Point2d p1, cv::Scalar color,
                      int thickness, int lineType, int shift)
{
    //  Need to verify radius of ellipse in graph coordinates for each axis in
    //  case axes have different ranges -> result is ellipse
    cv::Point2i radUV = _XYtoUV(cv::Point2i(rad, rad));
    //  Calculate radius of ellipse in X and Y direction (based on scale of graph)
    cv::Size2i axes = cv::Size2i(-_center.x+(radUV.x), _center.y-(radUV.y));

    cv::ellipse(_graphHolder, _XYtoUV(p1), axes, 0, 0, 360, color, thickness,
                lineType, shift);
}

/**
 *  Plot polynomial of form coefs[n-1]*x^(n-1)+coefs[n-2]*x^(n-2)+...+coefs[0]*x^(0)
 *  in range between xmin and xmax. In case xmin and xmax are equal function is
 *  plotted on interval visible in the graph
 *  @param coefs Vector of polynomial coefficients. Coefficient at index i
 *  is multiplied by x^i
 *  @param xmin Lower bound of x while computing y
 *  @param xmax Upper bound of x while computing y
 */
void OCVGraph::PolyN(std::vector<double>&coefs, double xmin, double xmax,
                     cv::Scalar color, int thickness, int lineType, int shift)
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
            cv::line(_graphHolder, imgC, lastP, color, thickness, lineType, shift);
        else if ((imgC.y >= 0) && (imgC.y < _graphHolder.rows) && !firstEntry)
            firstEntry = true;
        else
            continue;

        lastP = imgC;
    }
}

/**
 * @brief Add text to a specific point in graph
 * @param txt
 * @param p1 Origin point of text
 * @param
 */
void OCVGraph::Text(std::string txt, cv::Point2d p1, int fontFace, double scale,
                    cv::Scalar color, int thickness, int lineType, bool blo)
{
    cv::putText(_graphHolder, txt, _XYtoUV(p1), fontFace, scale, color,
                thickness, lineType, blo);
}

/**
 *  Plot data axes together with values
 *  @param xticks Number of steps between division lines on X-axis
 *  @param yticks Number of steps between division lines on Y-axis
 */
void OCVGraph::AddAxes(double xticks, double yticks, int xlarge, int ylarge)
{
    //  Boundary 1: (+Xmax, -Ymax)
    cv::Point2d bound1 = _UVtoXY(cv::Point2i(_graphHolder.cols, _graphHolder.rows));
    //  Boundary 2: (-Xmax, +Ymax)
    cv::Point2d bound2 = _UVtoXY(cv::Point2i(0, 0));

    //  Add axis lines
    LineCartesian(cv::Point2d(bound1.x, 0), cv::Point2d(bound2.x, 0));
    LineCartesian(cv::Point2d(0, bound1.y), cv::Point2d(0, bound2.y));

    double xmax = std::max(std::abs(bound1.x),std::abs(bound2.x));
    double ymax = std::max(std::abs(bound1.y),std::abs(bound2.y));

    int tickCounter = 0;
    for (double ix = 0; ix < xmax; ix+=xticks)
    {
        double len;

        //  Larger division line is added every 'xlarge' ticks
        if ((tickCounter % xlarge) == 0)
        {
            //  Holds number of current tick in string format
            std::string txt;
            int baseline = 0;
            double dec;

            //  Size of larger division line
            len = (2*_ylim[1])/100.0;   //len = (-_ylim[0]+_ylim[1])/100.0;

            //  Add label with number 'ix'
            txt = dtos(ix);
            //  Calculate dimensions that text occupies (in image coord.)
            cv::Size txtSize = cv::getTextSize(txt, 2, 0.4, 1, &baseline);
            //  Convert image coord. to graph coord. (to match scaling)
            cv::Point2d txtimc = _UVtoXY(cv::Point2i(_center.x+txtSize.width/2, _center.y-txtSize.height));
            //  Appends text
            Text(txt, (cv::Point2d(ix, -1.5*len)-txtimc), 2, 0.4);

            //  Add label with number '-ix'
            txt = dtos(-ix);
            //  Calculate dimensions that text occupies (in image coord.)
            txtSize = cv::getTextSize(txt, 2, 0.4, 1, &baseline);
            //  Convert image coord. to graph coord. (to match scaling)
            txtimc = _UVtoXY(cv::Point2i(_center.x+txtSize.width/2, _center.y-txtSize.height));
            //  Appends text
            Text(txt, (cv::Point2d(-ix, -1.5*len)-txtimc), 2, 0.4);
        }
        else
            //  Size of smaller division line
            len = (2*_ylim[1])/220.0;

        LineCartesian(cv::Point2d( ix, len), cv::Point2d( ix, -len));
        LineCartesian(cv::Point2d(-ix, len), cv::Point2d(-ix, -len));

        tickCounter++;
    }

    tickCounter = 0;
    for (double iy = 0; iy < ymax; iy+=yticks)
    {
        double len;

        //  Larger division line is added every 'ylarge' ticks
        if ((tickCounter % ylarge) == 0)
        {
            std::string txt;
            int baseline = 0;
            double dec;
            //  Size of larger division line
            len = (2*_xlim[1])/100.0;   //len = (-_xlim[0]+_xlim[1])/100.0;

            //  Add label with number 'iy'
            txt = dtos(iy);
            //  Calculate dimensions that text occupies (in image coord.)
            cv::Size txtSize = cv::getTextSize(txt, 2, 0.4, 1, &baseline);
            //  Convert image coord. to graph coord. (to match scaling)
            cv::Point2d txtimc = _UVtoXY(cv::Point2i(_center.x+txtSize.width/2, _center.y-txtSize.height/2));
            txtimc.y = -txtimc.y;
            //  Appends text
            Text(txt, (cv::Point2d(len, iy)+txtimc), 2, 0.4);

            //  Add label with number '-iy'
            txt = dtos(-iy);
            //  Calculate dimensions that text occupies (in image coord.)
            txtSize = cv::getTextSize(txt, 2, 0.4, 1, &baseline);
            //  Convert image coord. to graph coord. (to match scaling)
            txtimc = _UVtoXY(cv::Point2i(_center.x+txtSize.width/2, _center.y-txtSize.height/2));
            txtimc.y = -txtimc.y;
            //  Appends text
            Text(txt, (cv::Point2d(len, -iy)+txtimc), 2, 0.4);
        }
        else
            //  Size of smaller division line
            len = (2*_xlim[1])/220.0;

        LineCartesian(cv::Point2d(len,  iy), cv::Point2d(-len,  iy));
        LineCartesian(cv::Point2d(len, -iy), cv::Point2d(-len, -iy));

        tickCounter++;
    }
}

void OCVGraph::AddToLegend(uint8_t index, cv::Scalar color, std::string txt,
                            int fontFace, double scale)
{
    _legend.push_back(make_tuple(index, color, txt, fontFace, scale));
}

void OCVGraph::AppendLegend(legLoc location)
{
    static const auto colOffset = cv::Point2i(0, 15),
                      lineLength = cv::Point2i(20, 0),
                      linTxtDist = cv::Point2i(5, 0);
    cv::Point2i anchor;
    uint8_t legMaxLen;

    //  Sort _legend array based on first index
    cv::Size txtSize;
    legMaxLen = std::get<2>(_legend[0]).length();
    for (uint8_t i = 0; i < _legend.size(); i++)
        for (uint j = 1; j < _legend.size()-i; j++)
    {
        //  Note the longest length of the string
        int baseline=0;
        txtSize = cv::getTextSize(std::get<2>(_legend[j]),
                std::get<3>(_legend[j]), std::get<4>(_legend[j]), 1, &baseline);

        if (txtSize.width > legMaxLen)
            legMaxLen = txtSize.width;

        //  Swap elements
        if (std::get<0>(_legend[j-1]) > std::get<0>(_legend[j]))
        {
            legEntry tmp = _legend[j-1];
            _legend[j-1] = _legend[j];
            _legend[j] = tmp;
        }
    }

    //  Calculate anchor point based on where the legend needs to be
    switch(location)
    {
        case TopLeft:
            anchor = cv::Point2i(10, 10);
        break;

        case TopRight:
            anchor = cv::Point2i(_graphHolder.cols-10, 10);
            anchor = anchor - (lineLength+linTxtDist);
            anchor.x -= legMaxLen;
        break;

        case BottomLeft:
            anchor = cv::Point2i(10, _graphHolder.rows-10);
            anchor.y -= (colOffset.y+txtSize.height)*(_legend.size()-1)+txtSize.height;
        break;

        case BottomRight:
            anchor = cv::Point2i(_graphHolder.cols-10, _graphHolder.rows-10);
            anchor.y -= (colOffset.y+txtSize.height)*(_legend.size()-1)+txtSize.height;
            anchor = anchor - (lineLength+linTxtDist);
            anchor.x -= legMaxLen;
        break;
    }

    //	Add legend line-by-line
    for (auto X : _legend)
    {
        //  Calculate size of the text
        int baseline=0;
        cv::Size txtSize = cv::getTextSize(std::get<2>(X), std::get<3>(X),
                                           std::get<4>(X), 1, &baseline);

        //  Text halflength in each dimension
        cv::Point2i txtOff = cv::Point2i(0, txtSize.height/2);

        //	Draw short line of given color
        cv::line(_graphHolder, anchor, anchor+lineLength, std::get<1>(X), 2);

        //	Write text corresponding to the line
        cv::putText(_graphHolder, std::get<2>(X), anchor+lineLength+linTxtDist+txtOff,
                    std::get<3>(X), std::get<4>(X), COLOR_BLACK);

        anchor = anchor + colOffset;
    }
}

///-----------------------------------------------------------------------------
///                      Miscellaneous functions                        [PUBLIC]
///-----------------------------------------------------------------------------

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

/**
 *  Clear current graph and (if requested) reset its origin to center of image
 *  @param resetCenter Whether or not to reset current center of graph to center
 *  of underlying cv::Mat image file
 */
void OCVGraph::Clear(bool resetCenter)
{
    cv::Mat tmp (_graphHolder.rows, _graphHolder.cols, CV_8UC3, _background);

    tmp.copyTo(_graphHolder);

    if (resetCenter)
        _center = cv::Point2i(_graphHolder.cols/2, _graphHolder.rows/2);
}

///-----------------------------------------------------------------------------
///                      Operator definitions                           [PUBLIC]
///-----------------------------------------------------------------------------

OCVGraph& OCVGraph::operator=(OCVGraph& arg)
{
    _graphHolder = arg._graphHolder;
    _center = arg._center;
    _background = arg._background;
    _legend = arg._legend;
    memcpy(_xlim, arg._xlim, 2*sizeof(double));
    memcpy(_ylim, arg._ylim, 2*sizeof(double));

    return *this;
}

OCVGraph& OCVGraph::operator=(OCVGraph arg)
{
    _graphHolder = arg._graphHolder;
    _center = arg._center;
    _background = arg._background;
    _legend = arg._legend;
    memcpy(_xlim, arg._xlim, 2*sizeof(double));
    memcpy(_ylim, arg._ylim, 2*sizeof(double));

    return *this;
}

///-----------------------------------------------------------------------------
///                      Unit conversion functions                     [PRIVATE]
///-----------------------------------------------------------------------------

/**
 * @brief Convert U,V pixel coordinates into a corresponding X,Y graph coordinates
 * @param imageCoord Coordinates of pixel in image coordinates
 * @return cv::Point2i X,Y graph coordinates
 */
cv::Point2d OCVGraph::_UVtoXY(cv::Point2i imageCoord)
{
    //cv::Point2d retVal = cv::Point2i(imageCoord.x-_center.x, _center.y-imageCoord.y);
    cv::Point2d retVal;

    retVal = cv::Point2i(imageCoord.x-_center.x, _center.y-imageCoord.y);
    retVal.x = _xlim[1]*retVal.x/(double)_center.x;
    retVal.y = _ylim[1]*retVal.y/(double)_center.y;
//    if (imageCoord.x > _center.x)
//        retVal.x = _xlim[1]*(double)imageCoord.x/((double)_graphHolder.cols-(double)_center.x);
//    else if (imageCoord.x <= _center.x)
//        retVal.x = -_xlim[0]*((double)imageCoord.x-(double)_center.x)/((double)_center.x);
//
//    if (imageCoord.y < _center.y)
//        retVal.y = _ylim[1]*(double)imageCoord.y/((double)_graphHolder.rows-(double)_center.y);
//    else if (imageCoord.x >= _center.x)
//        retVal.x = -_ylim[0]*((double)imageCoord.y-(double)_center.y)/((double)_center.y);

    return retVal;
}

/**
 * @brief Convert graph coordinates to image coordinates
 * @param graphCoord
 * @return Image coordinates corresponding to graph coordinates
 */
cv::Point2i OCVGraph::_XYtoUV(cv::Point2d graphCoord)
{
    //cv::Point2i retVal = cv::Point2i(graphCoord.x+_center.x, _center.y-graphCoord.y);
    cv::Point2i retVal;

    retVal.x = (int)(((double)_center.x)*graphCoord.x/_xlim[1]+((double)_center.x));
    retVal.y = (int)(((double)_center.y)-((double)_center.y)*graphCoord.y/_ylim[1]);

//    if (graphCoord.x > 0)
//        retVal.x = _center.x + ((double)(_graphHolder.cols-_center.x))*graphCoord.x/_xlim[1];
//    else if (graphCoord.x <= 0)
//        retVal.x = _center.x - ((double)(_center.x))*graphCoord.x/_xlim[0];
//
//    if (graphCoord.y < 0)
//        retVal.y = _center.y + ((double)(_graphHolder.rows-_center.y))*graphCoord.y/_ylim[0];
//    else if (graphCoord.x >= 0)
//        retVal.y = _center.y - ((double)(_center.y))*graphCoord.y/_ylim[1];

    return retVal;
}
