# ocvGraph (V1.3.0)
Implementation of class for graph plotting in C++ through OpenCV 3.

Features:
+Set arbitrary range on X & Y axes
+Plot circles, lines(in Cartesian or polar), polynomials of any degree
+Plot in different colors
+Append legend describing which line represents which function

Current implementation developed and tested under OpenCV 3.2.0

For compiling sources in /examples from Linux terminal you can use the command below while in folder /examples:

`` g++ -std=c++11  main.cpp ../ocvgraph.cpp -o main `pkg-config --cflags --libs opencv-3.2.0-dev` ``

NOTE: Replace 'opencv-3.2.0-dev' with your version of OpenCV.


Results from example/main.cpp program:

![alt tag](https://hsr.duckdns.org/images/polarRadar.png)

![alt tag](https://hsr.duckdns.org/images/parabola.png)

![alt tag](https://hsr.duckdns.org/images/lines.png)

Change-log:

 *  V1.0.0 - 1.6.2017
   +Custom origin of coordinate system at any point of image. Support for
 *  plotting circles, lines (in Cartesian and polar coordinate system) and
 *  writing text.
 *  V1.1.0
   +Save plot to a given path
 *  V1.2.0 - 17.6.2017
   +Added clear function to reset the graph
   +Added option to plot N-degree polynomials in Cartesian, within given range
   +Show data axes, variable number of divisions for X and Y axes
 *  V1.3.0 -20.6.2017
   +Added macros for common colors
   +Can print legend based on the color of a function
   +Added interface for setting X and Y intervals for axes
   +Expanded PolyN function by adding more options for customization
