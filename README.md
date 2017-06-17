# ocvGraph
Basic implementation of graph plotting in C++ through OpenCV

Current implementation developed and tested under OpenCV 3.2.0

For compiling sources in /examples from terminal you can use the command below from terminal while in folder /examples

`` g++ -std=c++11  main.cpp ../ocvgraph.cpp -o main `pkg-config --cflags --libs opencv-3.2.0-dev` ``

NOTE: Replace 'opencv-3.2.0-dev' with your version of OpenCV.


Results from example/main.cpp program:

![alt tag](https://hsr.duckdns.org/images/polarRadar.png)

![alt tag](https://hsr.duckdns.org/images/parabola.png)

![alt tag](https://hsr.duckdns.org/images/lines.png)
