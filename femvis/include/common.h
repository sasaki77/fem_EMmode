#ifndef COMMON
#define COMMON


#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <time.h>

#define EPS_0 (8.854187817e-12)
#define CALC_EPS (1e-15)
#define LIGHT_SPEED (299792458)
#define MU0 (4*M_PI*1e-7)
#define EEC (1.60217646e-19)
using namespace std;

class Node;
class Element;
ostream& operator<<(std::ostream& os,const Node& obj);
ostream& operator<<(std::ostream& os,const Element& obj);

#include "classes.h"

#endif
