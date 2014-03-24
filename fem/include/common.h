#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define EPS_0       (8.854187817e-12)
#define CALC_EPS    (1e-15)
#define LIGHT_SPEED (299792458)
#define MU0         (4*M_PI*1e-7)
#define EEC         (1.60217646e-19)
#define BUFF_SIZE   255

using namespace std;

#include "classes.h"

ostream& operator<<(ostream& os,const Node& obj);
ostream& operator<<(ostream& os,/*const*/ Element& obj);
ostream& operator<<(ostream& os,/*const*/ Edge& obj);
#endif
