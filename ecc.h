#ifndef _ECC_H_
#define _ECC_H_

#include <utility>
#include <cmath>
#include <string>
#include <iostream>
#include <random>

#define DOUBLE long double

namespace ecc {

DOUBLE hcf(DOUBLE i, DOUBLE j);
DOUBLE modInverse(DOUBLE a, DOUBLE m);
DOUBLE randRange(DOUBLE min, DOUBLE max);

struct Params {
    DOUBLE p, a, b, x, y, n, h;
};

struct Point {
    DOUBLE x, y;
    bool valid;
};

static Params params;

void init();
Point add(Point p, Point q);
Point mult(Point p, DOUBLE n);
std::pair<DOUBLE, Point> deriveKeys();

};  

#endif
// 32670510020758816978083085130507043184471273380659243275938904335757337482424
// 32670510020758816978925480275470138904606275542450609660186942991101442654208