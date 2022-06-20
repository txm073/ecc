/*
C++ source module implementing Elliptic Curve Cryptography
*/
#ifndef _ECC_CPP_
#define _ECC_CPP_

#ifdef _NULL_
    #define ASSISTANT_MODULE_EXPORT
    #define MODULE_RELPATH "security/crypto"
    #define MODULE_NAME "ecc"
#endif

#include "ecc.h"


namespace ecc {

// Calculate the highest common factor of two integers
DOUBLE hcf(DOUBLE i, DOUBLE j)
{
    DOUBLE temp;
    while (true) {
        temp = fmod(i, j);
        if (temp == 0) {
            return j;
        }
        i = j;
        j = temp;
    }
}

// Calculate the multiplicative inverse of a under modulo m
DOUBLE modInverse(DOUBLE a, DOUBLE m)
{
    DOUBLE m0 = m;
    DOUBLE y = 0, x = 1;
 
    if (m == 1)
        return 0;
 
    while (a > 1) {
        DOUBLE q = a / m;
        DOUBLE t = m;
 
        m = fmod(a, m), a = t;
        t = y;
 
        y = x - q * y;
        x = t;
    }
 
    return fmod(x, m0);
}
    
// Get a random integer within a certain range
DOUBLE randRange(DOUBLE min, DOUBLE max) 
{
   std::uniform_real_distribution<DOUBLE> unif(min, max);
   std::default_random_engine re;
   return unif(re);
}

// Initialise parameters
void init() 
{
    char** endPtr;
    params.p = std::pow(2, 256) - std::pow(2, 32) - std::pow(2, 9) - std::pow(2, 8) - std::pow(2, 7) - std::pow(2, 6) - std::pow(2, 4) - 1;
    params.a = 0;
    params.b = 7;
    params.x = std::stold("55066263022277343669578718895168534326250603453777594175500187360389116729240");
    params.n = std::stold("115792089237316195423570985008687907852837564279074904382605163141518161494337");
    params.h = 1;
   
    std::string yStr = "32670510020758816978083085130507043184471273380659243275938904335757337482424";
    for (int i = 0; i < yStr.size(); ++i) {
        params.y += ((char(yStr[yStr.size() - i - 1]) - 48) * pow(10, i));
    }
}

// Calculate the addition of two points in a finite field
Point add(Point p, Point q) 
{  
    DOUBLE m, x, y;
    Point result;
    // Check for vertical line
    if (p.x == q.x && p.y != q.y) {
        result.valid = false;
        return result;
    }
    // Division operation x / y defined as x * modinv(y, params.p)
    
    if (p.x == q.x) {
        // Gradient calculated with implicit differentiation
        m = (3 * p.x * p.x + params.a) * (modInverse(2 * p.y, params.p));
    } else {
        // Gradient calculated with Δy / Δx
        m = (p.y - q.y) * modInverse(p.x - q.x, params.p);
    }
    x = m * m - p.x - q.x;
    y = p.y + m * (x - p.x);

    result.x = fmod(x, params.p);
    result.y = fmod(-1 * y, params.p);
    result.valid = true;
    return result;
}

// Multiply a point on a finite field by n
Point mult(Point p, DOUBLE n) 
{
    Point result, point;
    if (!n) {
        result.valid = false;
        return result;
    }
    point = p;
    result = add(point, point);
    while (n) {
        if (std::fmod(n, 1)) {
            result = add(result, point);       
        }
        point = add(point, point);
        n = std::floor(n / 2);
    }
    return result;
}

std::pair<DOUBLE, Point> deriveKeys()
{
    Point gen;
    gen.x = params.x;
    gen.y = params.y;
    DOUBLE privateKey = randRange(1, params.n);
    Point publicKey = mult(gen, privateKey);
    return std::pair<DOUBLE, Point>(privateKey, publicKey);
}

}

#endif