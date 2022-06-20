#include "ecc.cpp"

void printPoint(ecc::Point p)
{
    std::cout << "(" << p.x << ", " << p.y << ")\n";
}       

int main(int argc, char** argv)
{
    ecc::init();
    std::cout.precision(100);
    std::cout << "(" << ecc::params.x << ", " << ecc::params.y << ")\n";

    // std::pair<DOUBLE, ecc::Point> k1 = ecc::deriveKeys(), k2 = ecc::deriveKeys();
    // DOUBLE a1 = k1.first, b1 = k2.first;
    // ecc::Point a2 = k1.second, b2 = k2.second;
    // ecc::Point s1 = ecc::mult(a2, b1), s2 = ecc::mult(b2, a1);
    // std::cout << s1.x << "\n";
    // std::cout << s2.x << "\n";

    return 0;
}