#include "../src/Vec.h"
#include "../lib/Eigen/Core"
#include <iostream>


int main()
{
    Eigen::Vector2d eigenv;
    eigenv << -1.0,2.0;
    Vec v(eigenv);
    std::cout<<"constructor from Eigen::Vector" << v << std::endl;
    Vec normal = outernormal(v);
    std::cout<< "noraml" << normal << std::endl;
}