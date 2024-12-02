#define EIGEN_MALLOC_ALREADY_ALIGNED 0
#include "../lib/Eigen/Core"
#include "../lib/unsupported/Eigen/Splines"
#include "../src/Vec.h"
#include <iostream>

typedef Eigen::Spline<double, 2, 3> Spline2d;

//插值得到的Spline,参数是均匀分布的


int main()
{
     // 定义四个点
    Eigen::MatrixXd points(2, 4);
    points << 0, 1, 2, 3,
              0, 1, 0, -1;

    // 定义参数
    Eigen::VectorXd parameters(4);
    parameters << 0, 1.0/3, 2.0/3, 1;

    // 创建样条插值对象
    // Eigen::Spline<double, 2> spline = Eigen::SplineFitting<Eigen::Spline<double, 2>>::Interpolate(points, 3, parameters);
    Eigen::Spline<double, 2> spline = Eigen::SplineFitting<Eigen::Spline<double, 2>>::Interpolate(points, 3);

    // 输出插值结果
    for (double t = 0; t <= 1; t += 0.1) {
        Eigen::Vector2d point = spline(t);
        std::cout << "t = " << t << ", point = (" << point(0) << ", " << point(1) << ")" << std::endl;
    }
    Eigen::Vector2d point = spline(1.0/3);
    std::cout<<"t="<<1.0/3<<", point = (" << point(0) << ", " << point(1) << ")"<<std::endl;

    return 0;
}