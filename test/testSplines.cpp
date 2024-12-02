#define EIGEN_MALLOC_ALREADY_ALIGNED 0
#include "../lib/Eigen/Core"
#include "../lib/unsupported/Eigen/Splines"
#include "../src/Vec.h"
#include <iostream>

typedef Eigen::Spline<double, 2, 3> Spline2d;

int main() 
{
    // 给定一组2D点
    std::vector<Vec> points;
    points.push_back(Vec{0, 0});
    points.push_back(Vec{1, 2});
    // points.push_back(Vec{1, 2});
    // points.push_back(Vec{2, 3});
    points.push_back(Vec{3, 1});
    points.push_back(Vec{4, 0});
    // points.push_back(Eigen::Vector2d(1, 2));
    // points.push_back(Eigen::Vector2d(2, 3));
    // points.push_back(Eigen::Vector2d(3, 1));
    // points.push_back(Eigen::Vector2d(4, 0));

    // 将点转换为Eigen::MatrixXd格式
    Eigen::MatrixXd pts(2, points.size());
    for (int i = 0; i < points.size(); ++i) {
        pts(0, i) = points[i][0];
        pts(1, i) = points[i][1];
    }

    // 创建样条曲线
    Spline2d spline = Eigen::SplineFitting<Spline2d>::Interpolate(pts, 3);

    // 在参数u = 0到1之间均匀采样点
    for (double u = 0; u <= 1; u += 0.1) {
        Eigen::Vector2d point = spline(u);
        std::cout << "u = " << u << ", point = (" << point(0) << ", " << point(1) << ")" << std::endl;
    }

     double t = 0.5; // 参数t的取值范围为[0, 1]
    Eigen::Vector2d derivative = spline.derivatives(t, 1).col(1);

    // 计算法向量
    Eigen::Vector2d normal(-derivative[1], derivative[0]);
    normal.normalize();

    // 输出结果
    std::cout << "导数: " << derivative.transpose() << std::endl;
    std::cout << "法向量: " << normal.transpose() << std::endl;

    

    return 0;
}