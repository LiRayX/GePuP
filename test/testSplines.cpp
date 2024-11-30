

#include <Eigen/Core>
#include <iostream>

int main(int argc, char const* argv[]) {
    Eigen::VectorXd xvals(3);
    Eigen::VectorXd yvals(xvals.rows());
    xvals << 0, 15, 30;
    yvals << 0, 12, 17;

    typedef Eigen::Spline<double, 1> SplineType;
    SplineType spline = Eigen::SplineFitting<SplineType>::Interpolate(yvals.transpose());

    std::cout << spline(12.34) << std::endl;
}