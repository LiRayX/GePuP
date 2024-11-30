#include "../src/Grid.h"
#include "../unsupported/Eigen/Splines"

int main()
{
  Vec low{0, 0};
  Vec high{1, 1};
  Vec pos{0.249, 0.25};
  Grid g(low, high, 0.25);
  std::cout<<"grid size: "<<g.get_size()[0]<<" "<<g.get_size()[1]<<std::endl;
  Grid gd = g;
  std::cout << gd;
  Vec index = gd(1, 2);
  Vec c = gd.center(1, 2);
  std::cout << index << " " << c << std::endl;
  std::cout <<"simple index"<<std::endl;
  int testindex = gd.MultiToSingle(1, 2);
  std::cout <<"testindex" <<testindex << std::endl;
  MultiIndex testmulti = gd.LocateCell(pos);
    std::cout << testmulti[0] << " " << testmulti[1] << std::endl;
  // int testindex = gd.CellIndex(1, 2);
  // std::cout << testindex << std::endl;
  // std::cout << gd.IfOnLeftBound(1, 2) << std::endl;
  // std::cout << gd.IfOnRightBound(3, 2) << std::endl;

}