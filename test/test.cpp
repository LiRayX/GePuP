#include "../src/Grid.h"

int main()
{
  Vec low{0, 0};
  Vec high{1, 1};
  Grid g(low, high, 0.25);
  std::cout<<"grid size: "<<g.get_size()[0]<<" "<<g.get_size()[1]<<std::endl;

  MultiIndex index{1,1};
  Normal normal{1,0};
  Vec center = g.center(index[0], index[1]);
  std::cout<<"center: "<<center<<std::endl;
  Vec pos{0.5, 0.35};
  Vec pos1{0.25, 0.5};
  Vec pos2{0.5, 0.75};
  bool onface = g.OnFace(index, normal, pos, 1e-6);
  // bool onface1 = g.OnFace(1, 1, 1, 0, pos1, 1e-6);
  // bool onface2 = g.OnFace(1, 1, 1, 0, pos2, 1e-6);
  std::cout<<"onface: "<<onface<<std::endl;
  // std::cout<<"onface1: "<<onface1<<std::endl;
  // std::cout<<"onface2: "<<onface2<<std::endl;

  // std::cout << index << " " << c << std::endl;
  // std::cout <<"simple index"<<std::endl;
  // MultiIndex testmulti = gd.LocateCell(pos);
  // std::cout << testmulti[0] << " " << testmulti[1] << std::endl;
  // int testindex = gd.CellIndex(1, 2);
  // std::cout << testindex << std::endl;
  // std::cout << gd.IfOnLeftBound(1, 2) << std::endl;
  // std::cout << gd.IfOnRightBound(3, 2) << std::endl;

}