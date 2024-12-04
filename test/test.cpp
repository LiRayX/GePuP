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
  Vec pos1{0.49, 0.35};
  Vec pos2{0.5-1e-7, 0.35};

  bool onface = g.OnFace(index, normal, pos, 1e-6);
  double dis = g.SignDistance(index, normal, pos);
  double dis1 = g.SignDistance(index, normal, pos1);
  double dis2 = g.SignDistance(index, normal, pos2);
  // bool onface1 = g.OnFace(1, 1, 1, 0, pos1, 1e-6);
  // bool onface2 = g.OnFace(1, 1, 1, 0, pos2, 1e-6);
  std::cout<<"onface: "<<onface<<std::endl;
  std::cout<<"dis:"<<dis<<std::endl;
  std::cout<<"dis1:"<<dis1<<std::endl;
  std::cout<<"dis2:"<<dis2<<std::endl;
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