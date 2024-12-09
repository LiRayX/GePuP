#include <unordered_map>
#include <vector>
#include <iostream>
// 假设 MultiIndex 和 CutCellInfo 已经定义好
using MultiIndex = int; // 示例中使用 int 作为 MultiIndex
struct ParaInterval {
    // 假设 ParaInterval 已经定义好
};
struct CutCellInfo {
    std::vector<int> index_spline;
    std::vector<ParaInterval> parainterval;
};

using CutCellMap = std::unordered_map<MultiIndex, CutCellInfo>;

int main() {
    CutCellMap cutCellMap;

    // 创建一个 MultiIndex
    MultiIndex index = 1;

    // 如果键不存在，operator[] 会创建一个默认的 CutCellInfo 对象
    cutCellMap[index].index_spline.push_back(1);
    cutCellMap[index].index_spline.push_back(2);
    cutCellMap[index].index_spline.push_back(3);

    // 访问并打印 index_spline 的内容
    for (int value : cutCellMap[index].index_spline) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return 0;
}