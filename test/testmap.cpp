#include <unordered_set>
#include <algorithm>
#include <iterator>
#include <iostream>

// 并集
template <typename T>
std::unordered_set<T> setUnion(const std::unordered_set<T>& set1, const std::unordered_set<T>& set2) {
    std::unordered_set<T> resultSet = set1;
    resultSet.insert(set2.begin(), set2.end());
    return resultSet;
}

// 交集
template <typename T>
std::unordered_set<T> setIntersection(const std::unordered_set<T>& set1, const std::unordered_set<T>& set2) {
    std::unordered_set<T> resultSet;
    for (const auto& elem : set1) {
        if (set2.find(elem) != set2.end()) {
            resultSet.insert(elem);
        }
    }
    return resultSet;
}

// 差集
template <typename T>
std::unordered_set<T> setDifference(const std::unordered_set<T>& set1, const std::unordered_set<T>& set2) {
    std::unordered_set<T> resultSet = set1;
    for (const auto& elem : set2) {
        resultSet.erase(elem);
    }
    return resultSet;
}

// 示例代码
int main() {
    std::unordered_set<int> set1 = {1, 2, 3, 4};
    std::unordered_set<int> set2 = {3, 4, 5, 6};

    auto unionSet = setUnion(set1, set2);
    auto intersectionSet = setIntersection(set1, set2);
    auto differenceSet = setDifference(set1, set2);

    std::cout << "Union: ";
    for (const auto& elem : unionSet) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;

    std::cout << "Intersection: ";
    for (const auto& elem : intersectionSet) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;

    std::cout << "Difference: ";
    for (const auto& elem : differenceSet) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;

    return 0;
}