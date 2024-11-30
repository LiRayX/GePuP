#pragma once

#include "Vec.h"
#include "Grid.h"
#include <iostream>
#include <vector>

using VecList = std::vector<Vec>;

class BoundaryCurve
{
public:
    BoundaryCurve() = default;
    BoundaryCurve(const VecList& controlPoints, const std::vector<int>& endPoints)
        : ControlPoints(controlPoints), EndPoints(endPoints) {}

    VecList GetControlPoints() const { return ControlPoints; }
    std::vector<int> GetEndPoints() const { return EndPoints; }


protected:
    VecList ControlPoints;
    std::vector<int> EndPoints;
};