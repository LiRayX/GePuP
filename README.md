## 0. Envirment

This package can be found in lib/

+ Eigen
+ unsupported/Eigen/Splines
+ matplotlib-cpp: dependency of python matplotlib, numpy>=1.24, but numpy 2.x can't work

## 1. Function of Each File

####  Vec.h

Class of point/vector in R^2.

#### Grid.h

Structured Grid.

#### BoundaryCurve.h

Fit the boundary control points to smooth curve by piece-wise B-spline.  

#### CurveBelonging.h

For each spline, calculating which part of the curve belongs to which cell, intersections to face will be stored.

#### CellClassifier.h

Classifying all cells to dead cell, alive cell (cut cell, side cell, edge cell ,core cell).

#### ButcherTable.h

ButcherTable of ClassicalRK4 and ARK436L2SA

#### RungeKutta.h

Time integral for the space semi-discrete system.

#### numlib.h

 Numerical integral