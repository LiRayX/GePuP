## 0. Envirment

The following packages can be found in lib/

+ Eigen
+ unsupported/Eigen/Splines
+ matplotlib-cpp: dependency of python matplotlib, numpy>=1.24, but numpy 2.x can't work

## 1. Function of Each File

### Basis

####  Vec.h

Class of point/vector in R^2.

#### Grid.h

Structured Grid.

#### Function.h

Class of scalar and vector function

#### MultiIndexSet.h

Operations on unordered_set of MultiIndex.

#### ButcherTable.h

ButcherTable of ClassicalRK4 and ARK436L2SA


#### RungeKutta.h

Time integral for the space semi-discrete system.

#### numlib.h

Numerical integral

#### Plot.h

Plot by matplotlib-cpp.

### Case 1: Free Boundary 

#### BoundaryCurve.h

Fit the boundary control points to smooth curve by piece-wise B-spline.  

#### CurveBelonging.h

For each spline, calculating which part of the curve belongs to which cell, intersections to face will be stored.

#### CellClassifier.h

Classifying all cells to dead cell, alive cell (cut cell, side cell, edge cell ,core cell).


### Case 2: Free Fall of A Ball

#### CyclicCurve.h

Describe the ball in 2D.

#### CellDivision.h

Dividing the cell into different type, according to the releative position to the ball.

#### CellHandler.h

Computering volume and centroid of cut cell.

## 2.Idea
+ 4th order operators coefficients can be writen to a funtion by swith case -2:2
+ Ghost cell cofficients can be treated in the same way, but swith case 0:4




