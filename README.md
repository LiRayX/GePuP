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

Numerical integral, different version for scalar function and vector function

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

#### CurvedTriangle

To computer the volume and centriod of the alive region of cut-cell

+ General triangle
+ Triangle with one curved edge, given by part of CyclicCurve
+ Quadrilateral with one curved edge, dividing into a curved triangle and a general triangle


#### CellHandler.h

Computering volume and centroid of cut cell.




