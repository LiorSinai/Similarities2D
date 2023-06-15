# Ellipses

General algorithms relating to ellipses. See [en.wikipedia.org/wiki/Ellipse](https://en.wikipedia.org/wiki/Ellipse).

Ellipses are points which satisfy the equation: `(x'/a)^2 + (y'/b)^2 = 1`. 

An ellipse can be rotated and translated: `x = Rx' + t`. This leads to the general conic form:
```
a0 + a1*x + a2*y + a3*x^2 + a4*xy + a5*y^2 = 0
```
Where `a4^2 - 4a3*a5 < 0`.

## Algorithms

General ellipse formulas:
- Point in/on ellipse.
- Conic section coefficients (`a_i` in the conic form equation). 
- Arc midpoint between two angles.
- Segment area.

Between two ellipses:
- Intersection points.
- Areas of intersection. 

Reference: [Calculating ellipse overlap areas by Gary B. Hughes and Mohcine Chraibi (2011)](https://arxiv.org/abs/1106.3787).

## Test

```julia
include("src/ellipses.jl")
include("test/ellipse.jl")
```