//geometry.h -- the simulation box as sph.cpp and rocks.h both see it: the
//periodic x direction and the smooth wall repulsion. Self-contained.
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>

/******Periodic boundary in x******/
//When true, the x walls are removed: a particle leaving through x = x0+Lx
//re-enters at x0 and vice versa, and every pair separation uses the minimum
//image in x. Requires Lx >= 3*h so the neighbour grid has at least 3 distinct
//cells across x (otherwise the 3x3x3 stencil would visit a cell twice).
const bool periodic_x = true;
double box_x0 = 0.0, box_Lx = 0.0; //set by boundaries::set_boundaries

//map x back into [box_x0, box_x0+box_Lx)
inline double wrap_x(double x){
  if(!periodic_x) return x;
  x -= box_Lx*floor( (x-box_x0)/box_Lx );
  if(x >= box_x0+box_Lx) x = box_x0;//floor rounding can land exactly on the upper edge
  return x;
}

//minimum-image x separation: the nearest periodic copy of the pair
inline double min_image_x(double dx){
  if(!periodic_x) return dx;
  return dx - box_Lx*round(dx/box_Lx);
}

/******Wall repulsion******/
//Monaghan-style boundary repulsion, bounded (no singularity): ramps from 0 at
//d=r0 up to D at d=0, so it stays numerically safe under a fixed dt.
//Past the wall (d<=0) it holds at D rather than dropping back to 0 -- returning
//0 there left a particle that overshot with nothing pushing it back in.
//Used for the box walls (sph.cpp) and the rock faces (rocks.h); the range r0
//and strength D are the boundary_r0/boundary_D parameters in sph.cpp.
inline double wall_repulsion(double d, double r0, double D){
  if(d >= r0) return 0.0;
  if(d <= 0.0) return D;
  double x = 1.0 - d/r0;
  return D*x*x;
}

#endif
