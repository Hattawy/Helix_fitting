//#include <stdio.h>
#include <math.h>

/* return the points of intersection of two circles
   function value  returned is number of points of intersection
   The method is solve the two circle equations R^2= (x-xc)^2 + (y-yc)^2.
   Below, we solve first for the two values of y providing  a  solution.
   Each of those can yield two solutions for  x on each circle, but
   points on both circles can be identified and returned to the caller.*/

/*REVISION HISTORY
  9-FEB-07: hcf opened up 'small' from .001 to .01 to solve a precision problem 

  9-NOV-07: hcf implemented the following more elegant solution to eliminate
            the 'difference of two large numbers problem' I had originally
            coded-in to the solution. The following was found at 
            http://mathforum.org/library/drmath/view/51836.html
    
    Let the centers be: (a,b), (c,d)
    Let the radii be: r, s

     e = c - a                          [difference in x coordinates]
     f = d - b                          [difference in y coordinates]
     p = sqrt(e2 + f2)                [distance between centers]
     k = (p2 + r2 - s2)/(2p)         [distance from center 1 to line
					 joining points of intersection]
     x = a + ek/p + (f/p)sqrt(r2 - k2)
     y = b + fk/p - (e/p)sqrt(r2 - k2)
    OR
     x = a + ek/p - (f/p)sqrt(r2 - k2)
     y = b + fk/p + (e/p)sqrt(r2 - k2)
*/

int circle_int(float r1, float xc1, float yc1, float r2, float xc2,  
	       float yc2, float *xx1, float *yy1, float *xx2, float *yy2)
{
  double ee, ff, kk, mm, pp;
  /* Define a small floating-point number */
#define SMALL 0.0001

  /* Solve for two points of intersection using above formula: */
  ee= xc2 - xc1;
  ff= yc2 - yc1;
  pp= sqrt( ee*ee + ff*ff);  /* [distance between centers] */

  /* Verify that the circles do intersect! */
  if (pp > (r1+r2)) return (0); /* the two circles do not overlap at all */
  if (pp < (fabs(r1-r2))) return (0); /* small circle entirely inside large one */

  kk= (pp*pp + r1*r1 - r2*r2)/(2.0*pp);  
  mm= sqrt(r1*r1 - kk*kk);

  *xx1= xc1 + ee*kk/pp + (ff/pp)*mm;
  *yy1= yc1 + ff*kk/pp - (ee/pp)*mm;

  *xx2= xc1 + ee*kk/pp - (ff/pp)*mm;
  *yy2= yc1 + ff*kk/pp + (ee/pp)*mm;

  if ( (pow(*xx1-*xx2,2) + pow(*yy1-*yy2,2)) < SMALL)
    return(1); 
  else
    return(2); 
}
