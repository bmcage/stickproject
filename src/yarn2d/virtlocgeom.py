#
# Copyright (C) 2010  B. Malengier
# Copyright (C) 2010  P.Li
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

""" 
Module with precomputed functions needed for virt loc computation 
"""
from __future__ import division
#-------------------------------------------------------------------------
#
# Global Imports
#
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#
# Intersection points of small circle with large circle
#
#-------------------------------------------------------------------------

def intersect_circles(a,b,r,s):
    """ Intersection points of the circles 
    sympy.solve((x**2 + y**2 - r**2, 
                     (x - a)**2 + (y - b)**2 - s**2), x, y)
    @return: a list with two points [(x1,y1), (x2,y2)]
    """
    return [
(-(-2*b**3*(b**6/(2*a**2*b**2 + a**4 + b**4) + 4*a**4/(-4*a**2 - 4*b**2) 
+ 4*b**4/(-4*a**2 - 4*b**2) + 4*r**4/(-4*a**2 - 4*b**2) + 4*s**4/(-4*a**2 - 4*b**2) 
+ a**4*b**2/(2*a**2*b**2 + a**4 + b**4) + b**2*r**4/(2*a**2*b**2 + a**4 + b**4) 
+ b**2*s**4/(2*a**2*b**2 + a**4 + b**4) - 8*a**2*r**2/(-4*a**2 - 4*b**2) 
- 8*a**2*s**2/(-4*a**2 - 4*b**2) - 8*b**2*s**2/(-4*a**2 - 4*b**2) - 8*r**2*s**2/(-4*a**2 - 4*b**2) 
- 2*b**4*s**2/(2*a**2*b**2 + a**4 + b**4) + 2*a**2*b**4/(2*a**2*b**2 + a**4 + b**4) 
+ 2*b**4*r**2/(2*a**2*b**2 + a**4 + b**4) + 8*a**2*b**2/(-4*a**2 - 4*b**2) 
+ 8*b**2*r**2/(-4*a**2 - 4*b**2) - 2*a**2*b**2*s**2/(2*a**2*b**2 + a**4 + b**4) 
- 2*b**2*r**2*s**2/(2*a**2*b**2 + a**4 + b**4) + 2*a**2*b**2*r**2/(2*a**2*b**2 + a**4 + b**4))**(1/2) 
- 2*b*a**2*(b**6/(2*a**2*b**2 + a**4 + b**4) + 4*a**4/(-4*a**2 - 4*b**2) + 4*b**4/(-4*a**2 - 4*b**2) 
+ 4*r**4/(-4*a**2 - 4*b**2) + 4*s**4/(-4*a**2 - 4*b**2) + a**4*b**2/(2*a**2*b**2 + a**4 + b**4) 
+ b**2*r**4/(2*a**2*b**2 + a**4 + b**4) + b**2*s**4/(2*a**2*b**2 + a**4 + b**4) 
- 8*a**2*r**2/(-4*a**2 - 4*b**2) - 8*a**2*s**2/(-4*a**2 - 4*b**2) - 8*b**2*s**2/(-4*a**2 - 4*b**2) 
- 8*r**2*s**2/(-4*a**2 - 4*b**2) - 2*b**4*s**2/(2*a**2*b**2 + a**4 + b**4) 
+ 2*a**2*b**4/(2*a**2*b**2 + a**4 + b**4) + 2*b**4*r**2/(2*a**2*b**2 + a**4 + b**4) 
+ 8*a**2*b**2/(-4*a**2 - 4*b**2) + 8*b**2*r**2/(-4*a**2 - 4*b**2) - 2*a**2*b**2*s**2/(2*a**2*b**2 
+ a**4 + b**4) - 2*b**2*r**2*s**2/(2*a**2*b**2 + a**4 + b**4) 
+ 2*a**2*b**2*r**2/(2*a**2*b**2 + a**4 + b**4))**(1/2) 
- 2*a**2*b**2 - 2*b**2*r**2 + 2*b**2*s**2 - 2*b**4)/(-4*a*b**2 - 4*a**3) 
- (s**2 - a**2 - b**2 - r**2)/(2*a),

 -(b*a**2 + b*r**2 - b*s**2 + b**3)/(2*(-a**2 - b**2)) 
+ (-4*(-2*a**2*b**2 - 2*b**2*r**2 + 2*a**2*r**2 + 2*a**2*s**2 + 2*b**2*s**2 + 2*r**2*s**2 
- a**4 - b**4 - r**4 - s**4)/(-4*a**2 - 4*b**2) 
+ (b*a**2 + b*r**2 - b*s**2 + b**3)**2/(-a**2 - b**2)**2)**(1/2)/2),

 (-(2*b**3*(b**6/(2*a**2*b**2 + a**4 + b**4) + 4*a**4/(-4*a**2 - 4*b**2) 
+ 4*b**4/(-4*a**2 - 4*b**2) + 4*r**4/(-4*a**2 - 4*b**2) + 4*s**4/(-4*a**2 - 4*b**2) 
+ a**4*b**2/(2*a**2*b**2 + a**4 + b**4) + b**2*r**4/(2*a**2*b**2 + a**4 + b**4) 
+ b**2*s**4/(2*a**2*b**2 + a**4 + b**4) - 8*a**2*r**2/(-4*a**2 - 4*b**2) 
- 8*a**2*s**2/(-4*a**2 - 4*b**2) - 8*b**2*s**2/(-4*a**2 - 4*b**2) - 8*r**2*s**2/(-4*a**2 - 4*b**2) 
- 2*b**4*s**2/(2*a**2*b**2 + a**4 + b**4) + 2*a**2*b**4/(2*a**2*b**2 + a**4 + b**4) 
+ 2*b**4*r**2/(2*a**2*b**2 + a**4 + b**4) + 8*a**2*b**2/(-4*a**2 - 4*b**2) 
+ 8*b**2*r**2/(-4*a**2 - 4*b**2) - 2*a**2*b**2*s**2/(2*a**2*b**2 + a**4 + b**4) 
- 2*b**2*r**2*s**2/(2*a**2*b**2 + a**4 + b**4) + 2*a**2*b**2*r**2/(2*a**2*b**2 + a**4 + b**4))**(1/2) 
+ 2*b*a**2*(b**6/(2*a**2*b**2 + a**4 + b**4) + 4*a**4/(-4*a**2 - 4*b**2) + 4*b**4/(-4*a**2 - 4*b**2) 
+ 4*r**4/(-4*a**2 - 4*b**2) + 4*s**4/(-4*a**2 - 4*b**2) + a**4*b**2/(2*a**2*b**2 + a**4 + b**4)
 + b**2*r**4/(2*a**2*b**2 + a**4 + b**4) + b**2*s**4/(2*a**2*b**2 + a**4 + b**4) 
- 8*a**2*r**2/(-4*a**2 - 4*b**2) - 8*a**2*s**2/(-4*a**2 - 4*b**2) - 8*b**2*s**2/(-4*a**2 - 4*b**2) 
- 8*r**2*s**2/(-4*a**2 - 4*b**2) - 2*b**4*s**2/(2*a**2*b**2 + a**4 + b**4) 
+ 2*a**2*b**4/(2*a**2*b**2 + a**4 + b**4) + 2*b**4*r**2/(2*a**2*b**2 + a**4 + b**4) 
+ 8*a**2*b**2/(-4*a**2 - 4*b**2) + 8*b**2*r**2/(-4*a**2 - 4*b**2) 
- 2*a**2*b**2*s**2/(2*a**2*b**2 + a**4 + b**4) 
- 2*b**2*r**2*s**2/(2*a**2*b**2 + a**4 + b**4) + 2*a**2*b**2*r**2/(2*a**2*b**2 + a**4 + b**4))**(1/2) 
- 2*a**2*b**2 - 2*b**2*r**2 + 2*b**2*s**2 - 2*b**4)/(-4*a*b**2 - 4*a**3) 
- (s**2 - a**2 - b**2 - r**2)/(2*a), 

-(b*a**2 + b*r**2 - b*s**2 + b**3)/(2*(-a**2 - b**2)) 
- (-4*(-2*a**2*b**2 - 2*b**2*r**2 + 2*a**2*r**2 + 2*a**2*s**2 + 2*b**2*s**2 
+ 2*r**2*s**2 - a**4 - b**4 - r**4 - s**4)/(-4*a**2 - 4*b**2) 
+ (b*a**2 + b*r**2 - b*s**2 + b**3)**2/(-a**2 - b**2)**2)**(1/2)/2)]


#solve((x**2+y**2-0.2**2, (x+)**2+(y-0.133749497)**2-0.07142857**2),x,y)
print intersect_circles(-0.35810915869,-0.0213824327636,0.4,0.07142857)