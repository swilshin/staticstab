'''
Demonstration of the Polyped system for a spider like animal. Finds stance 
parameters at midstance and plots the foot positions. Has slightly 
odd middle feet to demonstrate that the convex hull is indeed found.

LICENSE
=======
  
  This file is part of Simon Wilshin's staticstab module.

  Simon Wilshin's staticstab module is free software: you can redistribute 
  it and/or modify it under the terms of the GNU General Public License as 
  published by the Free Software Foundation, either version 3 of the License, 
  or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.  


@author: Simon Wilshin
@contact: swilshin@rvc.ac.uk
@date: Jan 2015
'''

from staticstab import Polyped

from numpy import array,ones,pi
from pylab import plot,scatter,figure,savefig,axis,tight_layout

from os.path import join as osjoin

# Define polyped parameters
g = 2.0*pi*ones((7))/8.0 # Relative phase
d = 2.0*pi*ones((8))*(7.0/8.0) # Duty cycle
x = array([[-1,1],[1,-1],[-1,-1],[1,1],[-0.9,0.5],[1.1,-0.5],[-1.1,-0.5],[0.9,0.5]]) # Touchdown positions
lda = 0.3# Stride length

# Construct polyped
p = Polyped(g,d,x,lda)

# Do some calculations
print p.footPositions(pi) # Foot positions mid-cycle
print p.support(pi) # Polygon of support at mid-cycle
print p.margin(pi) # Stability margin at mid-cycle

# Scatter plot of foot position and polygon of support
figure()
plot(*p.support(pi).boundary.coords.xy,color='b',lw=4)
scatter(*array(p.footPositions(pi)).T,s=900,c='g')
axis('off')
tight_layout()
savefig(osjoin("figure","spiderfootplacements.png"))
savefig(osjoin("figure","spiderfootplacements.pdf"))
savefig(osjoin("figure","spiderfootplacements.svg"))