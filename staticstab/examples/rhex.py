'''
File rhexwork.py

Calculating the stability margin of the RHex robotic system and plotting the 
behaviour

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
@date: Nov 2014
'''

from staticstab import Polyped

from os.path import join as osjoin

from numpy import pi,array,sum,abs,linspace
from pylab import (figure, scatter,plot,clf,axis,tight_layout,savefig,legend,
  xlabel,xticks,xlim,ylabel,xlabel)

# Directories to save plots in
animDir = "anim"
figDir = "figure"

###############################
# Define our legged system
###############################
# Relative phases
g0 = 2.0*pi*array([1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0])
# Duty factors
d0 = 2.0*pi*array([5.0/6.0,5.0/6.0,5.0/6.0,5.0/6.0,5.0/6.0,5.0/6.0])
# Normalisation 40cm = 1.0 u
normL = 40
# Dimensionless stride lengths
lda0 = 0.50667 # 1.8288m in nine strides
# Dimensionless positions
l = 0.5
a = 0.33125
b = 0.43125
s = -lda0/(2.0*2.0)#-lda0/(2.0*6.0)
x0 = array([[-l+s,a],[s,-b],[l+s,a],[-l+s,-a],[s,b],[l+s,-a]])
# McGhee and Frank style parameter sets
k0 = (g0,d0,x0,lda0)
k1 = (-g0,d0,x0,lda0)

# Simulation parameters
# Number of frames to interpolate the cycle by
Ninterp = 240 # 120,240
# Number of strides to animate
Nstride = 5 #1,5

# Create walkers
p0 = Polyped(*k0)
p1 = Polyped(*k1)

M = list()
tau = linspace(0,Nstride*2*pi,Ninterp*Nstride)
for i,tm in enumerate(tau):
  t = tm %(2*pi) # Stride fractions must be between 0 and 2*pi
  
  # Calculate stability margins
  M.append([p0.margin(t),p1.margin(t)])

  # Plot the convex hull of this pattern of motion for the forward walk
  figure(1)
  clf()
  plot(*p0.support(t).boundary.coords.xy,color='b',lw=4)
  #cv = (t-p0.fon)/((p0.foff-p0.fon)%(2*pi))%1.0
  scatter(*array(p0.support(t).boundary.coords.xy)[:,:-1],s=900,c='m')
  axis([-1.5,1.5,-1.5,1.5])
  axis('off')
  tight_layout()
  savefig(osjoin(animDir,"forwardgait"+str(i).zfill(5)))
  
  # Plot the convex hull of this pattern of motion for the reverse walk
  figure(2)
  clf()
  plot(*p1.support(t).boundary.coords.xy,color='b',lw=4)
  scatter(*array(p1.support(t).boundary.coords.xy)[:,:-1],s=900,c='r')
  axis([-1.5,1.5,-1.5,1.5])
  axis('off')
  tight_layout()
  savefig(osjoin(animDir,"reversegait"+str(i).zfill(5)))

M = array(M)

# Plot the quasi-static stability margin for the forward and reverse walker 
# as a function of position in the cycle.
figure(3)
clf()
plot(tau,normL*M[:,0],c='m',label="forward",lw=3)
plot(tau,normL*M[:,1],c='r',label="reverse",lw=3,ls='--')
legend()
xlabel("phase of cycle (rad)")
xticks(
  [0,pi/2.0,pi,3.0*pi/2.0,2*pi],
  ["$0$","$\\pi/2$","$\\pi$","$3 \\pi/2$","$2\\pi$"]
)
xlim([0,2*pi])
ylabel("static stability margin (cm)")
savefig(osjoin(figDir,"rhexstabmargins.png"))
savefig(osjoin(figDir,"rhexstabmargins.pdf"))
savefig(osjoin(figDir,"rhexstabmargins.svg"))

