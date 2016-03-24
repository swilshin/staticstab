'''
Calculate the difference in stability margin as a function of the fraction 
of the gait cycle for a quadrupedal walker, with varying lateral separation 
of the feet.

@author: Simon Wilshin
@contact: swilshin@rvc.ac.uk
@date: Mar 2016
'''

from staticstab import Polyped

from numpy import array,ones,pi,linspace
from pylab import plot,scatter,figure,show,cm,savefig,xlabel,ylabel,title

from os.path import join as osjoin

Nsubcyc = 500 # Number of subcycles to calculate
Nlatsep = 10 # Number of lateral seperations to calculate
amax = 2.0
amin = 0.1
tau = linspace(0,2*pi,Nsubcyc)
  # Define polyped parameters
g = 2.0*pi*ones((3))/4.0 # Relative phase
d = 2.0*pi*ones((4))*(3.0/4.0) # Duty cycle

for a in linspace(amin,amax,Nlatsep):
  x = array([[-0.4,a],[1.2,-a],[-0.8,-a],[1,a]]) # Touchdown positions
  lda = 0.3 # Stride length
  
  # Create a forward and reverse walker
  p0 = Polyped(g,d,x,lda)
  p1 = Polyped(-g,d,x,lda)
    
  plot(tau,array([p0.margin(t)-p1.margin(t) for t in tau]),c=cm.hot((a-amin)/(amax-amin)))

xlabel("fraction of stride")
ylabel("stability margin difference")
title("stability margin differences between forward and reverse walker\n(hotter colour is wider body)")
savefig(osjoin("figure","quadconfigs.png"))