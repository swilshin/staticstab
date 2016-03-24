'''
File staticstable.py

OVERVIEW
========

  Analyse the static stability of polypedal animal gaits. Primarily used by 
  creating a polyped object and then interrogating this object for the static 
  stability margin, polygon of support and other quantities at different parts 
  of the cycle.

MAIN CLASSES
============
  
  Foot - Created with a kinematic gait formula, this allows us to calculate 
  where feet will be during a cycle.
  Quadruped - A wrapper for four feet, can calculate the minimal longitudinal 
  stability of that four feet system.

TYPICAL USAGE
=============

  >>> # According to MG-F this walker is stable, let's check
  >>> beta = 11.0/12.0 # Duty cycle
  >>> delta = 1.0/4.0 # Lateral position of foot
  >>> k = ( \
    beta,beta,beta,beta, \
    1,1,0,0, \
    delta,-delta,delta,-delta, \
    0.0,1.0/2.0,3.0/4.0,1.0/4.0 \
   )
  >>> f = Foot(0.75,0.3,0.3,0.0)
  >>> print f.R(0.3)
  (0.0, 0.3)
  >>> q = Quadruped(k)
  >>> print q.R(0.5)
  [(0.5, 0.25), (1.0, -0.25), (-0.75, 0.25), (-0.25, -0.25)]

THEORY OF OPERATION
===================

  This class implements the model of McGhee and Frank, based on longitudinal 
  static stability. The gait is specified in the same way.
  
  All phases are stored in fraction of cycle.

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
  
REFERENCES
==========

  McGhee, Robert B., and Frank, Andrew A. "On the stability properties of 
  quadruped creeping gaits." Mathematical Biosciences 3 (1968): 331-351.

@author: Simon Wilshin
@contact: swilshin@rvc.ac.uk
@date: Mar 2016
'''

from numpy import isnan,nan

class Foot(object):
  '''
  A foot for a walker with duty factor beta, fore-aft touch down position
  gamma, lateral position delta and phase offset phi.

  @ivar beta: duty factor
  @type beta: float
  @ivar gamma: fore-aft positions
  @type gamma: float
  @ivar delta: lateral position
  @type delta: float
  @ivar phi: phase offset
  @type phi: float
  '''
  def __init__(self,beta,gamma,delta,phi):
    '''
    A foot for a walker, full documentation in L{Foot}.
    
    @param beta: duty factor
    @type beta: float
    @param gamma: fore-aft positions
    @type gamma: float
    @param delta: lateral position
    @type delta: float
    @param phi: phase offset
    @type phi: float
    '''
    self.beta = beta
    self.gamma = gamma
    self.delta = delta
    self.phi = phi
  
  def X(self,t):
    '''
    Calculate the fore-aft position of the foot at time t. If foot is not on 
    the ground this returns nan.
    
    @param t: time in cycle (fraction of cycle)
    @type t: float    
    @return: fore-aft position of foot
    @rtype: float
    '''
    return(
      self.gamma - (t-self.phi)%1.0 if (t-self.phi)%1<=self.beta else nan
    )

  def R(self,t):
    '''
    Get the 2D position of the foot at time t. If the foot is not on the 
    ground this returns nan.
    
    @param t: time in cycle (fraction of cycle)
    @type t: float    
    @return: (x,y) position of foot
    @rtype: 2-tuple of floats
    '''
    x = self.X(t)
    return((nan,nan) if isnan(x) else (x,self.delta))  
  
  def __repr__(self):
    '''
    Generate a string representation of the foot that can be evaluated to 
    recreate the foot.
    
    @return: Rep of foot
    @rtype: string
    '''
    return("Foot("+
      "beta="+str(self.beta)+
      ",gamma="+str(self.gamma)+
      ",delta="+str(self.delta)+
      ",phi="+str(self.phi)+
    ")")
  
class Quadruped(object):
  '''
  Wrapper for the four feet which also allows us to calculate the 
  longitudinal stability margin through the limb cycle.
  
  @ivar feet: the feet of the quadruped
  @type beta: list of feet
  '''
  def __init__(self,k):
    '''
    A four footed system, specified via a kinematic gait formula k as 
    detailed in McGhee and Frank.
    
    @param k: kinematic gait formula for the quadruped
    @type k: list of float
    '''
    self.feet = [Foot(k[i],k[4+i],k[8+i],k[12+i]) for i in xrange(4)]
  
  def R(self,t):
    '''
    Feet positions at time t.
    
    @param t: time in cycle (fraction of cycle)
    @type t: float    
    @return: list of (x,y) position of feet
    @rtype: list of 2-tuples of floats
    '''
    return([f.R(t) for f in self.feet])
  
  def eta(self,t):
    '''
    Get the longitudinal stability margins of the different connections 
    between the foot contacts.
    
    @param t: time in cycle (fraction of cycle)
    @type t: float    
    @return: list of (x,y) position of feet
    @rtype: list of 2-tuples of floats
    '''
    x = self.R(t)
    return((
      x[0][0]+x[1][0],
      x[0][0]+x[3][0],
      x[2][0]+x[1][0],
      x[2][0]+x[3][0]
    ))

if __name__=="__main__":
  import doctest
  doctest.testmod()    
  
  try:
    # Make example figure

    from numpy import array,linspace,nanmax,nanmin,minimum
    from pylab import plot,hlines,figure,savefig,xlabel,ylabel
    
    N = 100
    beta = 11.0/12.0 # Duty cycle
    delta = 1.0/4.0 # Lateral position of feet
    k = ( 
      beta,beta,beta,beta, 
      1,1,0,0, 
      delta,-delta,delta,-delta, 
      0.0,1.0/2.0,3.0/4.0,1.0/4.0 
    )
    q = Quadruped(k)
    fp = array([q.R(t) for t in linspace(0,1,N)])
    
    # Plot the foot positions through time, should form solid bars either 
    # side of the body
    figure()
    for m in fp.transpose(1,2,0):
      plot(*m,alpha=0.5,lw=5)
    xlabel("fore-aft foot position (arb)")
    ylabel("lateral foot positoin (arb)")
    savefig("footpositionswalker1.png")

    # Plot the min and max static stabilities for this walker
    # It is stable as the lines do not cross zero.
    figure()
    eta = array([q.eta(t) for t in linspace(0,1,N)])
    plot(nanmax(eta,1))
    plot(nanmin(eta,1))
    hlines(0,0,len(eta))
    xlabel("time (frac. cyc.)")
    ylabel("min. and max. stability margin (arb)")
    savefig("staticstabilitieswalker1.png")

    # According to MG-F then this walker should just be stable
    beta2 = 3.0/4.0
    delta2 = 1.0/4.0
    a = beta2/2 + 1.0
    dg = [0.0,0.0,1.0,1.0]
    k2 = (
      beta2,beta2,beta2,beta2, # Duty cycle
      a+beta2/2.0,a+beta2/2,-a+beta2/2,-a+beta2/2, # Initial foot positions
      delta2,-delta2,delta2,-delta2,
      0.0,1.0/2.0,beta2,beta2-1.0/2.0
    )
    q2 = Quadruped(k2)
    eta2 = array([q2.eta(t) for t in linspace(0,1,N)])
    marg = minimum(nanmax(eta2,1),-nanmin(eta2,1))

    # Plot the longitudinal static stabilities of this walker, it is 
    # just barely stable as the lines touch zero.
    figure()
    plot(nanmax(eta2,1))
    plot(nanmin(eta2,1))
    xlabel("time (frac. cyc.)")
    ylabel("stability margin (arb)")
    hlines(0,0,len(eta2))
    plot(marg)
    savefig("staticstabilitieswalker2.png")
    
  except ImportError as e:
    print "Skipping example figures, this depends on numpy and pylab."