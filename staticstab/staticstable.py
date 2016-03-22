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
  
  Polyped - Instantiated with a McGhee and Frank style limb configuration 
  specification. uses L{footfallPattern} to calculate footfall patterns. 
  Can calculate the polygon of support and stability margin.

MAIN FUNCTIONS
==============

  footfallPattern - Calculate the foot on and off events and foot contact 
  patterns from a set of relative phases and duty factors.

TYPICAL USAGE
=============

  >>> from numpy import array,ones,pi
  >>> g = 2.0*pi*ones((3))/4.0
  >>> d = 2.0*pi*ones((4))*3.0/4.0
  >>> x = array([[-1,1],[1,-1],[-1,-1],[1,1]])
  >>> lda = 0.7
  >>> p = Polyped(g,d,x,lda)
  >>> p.footPositions(pi)
  array([[-1.46666667,  1.        ],
         [ 0.76666667, -1.        ],
         [ 0.53333333,  1.        ]])

THEORY OF OPERATION
===================

  This class implements a model very similar to that of McGhee and Frank, 
  although based on quasi-static stability rather than longitudinal static 
  stability. The gait is specified in exactly the same way, as in this paper 
  and uses the Shapely library to find the things like the convex hull and 
  the distance between polygons and points (like the convex hull and the 
  centre of mass.
  
  All phases are stored in radians

REFERENCES
==========

  McGhee, Robert B., and Frank, Andrew A. "On the stability properties of 
  quadruped creeping gaits." Mathematical Biosciences 3 (1968): 331-351.

@author: Simon Wilshin
@contact: swilshin@rvc.ac.uk
@date: Dec 2014
'''

from numpy import (zeros,hstack,cumsum,array,pi,array,where,outer,abs,sum,nan,
  isnan,inf,diff,vstack,hstack)

from scipy.spatial import ConvexHull
from shapely.geometry import MultiPoint,Point,Polygon,LineString

def footfallPattern(d,phi,eps=1e-10,fullOutput=False):
  '''
  Given a duty factor d, and a set of relative phases, phi this computes the 
  resulting footfall patterns. It returns two lists, the first contains the 
  footfall patterns, the second the onset phase for these patterns. Tolerance 
  can be set to remove very short events which result from floating point 
  precision errors. If fullOutput is true it returns four lists, the first 
  two being the foot on and foot off timings and the last pair as before.
  
  TYPICAL USAGE
  =============

    >>> from numpy import array,ones,pi
    >>> g = 2.0*pi*ones((3))/4.0
    >>> d = 2.0*pi*ones((4))*3.0/4.0
    >>> # Calculate the patterns the foot form
    >>> footfallPattern(d,g)[0]
    [(1, 0, 1, 1), (1, 1, 0, 1), (1, 1, 1, 0), (0, 1, 1, 1), (1, 0, 1, 1)]

  @param d: duty cycle of leg (radians)
  @type d: float, array
  @param phi: relative phases of legs (radians)
  @type phi: float
  @param eps: cut off for floating pont precision
  @type eps: float  
  @param fullOutput: return foot on / off events as well as foot patterns / 
    events
  @type fullOutput: bool
  @return: foot patterns and associated events, if fullOutput is true then 
    return foot on / off events and the aforementioned
  @rtype: tuple of lists
  '''
  # Calculate foot on and off phases
  fon = zeros(len(phi)+1)
  fon[1:] = cumsum(phi)%(2*pi)
  foff = (fon+d)%(2*pi)
  
  # Enumerate events by creating a duplicated list of them and culling said 
  # duplicates
  events = array(sorted(hstack([fon,foff])))
  events0 = None
  while events0==None or len(events0)!=len(events) or all(events0!=events):
    events0 = events
    events = hstack([events,events[0]])
    events = events[abs(events[1:]-events[:-1])%(2*pi)>eps]
  events = hstack([events[-1]-2*pi,events,events[0]+2*pi])
  
  # Calculate the foot configurations immediately after each event
  feet = list()
  for event in events:
    cfig = array([
      (fon-event-eps)%(2*pi)-2*pi+eps,
      (foff-event-eps)%(2*pi)-2*pi+eps
    ])
    cfig *= 1-abs(cfig)<eps
    feet.append(tuple([int(i) for i in cfig[0]>=cfig[1]]))
  
  # Trim superfluous events resulting from float point precision errors
  feet = array(feet)[events>-eps]
  events = events[events>-eps]
  feet = [tuple(f) for f in feet]
  
  # Return either all calculations or just feet configurations and events
  if fullOutput:
    return(fon,foff,feet,events)
  return(feet,events)

def shiftedFootfallPattern(d,g,eps=1e-10,fullOutput=False):
  q = footfallPattern(d,g,eps,fullOutput)
  f,p = q[:2]
  Lf = len(p)
  idx = where((diff(f,axis=0)==1)[:,0])[0][0]
  f = vstack([f,f[1:]])
  p = hstack([p,p[1:]+2*pi])
  if len(q)==2:
    return(
      [tuple(x) for x in f[idx:idx+Lf]],
      p[idx:idx+Lf]-p[idx],
    )
  else:
    return(
      [tuple(x) for x in f[idx:idx+Lf]],
      p[idx:idx+Lf]-p[idx],
      q[3],
      q[4]
    )
    

class Polyped(object):
  '''
  Represents the polyped locomotion system.
  
  Compute limb configurations, polygons of support and stability margins of a 
  McGhee and Frank parametrized polyped. Requires the relative phases of the 
  limbs, the duty cycle of the limbs, the the touchdown positions of the feet 
  and the stride lengths.
  
  TYPICAL USAGE
  =============

    >>> from numpy import array,ones,pi
    >>> g = 2.0*pi*ones((3))/4.0
    >>> d = 2.0*pi*ones((4))*3.0/4.0
    >>> x = array([[-1,1],[1,-1],[-1,-1],[1,1]])
    >>> lda = 0.7
    >>> p = Polyped(g,d,x,lda)
    >>> p.footPositions(pi)
    array([[-1.46666667,  1.        ],
           [ 0.76666667, -1.        ],
           [ 0.53333333,  1.        ]])
  
  @ivar g: relative phases
  @type g: float, array
  @ivar d: duty cycle in radians
  @type d: float, array
  @ivar x: foot touchdown position
  @type x: float, array
  @ivar lda: stride length
  @type lda: float
  @ivar fon: foot touch-down phases
  @type fon: float, array
  @ivar foff: foot lift-off phases
  @type foff: float, array
  @ivar feet: foot contact patterns associated with events
  @type feet: float, array
  @ivar events: phases of different event boundaries
  @type events: float, array
  @cvar origin: The location of the centre of mass, origin by default 
  @type origin: float, array
  @cvar e0: direction of movement, should be unit vector
  @type e0: float, array
  '''
  def __init__(self,g,d,x,lda,eps=1e-10):
    '''
    Constructor, see class documentation and module documentation for further 
    details.
    
    @param g: relative phases
    @type g: float, array
    @param d: duty cycle in radians
    @type d: float, array
    @param x: foot touchdown position
    @type x: float, array
    @param lda: stride length
    @type lda: float
    @param eps: cut off for floating pont precision
    @type eps: float    
    '''
    self.g = g
    self.d = d
    self.x = x
    self.lda = lda
    self.fon,self.foff,self.feet,self.events = footfallPattern(d,g,fullOutput=True)
  
  def footPositions(self,t):
    '''
    Get the position of the feet on the ground at phase t of this polyped.
    
    @param t: phase to calculate foot positions at
    @type t: float
    @return: foot positions at phase t
    @rtype: float, array
    '''
    i = (where(self.events>=t)[0][0]-1)%(len(self.events))
    f = array(self.feet[i],dtype=bool)
    return(self.x[f] - outer(
      (( 
        (t-self.fon[f])*(t>=self.fon[f]) +
        (2*pi-self.fon[f]+t)*(t<self.fon[f]) 
      )*self.lda)/self.d[f],Polyped.e0)
    )
  
  def support(self,t):
    '''
    Get the polygon of support at phase t.

    @param t: phase to calculate polygon of support at
    @type t: float
    @return: polygon of support at phase t
    @rtype: Polygon
    '''
    y = self.footPositions(t)
    if y.shape[0]==0:
      return(nan)
    return(MultiPoint(y).convex_hull)
  
  def margin(self,t):
    '''
    Get the stability margin at phase t, the shortest distance between the 
    origin (assumed centre of mass) and the polygon of support. Can change 
    the origin by changing L{Polyped.origin}
    
    @param t: phase to calculate stability margin at
    @type t: float
    @return: stability margin at phase t
    @rtype: Polygon
    '''
    h = self.support(t)
    if (
      not isinstance(h,Polygon)
        and
      not isinstance(h,Point)
        and
      not isinstance(h,LineString)
        and
      isnan(h)
    ):
      return(-inf)
    l = h.boundary
    return(((-1)**(h.contains(Polyped.origin)+1))*l.distance(Polyped.origin))
Polyped.origin = Point([0.0,0.0])
Polyped.e0 = Point([1.0,0.0])    
    

if __name__=="__main__":
  import doctest
  doctest.testmod()