
from scipy.optimize import fmin
from numpy import nan,sum,min,pi,array,linspace,where,diff,argmax,isnan
from staticstab import shiftedFootfallPattern
from pylab import figure,plot

def mfConfigs(x,eta):
  if sum(x)<2:
    return((nan,))
  
  # Six cases with two feet on the ground
  if x==(1,0,1,0):
    return((nan,))    
  if x==(0,1,1,0):
    return((eta[0],))
  if x==(0,0,1,1):
    return((eta[2],))
  if x==(1,1,0,0):
    return((eta[1],))
  if x==(1,0,0,1):
    return((eta[1]+eta[2]-eta[0],))
  if x==(0,1,0,1):
    return((nan,))
  if x==(1,1,0,0):
    return((eta[0],eta[1]))
  
  # Four cases with three feet on ground
  if x==(1,0,1,1):
    return((eta[2],eta[1]+eta[2]-eta[0]))
  if x==(0,1,1,1):
    return((eta[0],eta[2]))
  if x==(1,1,0,1):
    return((eta[1],eta[1]+eta[2]-eta[0]))
  if x==(1,1,1,0):
    return((eta[0],eta[1]))
  
  # One case with four feet on ground
  if x==(1,1,1,1):
    return((eta[0],eta[1]+eta[2]-eta[0]))

def mfMargin(x,eta,tau,lda,x0=0.0):
  c = x0+tau*lda/(2*pi) # CoM position
  ed = mfConfigs(x,eta) # Grab edges
  if len(ed)==1:
    return(-abs(ed[0]-c))
  else:
    assert ed[0]<ed[1]
    if c<ed[0]:
      return(-abs(ed[0]-c))
    if c>ed[1]:
      return(-abs(ed[1]-c))
    return(min([abs(ed[0]-c),abs(ed[1]-c)]))

def mfStepMarg(tau,eta,lda,g,d,nanPen=None,x0=0.0):
  f,t = shiftedFootfallPattern(d,g)
  i = (where(t>=tau)[0][0]-1)%(len(t))
  if nanPen==None:
    return(mfMargin(f[i],eta,tau,lda,x0))
  elif nanPen=='max':
    m = mfMargin(f[i],eta,tau,lda,x0)
    return(-abs(array(eta)-(x0+tau*lda/(2*pi))).max() if isnan(m) else m)
  else:
    m = mfMargin(f[i],eta,tau,lda,x0)
    return(nanPen if isnan(m) else m)

def etaToLatin(eta):
  return(
    -1,
    2*(eta[1]-eta[0]-0.5),
    2*(eta[0]+0.5),
    2*(eta[2]+0.5)
  )

def xToEta(x):
  return((x[0],x[0]+abs(x[1]),x[0]+abs(x[2])))

def omegaToEta(w):
  return((
    (w[1]+w[2])/2,
    (w[0]+w[1])/2,
    (w[2]+w[3])/2,
    (w[3]+w[0])/2,
  ))

d = 3*pi/2
g = (pi/2,pi,-pi/2)
lda = 1.0
ff = [array(shiftedFootfallPattern(d,g)[1])[(diff(array(shiftedFootfallPattern(d,g)[0]).T)==1)[i]][0] for i in xrange(4)]
w = array([1,-1,-1,1])+d/(2*pi)+array(ff)*lda/(2*pi)
eta = omegaToEta(w)
dg = array([1,0,-1])
x0 = -0.05
taus = linspace(0,2*pi,100)
xs = linspace(-1.0,1.0,100)
nanPen = 'max' #-10.0

def fscore(x,g):
  return(sum([mfStepMarg(tau,eta,lda,g,d,nanPen,x) for tau in taus]))

def gscore(x,g):
  return(min([mfStepMarg(tau,eta,lda,g,d,nanPen,x) for tau in taus]))

sc = [fscore(ix0,g) for ix0 in xs]
x0 = xs[argmax(sc)]
x1 = fmin(lambda x: -fscore(x,g),[x0],xtol=1e-9,ftol=1e-9)

figure()
plot(taus,[-mfStepMarg(tau,eta,lda,g,d,nanPen,x1) for tau in taus],c='r')
plot(taus,[-mfStepMarg(tau,eta,lda,g+0.1*dg,d,nanPen,x1) for tau in taus],c='g')
plot(taus,[-mfStepMarg(tau,eta,lda,g-0.1*dg,d,nanPen,x1) for tau in taus],c='k')

xn = [fmin(lambda x: -fscore(x,g+q*dg),[x1],xtol=1e-9,ftol=1e-9) for q in taus-pi]
figure()
for x in linspace(-1.0,1.0,20):
  plot(taus,cumsum([mfStepMarg(tau,eta,lda,g,d,nanPen,x) for tau in taus]))

figure()
plot(taus-pi,[-fscore(x,g+q*dg) for x,q in zip(xn,taus-pi)],c='r')
plot(taus-pi,[-fscore(x1,g+q*dg) for q in taus-pi],c='b')

# Minimal margin
sc = [gscore(ix0,g) for ix0 in xs]
z0 = xs[argmax(sc)]

z1 = fmin(lambda x: -gscore(x,g),[z0],xtol=1e-9,ftol=1e-9)
zn = [fmin(lambda x: -gscore(x,g+q*dg),[z1],xtol=1e-9,ftol=1e-9) for q in taus-pi]

figure()
plot(taus,[-mfStepMarg(tau,eta,lda,g,d,nanPen,z1) for tau in taus],c='r')
plot(taus,[-mfStepMarg(tau,eta,lda,g+0.1*dg,d,nanPen,z1) for tau in taus],c='g')
plot(taus,[-mfStepMarg(tau,eta,lda,g-0.1*dg,d,nanPen,z1) for tau in taus],c='k')

figure()
plot(taus-pi,[-gscore(x,g+q*dg) for x,q in zip(zn,taus-pi)],c='r')
plot(taus-pi,[-gscore(z1,g+q*dg) for q in taus-pi],c='b')
