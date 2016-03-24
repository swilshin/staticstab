'''
Make the plots for the dogs on rough terrain paper. 
'''
from staticstab import Quadruped

from numpy import (linspace,array,minimum,isnan,nanmax,nanmin,where,mean,nan,
  meshgrid,zeros_like,pi)
from pylab import figure,savefig,xlabel,ylabel
from matplotlib.colors import LinearSegmentedColormap

from os.path import join as osjoin

def mfmargins(beta,lda,delta,tpe="optimal",Nsubcyc=100,nanval=nan,mode='min',dg = (0.0,0.0,1.0,1.0)):
  '''
  Calculate a representative stability margin for the kinematic gait formula 
  specified by beta, lda, and delta as detailed in McGhee and Frank (1968)
  
  @param dg: # Direction of shift that lambda performs
  @type dg: 4-tuple of floats
  '''
  a = beta/2.0 + 1.0
  if tpe=="optimal":
    '''
    tpe controls if phase is adjusted with duty cycle for more optimal 
    results. This has no major impact on our core conclusions.
    '''
    k = (
      beta,beta,beta,beta, # Duty cycle
      a+beta/2.0,a+beta/2,-a+beta/2,-a+beta/2, # Initial foot positions
      delta,-delta,delta,-delta,
      0.0+dg[0]*lda,1.0/2.0+dg[1]*lda,beta+dg[2]*lda,beta-1.0/2.0+dg[3]*lda
    )
  elif tpe=="bad":
    k = (
      beta,beta,beta,beta, # Duty cycle
      a+beta/2.0,a+beta/2.0,-a+beta/2.0,-a+beta/2.0, # Initial foot positions
      delta,-delta,delta,-delta,
      0.0+dg[0]*lda,1.0/2.0+dg[1]*lda,0.75+dg[2]*lda,0.25+dg[3]*lda
    )
  else:
    raise "Invalid desired kinematic gait formula."
  q = Quadruped(k)
  eta = array([q.eta(t) for t in linspace(0,1,Nsubcyc)])
  
  marg = minimum(nanmax(eta,1),-nanmin(eta,1))
  
  idxOne = where(4-sum(isnan(eta),1)==1)[0]
  marg[idxOne] = -abs(marg[idxOne])
  if mode=="min":
    smarg = [nanval if isnan(m) else m for m in marg]    
    return(nanval if nan in smarg else min(smarg))
  elif mode=="average":
    mu = mean(marg)
    return(nanval if isnan(mu) else mu)
  else:
    raise "Invalid summary statistic for margin."

'''
Calculate the impact on longitudinal stability of changes in duty cycles 
and shifts in gait towards and away from trot.
'''
Ngrid = 100 # Number of different duty cycles and lambdas to consider
Nsubcyc = 1000 # How many subcycles to calculate to look for minima
delta2 = 1.0/4.0 # Lateral foot position (result is insensitive to this)
b0,b1 = 0.6, 0.9 # Lowest and highest duty cycles to consider
shifts = linspace(-0.5,0.5,Ngrid) # Shifts in lambda
betas = linspace(b0,b1,Ngrid) # Shifts in duty cycles
A,B = meshgrid(betas,shifts) # Make a grid

# Calculate the changes in the stability margins
Mbm = nan*zeros_like(A)
Mba = nan*zeros_like(A)
# These two version of the stability margin show the same behaviour, and 
# take into account the shift in limb phase necessary to restore optimal 
# stability with changing duty cycle.
#Mdm = nan*zeros_like(A)
#Mda = nan*zeros_like(A)


for i in xrange(Ngrid):
  for j in xrange(Ngrid):
    Mbm[i,j] = mfmargins(A[i,j],B[i,j],delta2,'bad',nanval=nan,Nsubcyc=Nsubcyc)
    Mba[i,j] = mfmargins(A[i,j],B[i,j],delta2,'bad',nanval=nan,mode='average',Nsubcyc=Nsubcyc)
    #Mdm[i,j] = mfmargins(A[i,j],B[i,j],delta2,nanval=nan,Nsubcyc=Nsubcyc)
    #Mda[i,j] = mfmargins(A[i,j],B[i,j],delta2,nanval=nan,mode='average',Nsubcyc=Nsubcyc)
    
    
# Plot  the stability margins
# Set up a green colour map for making stability plots
cmgreend = {'red':   ((0.0,  0.5, 0.5),
                   (1.0,  0.5, 0.5)),

         'green': ((0.0,  0.5, 0.5),
                   (1.0,  1.0, 1.0)),

         'blue':  ((0.0,  0.5, 0.5),
                   (1.0,  0.5, 0.5))}
cmgreen = LinearSegmentedColormap('GreenGrey', cmgreend)

# Contour levels for the two plots
V0 = linspace(-0.7,2.5,17)+0.1
V1 = linspace(-0.875,2.5,10)
fs = 18 # Font size

# Contour plot of the minimum stability margin
Bpp = -pi/2-B*2*pi
fig = figure(figsize=(12,6),dpi=100)
ax0 = fig.add_subplot(1,2,1)
ax0.contourf(Bpp,A,Mbm,V0,cmap=cmgreen)
C = ax0.contour(Bpp,A,Mbm,V0,colors='k',linewidths=3)
ax0.clabel(C,c='k',fontsize=fs)
ax0.set_axis_bgcolor('k')
ax0.set_xticks([-3*pi/2,-pi,-pi/2,0,pi/2])
ax0.set_xticklabels(["$-3\pi/2$","$-\pi$","$-\pi/2$","$0$","$\pi/2$"])
ylabel("duty factor (\%)")

# Plot the average stability margin
ax1 = fig.add_subplot(1,2,2)
ax1.contourf(Bpp,A,Mba,V1,cmap=cmgreen)
C = ax1.contour(Bpp,A,Mba,V1,colors='k',linewidths=3)
ax1.clabel(C,c='k',fontsize=fs)
ax1.set_axis_bgcolor('k')
ax1.set_xticks([-3*pi/2,-pi,-pi/2,0,pi/2])
ax1.set_xticklabels(["$-3\pi/2$","$-\pi$","$-\pi/2$","$0$","$\pi/2$"])
xlabel("projected distance to trot $\lambda$ (rad)")

savefig(osjoin("figure","stabilitycontours.pdf"))
savefig(osjoin("figure","stabilitycontours.png"))
savefig(osjoin("figure","stabilitycontours.svg"))
