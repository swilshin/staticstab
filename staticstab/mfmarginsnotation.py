
class Foot(object):
  def __init__(self,beta,gamma,delta,phi):
    self.beta = beta
    self.gamma = gamma
    self.delta = delta
    self.phi = phi
  
  def X(self,t):
    return(
      self.gamma - (t-self.phi)%1.0
        if (t-self.phi)%1<=self.beta else #self.phi%1<=t%1 and t%1<=(self.phi+self.beta)%1 else
      nan
    )

  def R(self,t):
    x = self.X(t)
    return((nan,nan) if isnan(x) else (x,self.delta))  
  
  def __repr__(self):
    return("Foot("+
      "beta="+str(self.beta)+
      ",gamma="+str(self.gamma)+
      ",delta="+str(self.delta)+
      ",phi="+str(self.phi)+
    ")")
  
class Quadruped(object):
  def __init__(self,k):
    self.feet = [Foot(k[i],k[4+i],k[8+i],k[12+i]) for i in xrange(4)]
  
  def R(self,t):
    return([f.R(t) for f in self.feet])
  
  def eta(self,t):
    x = self.R(t)
    return((
      x[0][0]+x[1][0],
      x[0][0]+x[3][0],
      x[2][0]+x[1][0],
      x[2][0]+x[3][0]
    ))

# According to MG-F this walker is stable
beta = 11.0/12.0
delta = 1.0/4.0
k = (
  beta,beta,beta,beta, # Duty cycle
  1,1,0,0, # Initial foot positions
  delta,-delta,delta,-delta,
  0.0,1.0/2.0,3.0/4.0,1.0/4.0
)

f = Foot(0.75,0.3,0.3,0.0)

q = Quadruped(k)
fp = array([q.R(t) for t in linspace(0,1,1000)])
for m in fp.transpose(1,2,0):
  plot(*m,alpha=0.5,lw=5)

figure()
eta = array([q.eta(t) for t in linspace(0,1,1000)])
plot(nanmax(eta,1))
plot(nanmin(eta,1))
hlines(0,0,len(eta))

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
eta2 = array([q2.eta(t) for t in linspace(0,1,1000)])
marg = minimum(nanmax(eta2,1),-nanmin(eta2,1))
#marg = [-10 if isnan(m) else m for m in marg]

figure()
plot(nanmax(eta2,1))
plot(nanmin(eta2,1))
hlines(0,0,len(eta2))
plot(marg)



M = list()
shifts = linspace(-1.0,1.0,1000)
for lda in shifts:
  k3 = (
    beta2,beta2,beta2,beta2, # Duty cycle
    a+beta2/2.0,a+beta2/2,-a+beta2/2,-a+beta2/2, # Initial foot positions
    delta2,-delta2,delta2,-delta2,
    0.0+dg[0]*lda,1.0/2.0+dg[1]*lda,beta2+dg[2]*lda,beta2-1.0/2.0+dg[3]*lda
  )
  q3 = Quadruped(k3)
  eta3 = array([q3.eta(t) for t in linspace(0,1,1000)])
  marg = minimum(nanmax(eta3,1),-nanmin(eta3,1))
  idxOne = where(4-sum(isnan(eta3),1)==1)[0]
  marg[idxOne] = -abs(marg[idxOne])
  smarg = [-100 if isnan(m) else m for m in marg]
  M.append((min(smarg) if not isnan(mean(marg)) else nan,mean(marg)))
figure()
plot(shifts,M)
xlim(-0.5,0.5)
hlines(0,shifts.min(),shifts.max())

# No entry is below -2 so use minus 2 as a minimum to get the plots to work.
def margins(beta2,lda):
  a = beta2/2 + 1.0  
  k3 = (
    beta2,beta2,beta2,beta2, # Duty cycle
    a+beta2/2.0,a+beta2/2,-a+beta2/2,-a+beta2/2, # Initial foot positions
    delta2,-delta2,delta2,-delta2,
    0.0+dg[0]*lda,1.0/2.0+dg[1]*lda,beta2+dg[2]*lda,beta2-1.0/2.0+dg[3]*lda
  )
  q3 = Quadruped(k3)
  eta3 = array([q3.eta(t) for t in linspace(0,1,1000)])
  marg = minimum(nanmax(eta3,1),-nanmin(eta3,1))
  idxOne = where(4-sum(isnan(eta3),1)==1)[0]
  marg[idxOne] = -abs(marg[idxOne])
  smarg = [-2 if isnan(m) else m for m in marg]
  return(min(smarg))


def marginsdfbad(beta2,lda):
  a = beta2/2 + 1.0  
  k3 = (
    beta2,beta2,beta2,beta2, # Duty cycle
    a+beta2/2.0,a+beta2/2,-a+beta2/2,-a+beta2/2, # Initial foot positions
    delta2,-delta2,delta2,-delta2,
    0.0+dg[0]*lda,1.0/2.0+dg[1]*lda,0.75+dg[2]*lda,0.25+dg[3]*lda
  )
  q3 = Quadruped(k3)
  eta3 = array([q3.eta(t) for t in linspace(0,1,1000)])
  marg = minimum(nanmax(eta3,1),-nanmin(eta3,1))
  idxOne = where(4-sum(isnan(eta3),1)==1)[0]
  marg[idxOne] = -abs(marg[idxOne])
  smarg = [-2 if isnan(m) else m for m in marg]
  return(min(smarg))

def marginsdfbadaverage(beta2,lda):
  a = beta2/2 + 1.0  
  k3 = (
    beta2,beta2,beta2,beta2, # Duty cycle
    a+beta2/2.0,a+beta2/2,-a+beta2/2,-a+beta2/2, # Initial foot positions
    delta2,-delta2,delta2,-delta2,
    0.0+dg[0]*lda,1.0/2.0+dg[1]*lda,0.75+dg[2]*lda,0.25+dg[3]*lda
  )
  q3 = Quadruped(k3)
  eta3 = array([q3.eta(t) for t in linspace(0,1,1000)])
  marg = minimum(nanmax(eta3,1),-nanmin(eta3,1))
  idxOne = where(4-sum(isnan(eta3),1)==1)[0]
  marg[idxOne] = -abs(marg[idxOne])
  #smarg = [-2 if isnan(m) else m for m in marg]
  mu = mean(marg)
  return(mu if not isnan(mu) else nan)#-2)

def marginsaverage(beta2,lda):
  a = beta2/2 + 1.0
  k3 = (
    beta2,beta2,beta2,beta2, # Duty cycle
    a+beta2/2.0,a+beta2/2,-a+beta2/2,-a+beta2/2, # Initial foot positions
    delta2,-delta2,delta2,-delta2,
    0.0+dg[0]*lda,1.0/2.0+dg[1]*lda,beta2+dg[2]*lda,beta2-1.0/2.0+dg[3]*lda
  )
  q3 = Quadruped(k3)
  eta3 = array([q3.eta(t) for t in linspace(0,1,1000)])
  marg = minimum(nanmax(eta3,1),-nanmin(eta3,1))
  idxOne = where(4-sum(isnan(eta3),1)==1)[0]
  marg[idxOne] = -abs(marg[idxOne])
  #smarg = [-2 if isnan(m) else m for m in marg]
  mu = mean(marg)
  return(mu if not isnan(mu) else -2)

N = 200
b0 = 0.6
b1 = 0.9
shifts = linspace(-0.5,0.5,N)
betas = linspace(b0,b1,N)
A,B = meshgrid(betas,shifts)

Mbm = nan*zeros_like(A)
for i in xrange(N):
  for j in xrange(N):
    Mbm[i,j] = marginsdfbad(A[i,j],B[i,j])

Mba = nan*zeros_like(A)
for i in xrange(N):
  for j in xrange(N):
    Mba[i,j] = marginsdfbadaverage(A[i,j],B[i,j])

Mda = nan*zeros_like(A)
for i in xrange(N):
  for j in xrange(N):
    Mda[i,j] = marginsaverage(A[i,j],B[i,j])

Mdm = nan*zeros_like(A)
for i in xrange(N):
  for j in xrange(N):
    Mdm[i,j] = margins(A[i,j],B[i,j])

    
# Interpolate and plot

def interpM(rat0,rat1,N,M):
  shiftsp = linspace(-0.5,0.5,N*rat0)
  betasp = linspace(b0,b1,N*rat1)
  Ap,Bp = meshgrid(betasp,shiftsp)
  Mdmp = zeros_like(Ap)

  s = RectBivariateSpline(shifts,betas,M,kx=1,ky=1,s=0)
  Mp = s(shiftsp,betasp)
  return(Mp,Ap,Bp)
  
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import RectBivariateSpline


rat0 = 2
rat1 = 2

#Mdmp,Ap,Bp = interpM(rat0,rat1,N,Mdm)
#Mdap,Ap,Bp = interpM(rat0,rat1,N,Mda)
#Mbmp,Ap,Bp = interpM(rat0,rat1,N,Mbm)
#Mbap,Ap,Bp = interpM(rat0,rat1,N,Mba)

from matplotlib.colors import LinearSegmentedColormap
cmgreend = {'red':   ((0.0,  0.5, 0.5),
                   (1.0,  0.5, 0.5)),

         'green': ((0.0,  0.5, 0.5),
                   (1.0,  1.0, 1.0)),

         'blue':  ((0.0,  0.5, 0.5),
                   (1.0,  0.5, 0.5))}
cmgreen = LinearSegmentedColormap('GreenGrey', cmgreend)
               
V0 = linspace(-0.7,2.5,17)+0.1
V1 = linspace(-0.875,2.5,10)
fs = 18
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

ax1 = fig.add_subplot(1,2,2)
ax1.contourf(Bpp,A,Mba,V1,cmap=cmgreen)
C = ax1.contour(Bpp,A,Mba,V1,colors='k',linewidths=3)
ax1.clabel(C,c='k',fontsize=fs)
ax1.set_axis_bgcolor('k')
ax1.set_xticks([-3*pi/2,-pi,-pi/2,0,pi/2])
ax1.set_xticklabels(["$-3\pi/2$","$-\pi$","$-\pi/2$","$0$","$\pi/2$"])
'''
ax2 = fig.add_subplot(2,2,3)
ax2.contourf(Bpp,Ap,Mbmp,V0,cmap=cmgreen)
C = ax2.contour(Bpp,Ap,Mbmp,V0,colors='k',linewidths=3)
ax2.clabel(C,c='k',fontsize=fs)
ax2.set_axis_bgcolor('k')
ax2.set_xticks([-3*pi/2,-pi,-pi/2,0,pi/2])
ax2.set_xticklabels(["$-3\pi/2$","$-\pi$","$-\pi/2$","$0$","$\pi/2$"])


ax3 = fig.add_subplot(2,2,4)
ax3.contourf(Bpp,Ap,Mbap,V1,cmap=cmgreen)
C = ax3.contour(Bpp,Ap,Mbap,V1,colors='k',linewidths=3)
ax3.clabel(C,c='k',fontsize=fs)
ax3.set_axis_bgcolor('k')
ax3.set_xticks([-3*pi/2,-pi,-pi/2,0,pi/2])
ax3.set_xticklabels(["$-3\pi/2$","$-\pi$","$-\pi/2$","$0$","$\pi/2$"])
'''
xlabel("projected distance to trot $\lambda$ (rad)")


savefig("stabilitycontours.pdf")
savefig("stabilitycontours.png")
savefig("stabilitycontours.svg")

idealMarg = Mbm[:,0]
val = B[:,0]

#
#
#
#
#
#
#
from csv import reader
from scipy.stats import gaussian_kde
from os.path import join as osjoin
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import PathPatch

# Set fonts
rc('text', usetex=True)
rc('font', **{'family':'sans-serif','sans-serif':["Computer Modern"],'size':24})

# ====================
# Load Data
# ====================
with open("strPROJroughGaitsWALKMichelleFILTERED.csv",'rb') as rgf:
  rg = reader(rgf)
  header = rg.next()
  dat = array([row for row in rg])
dat = dict(zip(header,dat.T))
dat["ProjTrot"] = array(dat["ProjTrot"],dtype=float)
for i in xrange(1,4):
  dat["Phi"+str(i)] = array(dat["Phi"+str(i)],dtype=float)
phi = array([dat["Phi"+str(i)] for i in xrange(1,4)])

dogs = set(dat["Dog"])
terrains = set(dat["Terrain.Type"])

markers = dict(zip(dogs,['x','+','*','o','|','_']))
cols = {"Flat":'b',"Rough":'r'}

pad = 0.1
theta = linspace(min(dat["ProjTrot"])-pad, max(dat["ProjTrot"])+pad,1000)
dens = dict()
rast = dict()
for dog in dogs:
  for terrain in terrains:
    rast[dog,terrain] = dat["ProjTrot"][
      logical_and(dat["Dog"]==dog,
      dat["Terrain.Type"]==terrain)
    ]
    dens[dog,terrain] = gaussian_kde(rast[dog,terrain])(theta)
dogOrder = zip(*sorted([(dog,mean(rast[dog,'Flat'])) for dog in dogs],key=lambda x: x[1]))
lr = dict(fliplr(array(zip(arange(len(dogs))%2,*dogOrder))[:,:2]))
 
for dog in dogs:
  lr[dog]=int(lr[dog])

# ====================
# KDE Plots
# ====================
fig = figure(figsize=(12,11),dpi=100)
peakDens = max([dens[dog,terrain].max() for dog in dogs for terrain in terrains])
rastShift = 0.15
ax=None
for i,dog in enumerate(dogs):
  ax = subplot(1,len(dogs),1+i,sharex=ax,sharey=ax)
  ax.boxplot([rast[dog,terrain] for terrain in terrains],positions=[-3.4,3.4])
  for j,terrain in enumerate(terrains):
    plot((-1)**(j+1)*dens[dog,terrain],theta,color=cols[terrain],label=dog,lw=3)
    pc = ax.fill_betweenx(theta,(-1)**(j+1)*dens[dog,terrain],color=cols[terrain],alpha=0.2)
  if i!=0:
    ax.get_yaxis().set_visible(False)
  else:
    yticks([-3*pi/4,-pi/2,-pi/4,0],["$-3\pi/4$","$-\pi/2$","$-\pi/4$","$0$"])
    ylabel("projected distance to trot $\lambda$ (rad)")
  #ylim(min(theta),max(theta))
  ylim(min(theta),0.5)
  if i<len(dogs)-1:
    ax.spines['right'].set_visible(False)
    if i>0:
      ax.spines['left'].set_visible(False)
    else:
      ax.spines['left'].set_visible(True)
    ax.tick_params(right="off")
  else:
    ax.spines['right'].set_visible(True)
    ax.tick_params(left="off",right="on",direction='in',which='major')
    ax.spines['left'].set_visible(False)
  ax.set_xlim([-4,4])
  ax.set_xticks([-4,-2,0,2,4])
ax.set_xticklabels([int(sign(itick)*itick) if itick%2==0 else "" for itick in ax.get_xticks()])
fig.text(.5, .03, 'frequency ($rad^{-1}$)', horizontalalignment='center')
savefig(osjoin("fig","kerneldensityplotsideon.pdf"))
savefig(osjoin("fig","kerneldensityplotsideon.png"))
savefig(osjoin("fig","kerneldensityplotsideon.svg"))


# ====================
# Make 3D Scatter plot
# ====================
fig = figure(figsize=(12,12),dpi=100)
ax = fig.add_subplot(111, projection='3d')
for terrain in terrains:
  for dog in dogs:
    theta = phi[:,where(logical_and(dat["Dog"]==dog,dat["Terrain.Type"]==terrain))[0]]
    ax.scatter(((theta[0]+pi)%(2*pi))-pi,theta[1],((theta[2]+pi)%(2*pi))-pi,marker=markers[dog],s=100,color=cols[terrain],alpha=0.5)

def plotCross(ax,c,cs=pi/10.0,col='k',txt=None,txtShift=pi/5.0):
  ax.plot([c[0]-cs,c[0]+cs],[c[1],c[1]],[c[2],c[2]],lw=2,c=col)
  ax.plot([c[0],c[0]],[c[1]-cs,c[1]+cs],[c[2],c[2]],lw=2,c=col)
  ax.plot([c[0],c[0]],[c[1],c[1]],[c[2]-cs,c[2]+cs],lw=2,c=col)
  if txt!=None:
    ax.text(c[0], c[1]+txtShift, c[2], txt)
walkPh = [pi/2,pi,-pi/2]
trotPh = [0,pi,0]
#pronkPhs = [[0,0,0],[0,2*pi,0]]
pacePhs = [[pi,pi,-pi],[-pi,pi,pi]]#[[pi,pi,-pi],[pi,pi,pi],[-pi,pi,pi],[-pi,pi,-pi]]
plotCross(ax,trotPh,txt='trot')
plotCross(ax,walkPh,pi/5.0,col='#505050',txt='walk')
#for pronkPh in pronkPhs:
#  plotCross(ax,pronkPh,txt='pronk')
for pacePh in pacePhs:
  plotCross(ax,pacePh,txt='pace')

#ax.scatter(pi/2,pi,3*pi/2,color='k',marker='8',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(0,pi,0,color='k',marker='D',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(2*pi,pi,0,color='k',marker='D',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(0,pi,2*pi,color='k',marker='D',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(2*pi,pi,2*pi,color='k',marker='D',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(2*pi,2*pi,2*pi,color='k',marker='p',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(0,2*pi,2*pi,color='k',marker='p',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(0,0,2*pi,color='k',marker='p',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(0,0,0,color='k',marker='p',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(2*pi,0,2*pi,color='k',marker='p',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(2*pi,2*pi,0,color='k',marker='p',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(2*pi,0,0,color='k',marker='p',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(0,2*pi,0,color='k',marker='p',s=1600,facecolor='none',lw=5,alpha=0.9)
#ax.scatter(pi,pi,pi,color='k',marker='h',s=1600,facecolor='none',lw=5,alpha=0.9)
ax.plot([-pi,pi],[pi,pi],[pi,-pi],c='k',lw=2)
ax.set_xlim(-pi,pi)
ax.set_ylim(0,2*pi)
ax.set_zlim(-pi,pi)
ax.set_xticks([-pi,-pi/2,0,pi/2,pi])
ax.set_xticklabels(["$-\pi$","$-\pi/2$","$0$","$\pi/2$","$\pi$"])
ax.set_yticks([0,pi/2,pi,3*pi/2,2*pi])
ax.set_yticklabels(["$0$","$\pi/2$","$\pi$","$3\pi/2$","$2\pi$"])
ax.set_zticks([-pi,-pi/2,0,pi/2,pi])
ax.set_zticklabels(["$-\pi$","$-\pi/2$","$0$","$\pi/2$","$\pi$"])
ax.set_xlabel("$\phi_1$")
ax.set_ylabel("$\phi_2$")
ax.set_zlabel("$\phi_3$")
tight_layout()
savefig(osjoin("fig","gaitusage3D.pdf"))
savefig(osjoin("fig","gaitusage3D.png"))
savefig(osjoin("fig","gaitusage3D.svg"))



# Load Rough gaits data
phaseDat = dict([(t,hstack([rast[(k,t)] for k in set(zip(*rast.keys())[0])])) for t in set(zip(*rast.keys())[1])])
fig = figure(figsize=(12,11),dpi=100)
for i in xrange(2):
  ax = fig.add_subplot(1,2,i+1)
  b = ax.boxplot([phaseDat[t] for t in phaseDat],False,vert=False,positions=[0,0.3])
  #ax.set_xlabel("projected distance to trot $\lambda$ (rad)")
  ax.set_xlim(-3*pi/2,pi/2)
  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_ylim((-0.2,0.5))
  cols=['r','b']
  for k in b:
    for i,ele in enumerate(b[k]):
      ele.set_color(cols[i<len(b[k])/2])
  ax.set_xticks([-3*pi/2,-pi,-pi/2,0,pi/2])
  ax.set_xticklabels(["$-3\pi/2$","$-\pi$","$-\pi/2$","$0$","$\pi/2$"])
savefig("usedgaitscollapsed.pdf")
savefig("usedgaitscollapsed.png")
savefig("usedgaitscollapsed.svg")

'''
fig = figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(Bp, Ap, Mdmp[:,:], rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
ax.set_xlabel("phase shift from ideal (fraction cycle)")
ax.set_ylabel("duty cycle (fraction cycle)")
ax.set_zlabel("minimum stability margin (stride lengths)")
ax.set_zlim3d([-2.0,2.2])
savefig("stabmarginvarydutyoptimalmin.png")

fig = figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(Bp, Ap, Mdap[:,:], rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
ax.set_xlabel("phase shift from ideal (fraction cycle)")
ax.set_ylabel("duty cycle (fraction cycle)")
ax.set_zlabel("average stability margin (stride lengths)")
ax.set_zlim3d([-2.0,2.2])
savefig("stabmarginvarydutyoptimalaverage.png")

fig = figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(Bp, Ap, Mbmp[:,:], rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
ax.set_xlabel("phase shift from ideal (fraction cycle)")
ax.set_ylabel("duty cycle (fraction cycle)")
ax.set_zlabel("minimum stability margin (stride lengths)")
ax.set_zlim3d([-2.0,2.2])
savefig("stabmarginbaddutyoptimalmin.png")

fig = figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(Bp, Ap, Mbap[:,:], rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
ax.set_xlabel("phase shift from ideal (fraction cycle)")
ax.set_ylabel("duty cycle (fraction cycle)")
ax.set_zlabel("average stability margin (stride lengths)")
ax.set_zlim3d([-2.0,2.2])
savefig("stabmarginbaddutyoptimalaverage.png")

dx,dy,dz = ax.get_xlim3d(),ax.get_ylim3d(),ax.get_zlim3d()
fig = figure()
ax = fig.gca(projection='3d')
ax.plot([0],[0],[0])
ax.set_xlabel("phase shift from ideal (fraction cycle)")
ax.set_ylabel("duty cycle (fraction cycle)")
ax.set_zlabel("average stability margin (stride lengths)")
ax.set_xlim3d(dx)
ax.set_ylim3d(dy)
ax.set_zlim3d(dz)
savefig("stabmarginblank.png")
'''