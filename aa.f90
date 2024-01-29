program aa
implicit none
integer i,imax
real fhb,fmax,d,del,dt,ksp,kgs,lam,lama,gam,x,xa,y,ya,p0,ax,axa,ay,ydot,aydot,m,s,aya,t
open(1,file='x') 
fmax=50.;d=60.9
del=1.;dt=1.e-2
imax=5000;m=1.
ksp=0.65;kgs=0.75
lam=0.1;lama=0.1;gam=0.1
x=-100.;xa=-120.;y=1.;ya=1.;s=0.65

do i=1,imax
p0=1./(1.-exp(-(x-xa)))
fhb=-gam*y-kgs*(x-xa-d*p0)
aydot=-gam*y+fhb/m
ay=ydot
ax=y
aya=(kgs/lama)*(x-xa-d*p0)-(fmax/lama)*(1.-s*p0)
axa=ya
ydot=ydot+aydot*dt
y=y+ay*dt
x=x+ax*dt
xa=xa+axa*dt
ya=ya+aya*dt
t=i*dt
write(1,*) x,y
enddo
end
