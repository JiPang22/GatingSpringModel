program aa
implicit none
integer i,imax,k,n
parameter(imax=5000)
real fhb,fmax,d,del,dt,ksp,kgs,lam,lama,gam,x,xa,y,p0,ax,axa,ay,m,s,t,sumi,sumr,om,amp
real xt(imax)
parameter(gam=0.3,dt=1.e-2,lam=2.8e-3,lama=10.e-3,d=60.9,s=0.65,m=1.,ksp=0.65,kgs=0.75,del=4.44)

open(1,file='aa') 
open(2,file='bb') 

!do n=1,300
x=1.;xa=1.;y=1.
!fmax=(n-1.)*0.2
fmax=200

do i=1,imax

p0=1./(1.+exp(10.)*exp(-(x-xa)))!probability
fhb=-gam*y-kgs*(x-xa-d*p0)-ksp*x!applied force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ay=-gam*y+fhb/m
ax=y

axa=(kgs/lama)*(x-xa-d*p0)-(fmax/lama)*(1.-s*p0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y=y+ay*dt!velocity of x
x=x+ax*dt!position of hair bundle
xa=xa+axa*dt!position of motor protein
write(1,*) i*dt,x
!write(1,*) x,xa
xt(i)=x
enddo

!DFT
do k=1,40
om=6.28*k/(imax*dt)
!reset
t=0.
sumr=0.
sumi=0.
!sum
do i=1,imax
t=i*dt
sumr=sumr+xt(i)*cos(om*t)*dt
sumi=sumi-xt(i)*sin(om*t)*dt
enddo
amp=sqrt(sumi**2+sumr**2)
write(2,*) om,amp
enddo
write(2,*) ''
!enddo
end
