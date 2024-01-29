program aa
implicit none
integer i,imax
real fhb,fmax,d,del,dt,ksp,kgs,lam,lama,gam,x,xa,y,p0,ax,axa,ay,m,s
open(1,file='aa') 
imax=5000;dt=1.e-2
lam=0.1;lama=0.1;gam=0.1
fmax=52.;d=60.9;s=0.65;m=1.
ksp=0.65;kgs=0.75;del=1.
x=-145.;xa=-135.;y=1.

do i=1,imax
p0=1./(1.+exp(-(x-xa)))!probability
fhb=-gam*y-kgs*(x-xa-d*p0)-ksp*x!applied force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ay=-gam*y+fhb/m
ax=y

axa=(kgs/lama)*(x-xa-d*p0)-(fmax/lama)*(1.-s*p0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y=y+ay*dt!velocity of x
x=x+ax*dt!position of hair bundle
xa=xa+axa*dt!position of motor protein
!write(1,*) i*dt,x
write(1,*) x,xa
enddo
end
