program shmmidpt
implicit none

integer(4):: i, n
real(8), dimension(:), allocatable:: x,v,a,k,u,e,t
real(8):: m, c, tau, pi, h, w

open(unit=100, file= 'shmxvt.dat')
open(unit=200, file= 'shmvvt.dat')
open(unit=300, file= 'shmavt.dat')
open(unit=400, file= 'shmuvt.dat')
open(unit=500, file= 'shmkvt.dat')
open(unit=600, file= 'shmevt.dat')

m=1.0
 c=1.0
a=1.0
w=sqrt(c/m)
h=0.1
pi=4*atan(1.0)
tau=2*pi*sqrt(m/c)
!nh=number of timestops, number of values -1 (initial)
n=floor(10*tau/h)

allocate (x(n), v(n), a(n), k(n), u(n), e(n), t(n))

t(1)=0.0
x(1)=1.0
v(1)=0.0
u(1)=(1/2)*c*x(1)**2
k(1)=(1/2)*m*v(1)**2
e(1)=k(1)+u(1)

do i= 1,(n-1)

!Program x, ... arrays
!(i) refers to array values programmed from do i= ...

a(i+1)=-(w**2)*x(i)
v(i+1)=v(i)+h*a(i)
x(i+1)=x(i)+h*((v(i)+v(i+1))/2)
t(i+1)=t(i)+h

!write (unit,*) x(i+1) ...

write(100,*) t(i), x(i)
write(200,*) t(i), v(i)
write(300,*) t(i), a(i)

u(i+1)=(0.5)*c*x(i+1)**2
k(i+1)=(0.5)*m*v(i+1)**2
e(i+1)=k(i+1)+u(i+1)

write(400,*) t(i), u(i)
write(500,*) t(i), k(i)
write(600,*) t(i), e(i)


enddo

endprogram shmmidpt
