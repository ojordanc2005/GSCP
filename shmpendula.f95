program shmpendula
implicit none

integer(4):: i, n
real(8), dimension(:), allocatable:: x,v,a,k,u,e,t
!rotational analogs
real(8):: m, tau, pi, h, l, g

open(unit=110, file= 'shmxvt.dat')
open(unit=210, file= 'shmvvt.dat')
open(unit=310, file= 'shmavt.dat')
open(unit=410, file= 'shmuvt.dat')
open(unit=510, file= 'shmkvt.dat')
open(unit=610, file= 'shmevt.dat')

m=1.0
!a=1.0
h=0.1
pi=4*atan(1.0)
n=floor((10.0*tau)/h)
l=1.0
g= 9.8
tau=2*pi*(l/g)
n=floor((10.0*tau)/h)

allocate (x(n), v(n), a(n), k(n), u(n), e(n), t(n))

t(1)=0.0
write(*,*) 'theta initial ='
read(*,*) x(1)
x(1)=(x(1)*pi)/180
v(1)=0.0
a(1)=(-g/l)*sin(x(1))
k(1)=(0.5)*m*(v(1)**2)
u(1)=m*g*l*(1-cos(x(1)))
e(1)=k(1)+u(1)

x(2)=x(1)+(h*v(1))+(0.5)*(h**2)*a(1)
t(2)=0.1

do i= 2,(n-1)

a(i)=(-g/l)*sin(x(i))
x(i+1)=(2.0*x(i))-x(i-1)+((h**2.0)*a(i))
v(i)=(x(i+1)-x(i-1))/(h*2.0)

write(110,*) t(i), x(i)
write(210,*) t(i), v(i)
write(310,*) t(i), a(i)

u(i)=m*g*l*(1-cos(x(i)))
k(i)=(0.5)*m*v(i)**2.0
e(i)=k(i)+u(i)

write(410,*) t(i), u(i)
write(510,*) t(i), k(i)
write(610,*) t(i), e(i)
t(i+1)=t(i)+h

enddo

end program shmpendula
