program orbitb
implicit none

real(8):: state(6), r(3), v(3), a(3), k1(6), k2(6), k3(6), k4(6), temp(6)
real(8):: K, P, E, x, y, z, t, h, pi, GM, m1
integer(4):: i

open(unit=100, file='orbitkvtb.dat')
open(unit=200, file='orbitpvtb.dat')
open(unit=300, file='orbitevtb.dat')
open(unit=400, file='orbitxyzb.dat')

pi= 4*atan(1.0)
h= 0.01
GM= 4*(pi**2.0)
m1= 1.0
r=0
r(1)=1.0
v=0
v(1)=1.0
v(2)=2.0*pi
state(1:3)= r
state(4:6)= v
do i=1,1000

!k1= d/dx(r0,v0)

a= -GM/((norm2(r))**3)*r
k1(1:3)= v
k1(4:6)= a

!k2=d/dx[r0,v0]+h/2*k
temp=state+(0.5*h*k1)
a= -GM/((norm2(temp(1:3)))**3)*temp(1:3)
k2(1:3)= temp(4:6)
k2(4:6)= a

!k3=d/dx[r0,v0]+h/2*k
temp=state+(0.5*h*k2)
a= -GM/((norm2(temp(1:3)))**3)*temp(1:3)
k3(1:3)= temp(4:6)
k3(4:6)= a

!k4=d/dx[r0,v0]+h*k
temp=state+(h*k3)
a= -GM/((norm2(temp(1:3)))**3)*temp(1:3)
k4(1:3)= temp(4:6)
k4(4:6)= a

!r(1:3)= r+((h/6)*(k1+(2*k2)+(2*k3)+k4))
state=state+((h/6)*(k1+(2*k2)+(2*k3)+k4))

!RK4: y1=y0+(h/6)*(k1+(2*k2)+(2*k3)+k4)


r=state(1:3)
v=state(4:6)

K=(0.5)*m1*((norm2(v))**2)
P=-(GM*m1)/(norm2(r)) 
E=K+P
a= -GM/((norm2(r))**3)*r
write(100,*) t, K
write(200,*) t, P
write(300,*) t, E
write(400,*) r(1),r(2),r(3)

t=t+h

enddo
 
endprogram orbitb
