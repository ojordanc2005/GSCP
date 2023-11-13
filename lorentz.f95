program lorentz
implicit none

real(8):: state(6), r(3), v(3), a(3), E(3), B(3), k1(6), k2(6), k3(6), k4(6), temp(6), Fe(3), Bv(3)
real(8):: x, y, z, t, h, pi, m, q
integer(4):: i

open(unit=100, file='lorentzxyz.dat')

pi= 4*atan(1.0)
h= 0.01
m= 1.0
q=1.0
r=0
v=0
v(1)=1.0
E=0
E(2)=0.01
E(3)=0.1
B=0
B(3)=0.1
a(1)= (q/m)*(E(1)+((v(2)*B(3))-(v(3)*B(2))))
a(2)= (q/m)*(E(2)-((v(1)*B(3))-(v(3)*B(1))))
a(3)= (q/m)*(E(3)+((v(1)*B(2))-(v(2)*B(1))))

!state(1:3)= r
!state(4:6)= v
do i=1,100000

!k1= d/dx(r0,v0)

!Fe=q*norm2(E)
!Fb=q*norm2(Bv) 
state(1:3)= r
state(4:6)= v
a(1)= (q/m)*(E(1)+((v(2)*B(3))-(v(3)*B(2))))
a(2)= (q/m)*(E(2)-((v(1)*B(3))-(v(3)*B(1))))
a(3)= (q/m)*(E(3)+((v(1)*B(2))-(v(2)*B(1))))
k1(1:3)= v
k1(4:6)= a

!k2=d/dx[r0,v0]+h/2*k
temp=state+(0.5*h*k1)
r= temp(1:3)
v= temp(4:6)
a(1)= (q/m)*(E(1)+((v(2)*B(3))-(v(3)*B(2))))
a(2)= (q/m)*(E(2)-((v(1)*B(3))-(v(3)*B(1))))
a(3)= (q/m)*(E(3)+((v(1)*B(2))-(v(2)*B(1))))
k2(1:3)= temp(4:6)
k2(4:6)= a

!k3=d/dx[r0,v0]+h/2*k
temp=state+(0.5*h*k2)
r= temp(1:3)
v= temp(4:6)
a(1)= (q/m)*(E(1)+((v(2)*B(3))-(v(3)*B(2))))
a(2)= (q/m)*(E(2)-((v(1)*B(3))-(v(3)*B(1))))
a(3)= (q/m)*(E(3)+((v(1)*B(2))-(v(2)*B(1))))
k3(1:3)= temp(4:6)
k3(4:6)= a

!k4=d/dx[r0,v0]+h*k
temp=state+(h*k3)
r= temp(1:3)
v= temp(4:6)
a(1)= (q/m)*(E(1)+((v(2)*B(3))-(v(3)*B(2))))
a(2)= (q/m)*(E(2)-((v(1)*B(3))-(v(3)*B(1))))
a(3)= (q/m)*(E(3)+((v(1)*B(2))-(v(2)*B(1))))
k4(1:3)= temp(4:6)
k4(4:6)= a

!r(1:3)= r+((h/6)*(k1+(2*k2)+(2*k3)+k4))
state=state+((h/6)*(k1+(2*k2)+(2*k3)+k4))

!RK4: y1=y0+(h/6)*(k1+(2*k2)+(2*k3)+k4)

r(1:3)=state(1:3)
v(1:3)=state(4:6)

write(100,*) r(1),r(2),r(3)

t=t+h

enddo

end program lorentz
