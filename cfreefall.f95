program cfreefall
implicit none
real(8):: y, v, t, a, h, g
real(8):: ys, vs, ts, as
real (8):: k, m

open(unit=120, file='yvtc.dat')
write(120,*) '#t-values, y-values'
open(unit=220, file='vvtc.dat')
write(220,*) '#t-values, v-values'
open(unit=320, file='avtc.dat')
write(320,*) '#t-values, a-values'

h= 0.001
y= 57.0
v= -3.8
t= 0.0
g= -9.8
k= 0.00132 !k= 0.5*dc*csa*p
m= .145

do while (y>=0)

vs=v
v=v+a*h
a=((k*(v**2))/m)+g
y=y+((vs+v)/2)*h

write(120,*) t, y
write(220,*) t, v
write(320,*) t, a
t=t+h

enddo
endprogram cfreefall



