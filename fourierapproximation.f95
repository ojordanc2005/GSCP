program fourierapproximation

implicit none

real(8):: x, y, pi
integer(4):: n
pi = 4.0*atan(1.0)

open(unit=200, file='fourierapproximation.dat')

write(200,*)   '# x-values       y-values'

x=0.0

do while (x<=2.0*pi)

y=0.0

do n=1,1000
y=y+((((-1.0)**(n+1.0))*2.0)/n)*sin(n*x)
enddo

write(200,*) x, y
x=x+0.01

enddo

end program fourierapproximation



