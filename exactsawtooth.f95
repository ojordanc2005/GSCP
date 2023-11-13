program exactsawtooth

implicit none

real(8):: y, x, pi
integer(4):: i

pi = 4*atan(1.0)
x=0.0

open(unit=100, file='exactsawtooth.dat')
write(100,*) "# x-values        y-values"

do while (x<=2.0*pi)
if(x<=pi) then
y=x
elseif (x>pi) then
y=x-(2.0*pi)

endif

x=x+0.01


write(100,*) x,y

enddo

close(100)

end program exactsawtooth
