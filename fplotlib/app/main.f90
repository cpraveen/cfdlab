program main
use fplotlib
implicit none
integer,parameter :: n = 500
integer :: i
double precision :: pi = 4.0d0 * atan(1.0d0)
double precision :: dx, x(n), y(n), y2(n)

dx = pi / (n - 1)
do i=1,n
   x(i) = i*dx
   y(i) = sin(x(i))
   y2(i) = cos(x(i))
enddo

call plot(x, y, "k-", label="sin")
call title("My plot")
call xlabel("x")
call ylabel("y")
call grid(.true.)
call legend()
call savefig("out1.svg")   ! or "out.png", "out.pdf", "out.eps"
call show()                ! Jupyter (LFortran) or writes fplotlib_show.svg

! subplots
call clf()
call subplot(2, 1, 1)
call plot(x, y, "b-", label="sin")
call title("top")
call subplot(2, 1, 2)
call plot(x, y2, "r--", label="cos")
call title("bottom")
call suptitle("figure title")
call savefig("out2.svg")   ! or "out.png", "out.pdf", "out.eps"
call show()                ! Jupyter (LFortran) or writes fplotlib_show.svg
end program
