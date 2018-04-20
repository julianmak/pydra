module test_mod

contains

! sample print

subroutine print_murp(n)

real, intent(in) :: n

print *, "call murp return: ", n

end subroutine

!!!!!!!!!!!!! sample sum

subroutine add_murp(n, m, p)

real :: n, m, p

! either: write the intent into the declaration in fortran (inout -> in,out)
!     or: add in a comment as below (inout -> in,out) to declare it as appropriate
! the latter seems to be safer

!f2py intent(in)  :: n, m
!f2py intent(out) :: p

p = n + m

call print_murp(p)

end subroutine

!!!!!!!!!!!!! sample add to itself

subroutine self_murp(n)

real :: n

!f2py intent(in, out) :: n

n = n + n

call print_murp(n)

end subroutine

!!!!!!!!!!!!! sample to have multiple outputs

subroutine multiple_murp(n, m, p, q)

real :: n, m, p

!f2py intent(in)  :: n, m
!f2py intent(out) :: p, q

p = n + m
q = n * m

call print_murp(p)
call print_murp(q)

end subroutine

end module
