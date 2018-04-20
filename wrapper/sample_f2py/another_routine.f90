subroutine test_routine

! to test defining arrays using variables from other modules

  use another_mod
  
  real :: a(n)
  
  print *, n
  
  a(:) = 5

end subroutine
