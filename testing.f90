program testing

use stafft

integer:: n,factors(5),ierr, i
double precision, allocatable:: trig(:)

n = 128

allocate(trig(2*n))

call initfft(n, factors, trig)

write(*, "(f13.10)") trig
  
end program
