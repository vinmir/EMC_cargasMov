!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Programa principal:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use :: constantes
use :: vector_sub
use :: sim_mod
implicit none

! Execução:
call read_const()

print 100, x0, v0, h, t0
100 format (/,"x0: ",ES12.5,/,"v0: ", ES12.5,/,"h: ",ES12.5,/,"t0: ",ES12.5,/)

call write_data()
print("(/,'Dados salvos. Fim de execução.')")

end program main