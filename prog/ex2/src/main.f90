!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Programa principal:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use :: constantes
use :: vector_sub
use :: sim_mod
implicit none

! Vars
! real(dp) :: tr
! real(dp), dimension(3) :: r = [1._dp,40._dp,10._dp]
! logical  :: sol

! Execução:
call read_const()

print 100, x0, v0,a, h, t0
100 format (/,"x0: ",ES12.5,/,"v0: ", ES12.5,/,"a: ",ES12.5,/,"h: ",ES12.5,/,"t0: ",ES12.5,/)

! Se necessário, execute a linha abaixo com valores de r sobre o eixo x.
! call comparison([12._dp,0*h,0*h], t0)

call write_data()
! call find_root(r,t0,tr,sol)
! print*, tr, sol


print("(/,'Dados salvos. Fim de execução.')")

end program main