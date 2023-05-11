!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Módulo constantes:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module constantes
! use :: iso_fortran_env, only: sp=>real32, dp=>real64 !!!! Caso esteja em Fortran 2008+, utilize isto para sp, dp.
implicit none

! Declaração de variáveis do módulo:
    ! Variáveis que denotam precisão simples e dupla pelo método <=F90:
integer, parameter :: sp = selected_real_kind(6,37)    
integer, parameter :: dp = selected_real_kind(15,307)
    ! Variáveis de uso geral:
real(dp), protected :: v0 ! Velocidade inicial (v0x) da carga
real(dp), protected :: x0 ! Posição inicial da carga
real(dp), protected :: h  ! Discretização da malha 3D. Cada direção tem, por definição, 20*x0 de comprimento.
real(dp), protected :: t0 ! Instante em que a carga está na origem.
integer, parameter :: m_pontos=10 ! A quantidade de pontos da malha, em cada direção, é N = 2*m+1 

! Declaração de rotinas do módulo (leitura das constantes, etc.):
contains
!**********************************************************************************************************************************!
subroutine read_const()
!! Subrotina que faz a leitura das constantes de uso geral.

! Local
integer, parameter :: N = 2*m_pontos+1 ! Determina a quantidade de pontos em cada direção da malha.

! Execução
print("(A)"), "Projeto 4, Exercício 1."
print("(A)"), "Insira o valor de v0 (-1 < v0 < 0):"
read(*,*) v0
print("(A)"), "Insira o valor de x0 (x0 > 0):"
read(*,*) x0

if (v0 < -1 .or. v0 > 0 .or. x0 < 0) then
    print("(A,//)"), "Constantes incompatíveis. Reajustando v0=-0.8, x0=2."
    v0=-0.8_dp
    x0=2._dp
    ! stop
end if

h = 20._dp*x0/(N-1)
t0 = -x0/v0

end subroutine
!**********************************************************************************************************************************!

end module constantes