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
real(dp), protected :: a  ! Aceleração (ax) da carga
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
real(dp) :: v_a ! Velocidade da carga quando esta passa pela origem; v_a = v(t0).

! Execução
! print("(A)"), "Projeto 4, Exercício 2: MRUV."
! print("(A)"), "Insira o valor de x0 (x0 > 0):"
! read(*,*) x0
! print("(A)"), "Insira o valor de a (-0.1 < a < 0):"
! read(*,*) a
! print("(A)"), "Insira o valor de v_a, vel. da carga em w(t0)=0 (-1 < v0 < 0):"
! read(*,*) v_a


! if (x0 < 0 .or. a < -0.1 .or. a > 0) then
    print("(/,A)"), "Constantes fixadas no módulo `mod_const`, por conveniência. x0=2.0, a=-0.05, v(t0)=-0.8"
    x0=2._dp
    a=-0.05_dp
    v_a=-0.8
! end if

! Define o espaçamento da malha:
h = 20._dp*x0/(N-1)

! Resta calcular o instante t0 onde a carga está na origem. Verifico também se t0 é real.
! Fixados v_a, a e x0, as expressões para t0 e v0 são:
! t0 = (-v0 - sqrt(v0**2 - 2*a*x0))/a (t0 > 0 e a<0, então a solução com subtração é a que fornece t0 > 0)
! v_a = v0 + a*t0 ou, por Torricelli, v_a^2 = v0^2 + 2*a*(0-x0) => v0 = -sqrt(v_a**2 + 2*a*x0). v0<0 por restrição.
! Finalmente, t0 = (v_a-v_0)/a.

if (v_a**2 + 2*a*x0 < 0._dp) then
    print*, "Parâmetros geram t0 complexo. Encerrando."
    stop
end if

v0 = -sqrt(v_a**2 + 2*a*x0)
t0 = (v_a-v0)/a

end subroutine
!**********************************************************************************************************************************!

end module constantes