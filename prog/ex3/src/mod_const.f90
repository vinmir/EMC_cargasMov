!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Módulo constantes:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module zmod_const
! use :: iso_fortran_env, only: sp=>real32, dp=>real64 !!!! Caso esteja em Fortran 2008+, utilize isto para sp, dp.
implicit none

! Declaração de variáveis do módulo:
    ! Variáveis que denotam precisão simples e dupla pelo método <=F90:
integer, parameter :: sp = selected_real_kind(6,37)    
integer, parameter :: dp = selected_real_kind(15,307)
    ! Variáveis de uso geral:
real(dp), protected :: omega ! Omega será acessível a todos os programas. Atenção! Omega deverá variar para o item 3b).
real(dp), protected :: A     ! Constante de amplitude de w(t).
real(dp), protected :: periodo ! Período de oscilação.
real(dp), parameter :: pi=4*atan(1._dp) ! pi
! Declaração de rotinas do módulo (leitura das constantes, etc.):
contains
!**********************************************************************************************************************************!
subroutine read_const(new_omega)
!! Subrotina que faz a leitura das constantes de uso geral.

! Local
real(dp), intent(in) :: new_omega

! Inicialização das variáveis
omega = new_omega ! Ajusta o novo omega da simulação.

! Define A tal que, conforme dita o enunciado, A << c/omega, ou A = (1/omega)*10^(-6):
! A = (1.0e-6_dp/omega)
A = 1.e-6_dp
! Define o período:
periodo = 2*pi/omega

end subroutine
!**********************************************************************************************************************************!

end module