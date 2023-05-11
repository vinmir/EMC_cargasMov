!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Módulo principal (subrotinas):
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module main_mod
use :: constantes
implicit none

! Declaração de variáveis do módulo:

! Declaração de procedures do módulo:
contains
!**********************************************************************************************************************************!
subroutine main_sub(x)
!! Subrotina de teste

! Declaração de variáveis:
real(dp), dimension(3), intent(out) :: x
real(dp), dimension(0:10,3) :: y
integer :: i

! Inicialização:
! Execução:
do i=0,10
    y(i,:) = x*i
end do
do i=1,3
    x(i) = sum(y(:,i))
end do

end subroutine 
!**********************************************************************************************************************************!



end module main_mod