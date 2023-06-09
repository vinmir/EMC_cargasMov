!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Módulo de subrotinas de vetores:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module vector_sub
use :: constantes
implicit none

! Declaração de variáveis do módulo:

! Declaração de procedures do módulo:
contains
!**********************************************************************************************************************************!
function cross(vec_a,vec_b) result(c)
!! Função que calcula o produto vetorial de dois vetores.

! IN
real(dp), dimension(3), intent(in) :: vec_a,vec_b
! OUT:
real(dp), dimension(3) :: c

! Execução:
c(1) = vec_a(2)*vec_b(3)-vec_a(3)*vec_b(2)
c(2) = vec_a(3)*vec_b(1)-vec_a(1)*vec_b(3)
c(3) = vec_a(1)*vec_b(2)-vec_a(2)*vec_b(1)

end function 
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
function norm(vec)
!! Função que calcula a norma de um vetor.

! IN
real(dp), dimension(3), intent(in) :: vec
! OUT:
real(dp) :: norm

! Execução:
norm = sqrt(dot_product(vec,vec))

end function 
!**********************************************************************************************************************************!


!**********************************************************************************************************************************!
function normalize(vec)
!! Função que normaliza um vetor.

! IN
real(dp), dimension(3), intent(in) :: vec
! OUT:
real(dp), dimension(3) :: normalize

! Execução:
normalize = vec/norm(vec)

end function 
!**********************************************************************************************************************************!

end module
