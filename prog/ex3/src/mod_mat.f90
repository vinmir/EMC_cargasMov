!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Módulo de rotinas matemáticas:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module zmod_mat
use :: zmod_const
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

!**********************************************************************************************************************************!
function simpson_int(func,ta,tb) result(integral)
!! Função que calcula a integral de uma função real em um intervalo [a,b] pelo método de Simpson.
!! A função poderá ter qualquer dimensão. Caso esteja usando o fortls, ou Fortran IntelliSense, saiba
!! que há um falso-problema. A interface funciona corretamente e o programa é devidamente compilado.

! IN:
real(dp), intent(in) :: ta,tb ! Limites de integração.
! real(dp), external   :: func  ! Função f(t) a ser integrada. Pode ser vetor.
interface
  function func(t) result(f)
    use zmod_const
    real(dp), intent(in) :: t
    real(dp), allocatable, dimension(:) :: f      ! Use isto, é muito superior. Assumed-shape arrays são o futuro.
    ! real(dp), dimension(3) :: f       ! Este é o método antigo
  end function
end interface

! OUT:
real(dp), allocatable, dimension(:) :: integral ! A integral final.

! Local:
real(dp) :: dt ! Espaçamento dos pontos
real(dp) :: t_i ! Variáveis temporárias da função, ip=i+1
integer, parameter  :: N=100  ! Pontos (f_0, f_1, ..., f_N): Total de pontos = N+1.
                                    ! Para o método de Simpson, o total de regiões deve ser PAR.
                                    ! Como (Regiões) = (Pontos-1) = (N+1-1) = N, N deve ser par.
integer  :: i  ! Índice de loop

! Inicialização:
dt = (tb-ta)/N ! Dado dt, t_i = a + dt*i, i=1,...,N
integral = 0._dp*func(ta) ! Não sei a quantidade de elementos em func(t). Esta é uma forma de inicializar o array.


! Execução:
! O algoritmo do trapézio consta na bibliografia.
! Em resumo: integral = (dt/3)(f_0+f_N) + (2*dt/3)*(f_2+f_4+...+f_(N-2)) + (4*dt/3)*(f_1+f_3+...+f_(N-1))
    ! Lembre-se: f(a) = f_0; f(b) = f_N
integral = integral + (dt/3._dp)*(func(ta)+func(tb))  ! Primeira parte da integração
do i=2,N-2,2 ! Segunda parte da integração
    t_i = ta+i*dt
    integral = integral + (2._dp*dt/3._dp)*func(t_i)
    ! print*, func(t_i)
end do
do i=1,N-1,2 ! Terceira parte da integração
    t_i = ta+i*dt
    integral = integral + (4._dp*dt/3._dp)*func(t_i)
end do

end function 
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
function simpson_array_int(f_vec,dt) result(integral)
!! Implimenta a integração de Simpson mas com um array f=[f_0,f_1,...,f_N] e um intervalo fixo de discretização dt.

! IN
real(dp), intent(in) :: dt
real(dp), dimension(0:), intent(in) :: f_vec ! O array terá bounds (0:N) ao invés de (1,N+1).
! OUT
real(dp) :: integral
! Local
integer :: i ! Índice de loop
integer :: N ! Total de pontos = (N+1)

! Inicialização:
integral=0._dp
N = size(f_vec) - 1
if (mod(N,2) == 1) then
    print("(A)"), "Array de tamanho par de pontos (equiv.: tamanho ímpar de seções) &
    & entrou na rotina de integração de Simpson. Encerrando o programa."
    stop
end if

! Execução:
! Primeira parte da integração:
integral = integral + (dt/3._dp)*(f_vec(0)+f_vec(N))
! Segunda Parte da integração:
do i=2,N-1,2
    integral = integral + (2._dp*dt/3._dp)*f_vec(i)
end do
! Terceira parte da integração:
do i=1,N-1,2 ! Terceira parte da integração
    integral = integral + (4._dp*dt/3._dp)*f_vec(i)
end do

end function
!**********************************************************************************************************************************!



end module
