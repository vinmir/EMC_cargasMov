!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Módulo de subrotinas de simulação
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module sim_mod
use :: constantes
use :: vector_sub
implicit none

! Declaração de variáveis do módulo:

! Declaração de procedures do módulo:
contains
!**********************************************************************************************************************************!
subroutine sub_test()
!! Subrotina de teste

! IN

! OUT:

! Local:

! Execução:
    print*, ""

end subroutine 
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
function root_func(r,t,tr) result(f)
!! Função auxiliar ao algoritmo de cálculo do tempo atrasado.
!! r: posição da malha
!! t: instante de tempo atual
!! tr: chute para tr

! IN
real(dp), dimension(3), intent(in) :: r
real(dp), intent(in) :: t, tr
! OUT:
real(dp) :: f
! Local:
real(dp), dimension(3) :: w
! Inicialização:
w = (/ v0*tr+x0, 0._dp, 0._dp /) ! Inicializo o vetor w no tempo tr indicado

! Execução:
f = norm(r-w) - 1*(t-tr)

end function
!**********************************************************************************************************************************!


!**********************************************************************************************************************************!
function find_root(r,t) result(tr_final)
!! Função que calcula o valor de tr para um dado r no instante de tempo t0.
!! r: pos. da malha
!! t: instante de tempo
! A fução aparenta funcionar perfeitamente. Comparei com valores x0,v0 -> 0, nos quais o Griffiths fornece que
! tr = t - r/c -> t - r, já que c = 1. Os resultados foram coerentes.

! IN:
real(dp), dimension(:), intent(in) :: r ! Ponto da malha a ser analisado.
real(dp), intent(in) :: t

! OUT:
real(dp) :: tr_final ! Tempo atrasado.

! Local:
real(dp) :: tr_a,tr_b,tr_m ! Variáveis para a bisecção
real(dp) :: dt ! Espaçamento inicial.

! Execução:
! Inicialmente, tr <= t (). Logo, um chute inicial para tr_b é t.
tr_b = t

! Primeiro, vejamos se tr = t (isto acontece quando r = w, claro):
if (abs(root_func(r,t,t)) < 1e-5_dp) then
    tr_final = t
else
    ! t não é raiz. Nesse caso, continuemos:
    ! O procedimento inicial é encontrar a faixa onde a raiz está. Portanto, recuarei t_a em intervalos de t0/5 até determinar
    ! a faixa de tr onde a função tem uma raiz:
    
    dt = t0/5
    tr_a = tr_b-dt

    do
        if (root_func(r,t,tr_a)*root_func(r,t,tr_b) < 0._dp) then
            exit
        end if

        tr_a = tr_a-dt/5

    end do

    ! print("(A,g0,A,g0)"), "tr_a: ",tr_a, "tr_b: ", tr_b
    ! Encontrei a faixa correta para a raiz. Agora, é simples bissecção:
    do
        tr_m = (tr_a+tr_b)/2._dp
        ! Vejamos se a precisão para tr_m é boa o suficiente:
        if (abs(root_func(r,t,tr_m)) < 1.e-5_dp ) then
            tr_final = tr_m
            ! print*, "tr_final = ",tr_final
            ! print*, "func(t,tr_final) = ", root_func(r,t,tr_m) 
            exit
        end if

        if (root_func(r,t,tr_m)*root_func(r,t,tr_b) < 0._dp) then
            ! Neste caso, a raiz está entre tr_m e tr_b
            tr_a = tr_m
        else
            ! Neste caso, a raiz está entre tr_a e tr_m
            tr_b = tr_m
        end if



    end do
end if

end function
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
subroutine calc_fields(r,t,E,B)
!! Subrotina que calcula os campos E(r,t) e B(r,t).

! IN
real(dp), dimension(3), intent(in) :: r
real(dp) :: t

! OUT:
real(dp), dimension(3), intent(out) :: E,B

! Local:
real(dp) :: tr
real(dp), dimension(3) :: delta_r, w, u ! w é o w(tr); delta_r = r-w(tr); u = 1*r_versor - (v0,0,0)

! Necessários para testar tr_ana:
! real(dp), dimension(3) :: new_r
! real(dp) :: tr_ana

! Inicialização:
! Primeiro, calculemos o tempo retardado associado ao par (r,t):
tr = find_root(r,t)
! Em seguida, definamos o vetor da posição da carga no tempo retardado:
w = (/x0+v0*tr,0._dp,0._dp /)
! Finalmente, resta definir u e delta_r
delta_r = r - w
u = normalize(delta_r) - (/ v0,0._dp,0._dp /)

! Execução:
E = norm(delta_r)/(dot_product(delta_r,u))**3 * (1-v0**2)*u
B = cross(normalize(delta_r),E)

! Teste para o tr:
! new_r = r - [x0,0._dp,0._dp]
! tr_ana = ((t-new_r(1)*v0) - sqrt((t-new_r(1)*v0)**2+(1-v0**2)*(norm(new_r)**2-t**2)))/(1-v0**2)
! if (abs(tr-tr_ana) > 1.e-3_dp) then
! print*, "t = ", t, "tr = ", tr, "tr_ana = ", tr_ana
! stop
! end if

end subroutine 
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
subroutine write_data()
!! Subrotina que escreve E e B em uma discretização no espaço.

! Local:
real(dp) :: t ! Instante real de tempo.
real(dp), dimension(3) :: r, E, B ! Posição real na malha e campos.
integer :: i,j,k ! Índices de loop.
character(len=*), parameter :: path_E="./out/E_field.dat", path_B="./out/B_field.dat" ! Caminhos dos arquivos de saída
integer, parameter :: unit_E=10, unit_B=11

! Inicialização:
t = t0 ! Todos os campos serão avaliados em t0.

! Execução:
open(unit=unit_E,file=path_E,status="replace",action="write")
write(unit_E,*) "x    y    z    Ex    Ey    Ez"
open(unit=unit_B,file=path_B,status="replace",action="write")
write(unit_B,*) "x    y    z    Bx    By    Bz"

! Loop de itreação para os pontos da malha. Lembre-se: m_pontos é definido no módulo de constantes.
do i=-m_pontos,m_pontos
    do j=-m_pontos,m_pontos
        do k=-m_pontos,m_pontos
            r = (/ i*h, j*h, k*h /)
            call calc_fields(r,t,E,B)
            write(unit_E,*) r, E
            write(unit_B,*) r, B
        end do
    end do
end do


end subroutine 
!**********************************************************************************************************************************!



end module
