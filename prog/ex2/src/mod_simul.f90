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
w = (/ a*tr**2/2._dp + v0*tr + x0, 0._dp, 0._dp /) ! Inicializo o vetor w no tempo tr indicado

! Execução:
f = norm(r-w) - 1*(t-tr)

end function
!**********************************************************************************************************************************!


!**********************************************************************************************************************************!
subroutine find_root(r,t,tr_final,solution)
!! Função que calcula o valor de tr para um dado r no instante de tempo t0.
!! r: pos. da malha
!! t: instante de tempo
! No caso do MRUV, não há garantia de haver solução real para (r,t) arbitrários, pois a equação é em potência 4 em tr.
! Existirão valores de r para os quais a solução é complexa. Isto é diretamente associado ao valor da velocidade
! da carga em tr, que deve ser estritamente menor do que 1 em módulo.
! Como o tempo será negativo, a velocidade com o retrocesso em t_a aumenta em +x, até atingir o máximo teórico (v(t_a_min) = +1).
! O t_a_min é, portanto, +1 = v0 +a*t_a_min, 
!   t_a_min = (1-v0)/a
! Se o algoritmo não encontrar uma solução t_a antes de atingir t_a_min, nunca encontrará.
! Precisei converter de função para subrotina. Do contrário, não teria como determinar que
! E = B = 0 em (r,t).

! IN:
real(dp), dimension(3), intent(in) :: r ! Ponto da malha a ser analisado.
real(dp), intent(in) :: t

! OUT:
real(dp), intent(out) :: tr_final ! Tempo atrasado.
logical, intent(out)  :: solution ! Determina se houve solução.

! Local:
real(dp) :: tr_a,tr_b,tr_m ! Variáveis para a bisecção
real(dp) :: dt ! Espaçamento inicial.
real(dp) :: t_a_min ! Veja informações no topo.

! Inicialização:
solution = .true. ! Deve iniciar como verdadeiro. Será falso em caso de não convergência.
tr_final = 0._dp  ! Valor qualquer apenas para inicializar tr_final.

! Execução:
! Inicialmente, tr <= t (). Logo, um chute inicial para tr_b é t.
tr_b = t
t_a_min = (1-v0)/a

! Primeiro, vejamos se tr = t (isto acontece quando r = w, claro):
if (abs(root_func(r,t,t)) < 1e-5_dp) then
    tr_final = t
else
    ! t não é raiz. Nesse caso, continuemos:
    ! O procedimento inicial é encontrar a faixa onde a raiz está. Portanto, recuarei t_a em intervalos de t0/5 até determinar
    ! a faixa de tr onde a função tem uma raiz:
    
    dt = t0/5
    tr_a = tr_b-dt

    tr_a_loop: do
        if (root_func(r,t,tr_a)*root_func(r,t,tr_b) < 0._dp) then
            exit tr_a_loop
        end if

        tr_a = tr_a-dt
        if (tr_a < t_a_min) then
            ! print*, "Solução impossível."
            solution = .false.
            exit tr_a_loop
        end if 

    end do tr_a_loop

    ! Encontrei a faixa correta para a raiz. Agora, é simples bissecção, contanto que haja t_a.
    tr_final_loop: do
        ! Verifica se uma solução é possível:
        if (solution .eqv. .false.) then
            exit tr_final_loop
        end if 
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



    end do tr_final_loop
end if

end subroutine
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
subroutine calc_fields(r,t,E,B)
!! Subrotina que calcula os campos E(r,t) e B(r,t).
! Utilizei o exercício 10.20 do Griffiths (4ª ed.) para confirmar. Sim, a simulação aparenta estar correta.
! Criarei uma subrotina melhor para isto.

! IN
real(dp), dimension(3), intent(in) :: r
real(dp), intent(in) :: t

! OUT:
real(dp), dimension(3), intent(out) :: E,B

! Local:
real(dp) :: tr, v_tr ! Explicativas
logical :: solution  ! Utilizado na subrotina de solução
real(dp), dimension(3) :: delta_r, w, u ! w é o w(tr); delta_r = r-w(tr); u = 1*r_versor - (v0,0,0)

! Necessários para testar tr_ana:
! real(dp), dimension(3) :: new_r
! real(dp) :: tr_ana

! Inicialização:
! Primeiro, calculemos o tempo retardado associado ao par (r,t):
call find_root(r,t,tr,solution)
! Em seguida, definamos o vetor da posição da carga no tempo retardado, assim como v_tr:
w = (/x0+v0*tr+a*tr**2/2, 0._dp, 0._dp /)
v_tr = v0 + a*tr
! Finalmente, resta definir u e delta_r
delta_r = r - w
u = normalize(delta_r) - (/ v0 + a*tr,0._dp,0._dp /) ! Faltava atualizar isto.

! Execução:
if (solution .eqv. .false.) then
    E = 0._dp
    B = 0._dp
    ! print*, "Solução impossível."
else
    ! Obviamente, faltava completar a expressão de E. Ademais, a velocidade aqui estava errada.
    E = norm(delta_r)/(dot_product(delta_r,u))**3 *( (1-v_tr**2)*u + cross(delta_r,cross(u,(/a,0._dp,0._dp/))))
    B = cross(normalize(delta_r),E)
end if




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

! integer :: sol=0, nsol=0

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
            ! read(*,*)
            ! print("('i=',i0,' j=',i0,' k=',i0)"),i,j,k 
            r = (/ i*h, j*h, k*h /)
            call calc_fields(r,t,E,B)
            write(unit_E,*) r, E
            write(unit_B,*) r, B

        end do
    end do
end do


end subroutine 
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
subroutine comparison(r,t)
!! Subrotina que compara a solução numérica do campo elétrico com a solução analítica
!! do problema 10.20, 4ª ed. do Griffiths.

! IN:
real(dp), dimension(3), intent(in) :: r
real(dp) :: t

! Local:
real(dp), dimension(3) :: E_n, B_n, w, delta_r
logical :: sol
real(dp) :: tr, E_ana, v_tr

! Execução:
call find_root(r,t,tr,sol)

if (sol .eqv. .true.) then
    w = [a*tr**2/2 + v0*tr + x0,0*h,0*h]
    v_tr = v0 + a*tr
    ! print*, v_tr
    delta_r = r-w
    call calc_fields(r,t,E_n,B_n)
    if (delta_r(1)*v_tr > 0) then
        print*, "Paralelo"
        E_ana = -1/norm(delta_r)**2*(1+abs(v_tr))/(1-abs(v_tr)) ! Sinal trocado por causa da re-orientação do eixo
                                                                ! com relação ao exercício 10.20.
    else
        print*, "Anti-paralelo"
        E_ana = +1/norm(delta_r)**2*(1-abs(v_tr))/(1+abs(v_tr)) ! Leia obs. acima.
    end if
    print*, "E_n = ", E_n(1)
    print*, "E_ana = ", E_ana
else
    print*, "Solução impossível."
end if

end subroutine 
!**********************************************************************************************************************************!


end module
