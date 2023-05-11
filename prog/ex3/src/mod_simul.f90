!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Módulo de subrotinas de simulação
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module zmod_simul
use :: zmod_const
use :: zmod_mat
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
w = pos_q(tr) ! Inicializo o vetor w no tempo tr indicado

! Execução:
f = norm(r-w) - 1*(t-tr)

end function
!**********************************************************************************************************************************!


!**********************************************************************************************************************************!
function find_root(r,t) result(tr_final)
!! Função que calcula o valor de tr para um dado r no instante de tempo t0.
!! r: pos. da malha
!! t: instante de tempo
! Criei uma função numérica no Mathematica para calcular tr (NSolve). As soluções são idênticas, então este método
! aparenta funciona corretamente.

! IN:
real(dp), dimension(3), intent(in) :: r ! Ponto da malha a ser analisado.
real(dp), intent(in) :: t

! OUT:
real(dp) :: tr_final ! Tempo atrasado.

! Local:
real(dp) :: tr_a,tr_b,tr_m ! Variáveis para a bisecção
real(dp) :: dt ! Espaçamento da rotina inicial de determinação de [t_a, t_b].
! real(dp), dimension(3) :: vel_t ! Determina a velocidade da carga no instante de teste. Se for superior a 1,
!                                 ! não haverá convergência

! Inicialização:
tr_final = 0._dp  ! Valor qualquer apenas para inicializar tr_final.

! Execução:
! Inicialmente, tr <= t (). Logo, um chute inicial para tr_b é t.
tr_b = t

! Primeiro, vejamos se tr = t (isto acontece quando r = w, claro):
if (abs(root_func(r,t,t)) < 1e-7_dp) then
    tr_final = t
else
    ! t não é raiz. Nesse caso, continuemos:
    ! O procedimento inicial é encontrar a faixa onde a raiz está. Portanto, recuarei t_a em intervalos de T até determinar
    ! a faixa de tr onde a função tem uma raiz, onde T é um período de oscilação.
    ! O motivo para um intervalo grande de recuo? Fiz análises numéricas pelo Mathematica, e aqui há apenas uma solução
    ! real ao invés de (potencialmente) duas. Um recuo de T/20 é extremamente custoso e demora a simular.
    ! Diferentemente dos outros casos, aqui vale mais a pena usar um chute inicial para t_a.
    ! Como |r-w(tr)| ~ |r|, dado que w = A cos(), A << 1, é só usar |r| = (t-t_a), ou t_a = t - |r|.
    
    dt = 1*periodo
    tr_a = t - norm(r)

    tr_a_loop: do
        if (root_func(r,t,tr_a)*root_func(r,t,tr_b) < 0._dp) then
            exit tr_a_loop
        end if

        tr_a = tr_a-dt

    end do tr_a_loop

    ! Encontrei a faixa correta para a raiz. Agora, é simples bissecção, contanto que haja t_a.
    tr_final_loop: do
        tr_m = (tr_a+tr_b)/2._dp
        ! Vejamos se a precisão para tr_m é boa o suficiente:
        if (abs(root_func(r,t,tr_m)) < 1.e-10_dp ) then
            ! print*, abs(root_func(r,t,tr_m))
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

end function
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
subroutine calc_fields(r,t,E_rad,S_rad)
!! Subrotina que calcula os campos E_rad(r,t) e S_rad(r,t).

! IN
real(dp), dimension(3), intent(in) :: r ! Posição r
real(dp), intent(in) :: t ! Instante t

! OUT:
real(dp), dimension(3), intent(out) :: E_rad,S_rad ! Campos E_rad(r,t) e S_rad(r,t) a serem calculados

! Local:
real(dp) :: tr       ! Tempo retardado
real(dp),parameter :: mu_0 = 4*3.14159265_dp ! Nas unidades normalizadas, 4*pi*e0=1, mu0*e0=1/c^2 = 1, mu0 = 1/e0 = 4*pi
real(dp), dimension(3) :: w_tr, v_tr, a_tr ! Vetores de pos, vel, acel de q no tempo tr.
real(dp), dimension(3) :: delta_r, u       ! delta_r e u, conforme documentado no enunciado
! real(dp), dimension(3) :: B_rad            ! Campo B_rad(r,t)

! Inicialização:
! Primeiro, calculemos o tempo retardado associado ao par (r,t):
tr = find_root(r,t)
! Em seguida, calculemos os vetores w(tr), v(tr), a(tr) no tempo retardado:
w_tr = pos_q(tr)
v_tr = vel_q(tr)
a_tr = acel_q(tr)
! Finalmente, resta definir u e delta_r
delta_r = r - w_tr
u = normalize(delta_r) - v_tr

! Execução:
! Obviamente, faltava completar a expressão de E. Ademais, a velocidade aqui estava errada.
E_rad = norm(delta_r)/(dot_product(delta_r,u))**3 * cross(delta_r,cross(u,a_tr)) ! Apenas a parcela de radiação será
                                                                                 ! calculada.
S_rad = (1/mu_0)*norm(E_rad)**2*normalize(delta_r)

end subroutine 
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
subroutine average_fields_v3(r, E_rad, S_rad)
!! Subrotina que fornece |E_rad|^2, |S_rad| médios em uma posição r.
!! Esta versão tenta ser mais eficiente que a anterior, modificando métodos de integração.
! Esta é a subrotina correta a se usar. As demais não são fornecem os resultados esperados.


! IN
real(dp), dimension(3), intent(in) :: r
! OUT
real(dp), intent(out) :: E_rad,S_rad
! Local
integer, parameter :: N=100 ! Este é o mesmo parâmetro N a ser usado na função de integração por arrays
real(dp), dimension(0:N) :: E_norms=0._dp,S_norms=0._dp ! Cada linha de E_temp será um |E(r,t)|^2 diferente;
                                                      ! cada linha de S_temp será um |S(r,t)| diferente (sem ^2)
real(dp), dimension(3) :: E_temp, S_temp ! Vetores temporários para chamar a subrotina
integer  :: i       ! Índice de loop.
real(dp) :: t_i,dt  ! Variáveis t_i = t_a + dt*i, i=0,1,...,N; dt = (t_b-t_a)/N

! Inicialização:
dt = (periodo-0._dp)/N
! Execução:
! Para cada instante de tempo, calcularei as linhas de E_temp e S_temp:
do i=0,N
    t_i = 0._dp + i*dt
    call calc_fields(r,t_i,E_temp,S_temp)
    E_norms(i) = dot_product(E_temp,E_temp) ! |E(r,t)|^2
    S_norms(i) = norm(S_temp)               ! |S(r,t)|
end do
! Finalmente, resta integrar cada um dos arrays.
do i=1,3
    E_rad = simpson_array_int(E_norms,dt)/periodo
    S_rad = simpson_array_int(S_norms,dt)/periodo
end do

end subroutine
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
subroutine average_fields_v2(r, E_rad, S_rad)
!! Subrotina que fornece E_rad, S_rad médios em uma posição r.
!! Esta versão tenta ser mais eficiente que a anterior, modificando métodos de integração.
! average_fields_v2 utiliza a função simpson_array_int para otimizar a execução.


! IN
real(dp), dimension(3), intent(in) :: r
! OUT
real(dp), dimension(3), intent(out) :: E_rad,S_rad
! Local
integer, parameter :: N=100 ! Este é o mesmo parâmetro N a ser usado na função de integração por arrays
real(dp), dimension(0:N,3) :: E_matrix=0._dp,S_matrix=0._dp ! Cada linha de E_matrix será um E(r,t) diferente; o mesmo para S_temp
real(dp), dimension(3) :: E_temp, S_temp ! Vetores temporários para chamar a subrotina
integer  :: i       ! Índice de loop.
real(dp) :: t_i,dt  ! Variáveis t_i = t_a + dt*i, i=0,1,...,N; dt = (t_b-t_a)/N

! Inicialização:
dt = (periodo-0._dp)/N
! Execução:
! Para cada instante de tempo, calcularei as linhas de E_temp e S_temp:
do i=0,N
    t_i = 0._dp + i*dt
    call calc_fields(r,t_i,E_temp,S_temp)
    E_matrix(i,:) = E_temp
    S_matrix(i,:) = S_temp
end do
! Finalmente, resta executar outro loop integrando cada uma das três colunas de E_temp e S_temp, que correspondem
! às direções x,y,z:
do i=1,3
    E_rad(i) = simpson_array_int(E_matrix(:,i),dt)/periodo
    S_rad(i) = simpson_array_int(S_matrix(:,i),dt)/periodo
end do

end subroutine
!**********************************************************************************************************************************!


!**********************************************************************************************************************************!
subroutine average_fields_v1(r, E_rad, S_rad)
!! Subrotina que fornece E_rad, S_rad médios em uma posição r.
! Para manter a estrutura da integral de Simpson, decidi criar funções internas que fornecem E(r,t) e S(r,t)
! com r fixo e t como entrada.
! Este método é altamente ineficiente, porque recalcula E(r,t) e B(r,t) duas vezes cada.

!IN
real(dp), allocatable, dimension(:), intent(in) :: r

!OUT
real(dp), allocatable, dimension(:), intent(out) :: E_rad,S_rad

! Execução:
E_rad = simpson_int(func_E_t,0._dp,periodo)/periodo ! <E> = int(E dt)/Período
S_rad = simpson_int(func_S_t,0._dp,periodo)/periodo ! <S> = int(E dt)/Período

! Funções internas:
contains
    function func_E_t(t) result(E_t)
    !! Função interna que fornece E(r,t)
        ! IN
        real(dp), intent(in) :: t
        ! OUT
        real(dp),allocatable,dimension(:) :: E_t
        ! Local
        real(dp),allocatable,dimension(:) :: S_t

        call calc_fields(r,t,E_t,S_t) ! Lembre-se: a variável "r" está definida no escopo superior.
        ! Após a execução da subrotina, E_t e S_t estarão devidamente calculados.
    end function

    function func_S_t(t) result(S_t)
    !! Função interna que fornece S(r,t)

        ! IN
        real(dp), intent(in) :: t
        ! OUT
        real(dp),allocatable,dimension(:) :: E_t
        ! Local
        real(dp),allocatable,dimension(:) :: S_t

        call calc_fields(r,t,E_t,S_t) ! Lembre-se: a variável "r" está definida no escopo superior.
    end function

end subroutine
!**********************************************************************************************************************************!


!**********************************************************************************************************************************!
function pos_q(t) result(w)
!! Função da posição w(t) = (wx,wy,wz) da carga q.

! IN
real(dp), intent(in) :: t ! Instante de tempo.

! OUT
real(dp),dimension(3) :: w ! w(t), vetor

! Execução
w = (/A*cos(omega*t), 0._dp, 0._dp/)
end function
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
function vel_q(t) result(v)
!! Função da velocidade v(t) = (vx,vy,vz) da carga q.

! IN
real(dp), intent(in) :: t ! Instante de tempo.

! OUT
real(dp),dimension(3) :: v ! v(t), vetor

! Execução
v = (/ -A*omega*sin(omega*t), 0._dp, 0._dp/)

end function
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
function acel_q(t) result(acel)
!! Função da aceleração a(t) = (ax,ay,az) da carga q.

! IN
real(dp), intent(in) :: t ! Instante de tempo.

! OUT
real(dp),dimension(3) :: acel ! v(t), vetor

! Execução
acel = (/ -A*omega**2*cos(omega*t), 0._dp, 0._dp/)

end function
!**********************************************************************************************************************************!

!**********************************************************************************************************************************!
subroutine write_data()
!! Esta subrotina escreve os dados <|E^2|> e <|S|> como função de r.

! Local
integer :: E_unit=10, S_unit=11
character(len=*), parameter :: E_path="./out/E_avg.dat", S_path="./out/S_avg.dat"
real(dp), dimension(3), parameter :: direction = [1._dp,1._dp,1._dp]/sqrt(3._dp) ! Indica a direção de r
real(dp), dimension(3) :: r             ! r, como especificado em average_fields_v3
real(dp)               :: E_rad, S_rad  ! E_rad e S_rad, como especificado em average_fields_v3
real(dp) :: r_norm ! Norma de r, norm_r = norm_0 + dr*i
real(dp) :: dr, norm_0
real(dp),parameter :: omega_0=2._dp, d_omega=1._dp ! omega = omega_0 + d_omega*i
integer :: i ! Loop
integer,parameter :: N=25 ! Quantidade de pontos = N+1

! Inicialização:
! Primeiro, devo me certificar de usar omega e A corretos:
call read_const(omega_0)

! Como o enunciado sugere |r| >> c/omega = 1/omega, bons valores para dr e norm_0 são:
dr = 100._dp/omega
norm_0 = 1000._dp/omega
! Abre os documentos de texto
open(unit=E_unit,action="write",status="replace",file=E_path)
open(unit=S_unit,action="write",status="replace",file=S_path)
! Cria a linha inicial para orientação:
write(E_unit,*) "|r|          <|E(r,t)|^2>"
write(S_unit,*) "omega          <|S(r,t)|>"


! Execução:
! Loop de iteração, inicialmente, sobre E(r,t) contra o omega inicial do programa.
do i=0,N
    ! print*, "Loop 1, i = ",i, "omega = ", omega
    r_norm = norm_0 + dr*i
    r = direction*r_norm
    call average_fields_v3(r,E_rad,S_rad)
    write(E_unit,*) r_norm, E_rad
end do
! Em seguida, loop de iteração sobre differnetes valores de omega para um mesmo r fixo
! O loop anterior já utiliza um r grande, mas é bom reiniciar para um r que seja >> que qualquer um dos
! omegas a seguir.
r = direction*1.e4_dp
do i=0,N
    call read_const(omega_0 + d_omega*i) ! Ajusta o novo omega
    ! print*, "Loop 2, i = ", i, "omega = ", omega, "A = ", A
    call average_fields_v3(r,E_rad,S_rad)
    write(S_unit,*) omega, S_rad
end do

end subroutine
!**********************************************************************************************************************************!

end module
