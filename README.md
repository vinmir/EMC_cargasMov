# EMC_cargasMov
Campos eletromagnéticos provocados por cargas em movimento.

# Referências
+ Gould: Computer Simulation Methods
+ Grifiths: Electrodynamics
+ Nakanishi: Computational Physics

# Imagens
## Carga em movimento retilíneo uniforme
Uma carga se move ao longo do eixo $x$ com rapidez constante. Os vetores campo elétrico e campo magnético terão dependência temporal e espacial devido ao movimento da carga. Construindo uma malha cúbica (*lattice*) discretizada em $x,y,z$, pode-se determinar os campos. 
+ Campo elétrico na malha cúbica quando a carga estiver na origem (fixando-se $t=t_0\ne0$):
[a.pdf](https://github.com/vinmir/EMC_cargasMov/files/11450374/a.pdf)

Dada a simetria do sistema de coordenadas perante o movimento da partícula, a representação é a mesma nos planos $x,z$ ou $y,z$, confirmando o resultado analítico (Griffiths).

+ Módulo do campo magnético
[b.pdf](https://github.com/vinmir/EMC_cargasMov/files/11450428/b.pdf)

O campo magnético deve ser azimutal, pois a corrente efetiva está no sentido do vetor velocidade. Em função disso, basta calcular $B_z(0,y,0,t_0)$, como no gráfico acima. Nota-se um decaimento do campo com um afastamento da carga.

## Carga em movimento retilíneo acelerado
+ Campo elétrico: [a.pdf](https://github.com/vinmir/EMC_cargasMov/files/11450445/a.pdf)
  + Como esperado, há uma "cauda" devido ao atraso da informação da carga (Griffiths)  
+ Campo magnético: [b.pdf](https://github.com/vinmir/EMC_cargasMov/files/11450454/b.pdf)
  + Continua azimutal, decaimento mais rápido conforme o afastamento do eixo $x$

## Radiação de carga em movimento harmônico (MHS): campo elétrico
Na radiação assintótica de dipolo, $\langle E^2 \rangle$ tem dependência $1/r^2$ ao longo de um período de oscilação da carga. Este comportamento é constatado com clareza ao discretizar numericamente $\mathbf{E}(\mathbf{r},t)$ e $\mathbf{B}(\mathbf{r},t)$ na malha cúbica dos itens anteriores, resolvendo o tempo retardado devido ao atraso das informações da carga e computando vetorialmente cada um dos campos. Após regressão estatística ([a.pdf](https://github.com/vinmir/EMC_cargasMov/files/11450515/a.pdf)), constata-se o comportamento assintótico provido pela teoria.

## Radiação de carga em MHS: intensidade da radiação
A intensidade da radiação $|\mathbf{S}|$, para uma dada posição da malha, depende da frequência angular de oscilação da carga. De fato, a regressão estatística confirma o comportamento $S \propto \omega^4$, assim como o valor numérico do coeficiente multiplicativo ([b.pdf](https://github.com/vinmir/EMC_cargasMov/files/11450545/b.pdf)).
