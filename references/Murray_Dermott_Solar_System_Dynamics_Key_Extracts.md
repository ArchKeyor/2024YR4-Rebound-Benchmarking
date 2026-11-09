# Compêndio Teórico: Murray & Dermott (1999) - Solar System Dynamics

> **Referência Bibliográfica:**  
> Murray, C. D., & Dermott, S. F. (1999). *Solar System Dynamics.*  
> Cambridge University Press, Cambridge, UK, 592 pp.  
> **ISBN:** 978-0-521-57597-3 (Paperback: 978-0-521-57295-8)  
> **DOI:** [10.1017/CBO9781139174817](https://doi.org/10.1017/CBO9781139174817)  
> **Website Oficial do Livro:** [ssdbook.maths.qmul.ac.uk](http://ssdbook.maths.qmul.ac.uk)

---

## 1. Visão Geral da Obra

Publicado pela Cambridge University Press, o livro **"Solar System Dynamics"** de Carl D. Murray (Queen Mary University of London) e Stanley F. Dermott (University of Florida) é universalmente aclamado como o livro-texto de referência para a dinâmica orbital e mecânica celeste moderna em cursos de pós-graduação no mundo todo.

O livro unifica o formalismo analítico clássico da mecânica celeste (Lagrange, Laplace, Gauss, Poincaré) com métodos computacionais modernos de simulação numérica de $N$ corpos.

---

## 2. Conceitos Centrais Aplicados ao OrbitREBOUND

### 2.1 O Problema dos Dois Corpos (Capítulo 2)

O movimento relativo de um corpo de massa $m_2$ em torno de $m_1$ sob a atração gravitacional newtoniana obedece:
$$\mathbf{\ddot{r}} + \mu \frac{\mathbf{r}}{r^3} = 0$$
onde $\mu = G(m_1 + m_2)$.

#### A Equação da Força Viva (*Vis-Viva*)
A velocidade escalar $v$ em qualquer distância radial $r$ em uma órbita elíptica de semieixo maior $a$ é dada por:
$$v^2 = \mu \left( \frac{2}{r} - \frac{1}{a} \right)$$
- Para órbita elíptica: $a > 0$ e $E < 0$.
- Para trajetória parabólica (fuga): $a \to \infty$ e $E = 0$.
- Para trajetória hiperbólica (encontro aberto): $a < 0$ e $E > 0$, com velocidade no infinito $v_\infty^2 = -\mu/a$.

#### Elementos Orbitais Osculantes e Vetores de Estado
Um estado orbital no espaço tridimensional é univocamente descrito ou por 6 variáveis cartesianas $(\mathbf{r}, \mathbf{v}) = (x, y, z, v_x, v_y, v_z)$ ou por 6 elementos orbitais keplerianos:
1. $a$ — **Semieixo Maior (*Semi-major axis*):** Define o tamanho da órbita e o período orbital via 3ª Lei de Kepler ($P^2 = 4\pi^2 a^3 / \mu$).
2. $e$ — **Excentricidade (*Eccentricity*):** Define a forma da cônica ($e=0$: círculo; $0 < e < 1$: elipse; $e=1$: parábola; $e>1$: hipérbole).
3. $i$ — **Inclinação (*Inclination*):** Ângulo entre o plano orbital e o plano fundamental (eclíptica J2000).
4. $\Omega$ — **Longitude do Nó Ascendente (*Longitude of Ascending Node*):** Posição onde o corpo cruza a eclíptica de sul para norte.
5. $\omega$ — **Argumento do Periélio (*Argument of Pericenter*):** Orientação da elipse em seu próprio plano, medido a partir do nó ascendente.
6. $M$ — **Anomalia Média (*Mean Anomaly*):** Parâmetro angular temporal que localiza o corpo na órbita, relacionado à Anomalia Excêntrica $E$ pela **Equação de Kepler**:
   $$M = E - e \sin E$$

---

### 2.2 O Problema Restrito dos Três Corpos (CR3BP) (Capítulo 3)

Quando um asteroide ou cometa de massa infinitesimal ($m_3 \to 0$) orbita sob a influência mútua de dois corpos massivos (ex: Sol com massa $m_1$ e Terra com massa $m_2$), o sistema é modelado pelo problema circular restrito.

#### O Sistema de Coordenadas Rotativo e a Integral de Jacobi
No referencial que gira solidário com os dois corpos primários com velocidade angular $n$, existe uma integral de movimento exata, a **Integral de Jacobi ($C_J$)**:
$$C_J = 2\Omega(x,y) - (u^2 + v^2 + w^2)$$
onde $u, v, w$ são as componentes de velocidade no referencial girante e o potencial efetivo corrigido pela força centrífuga é:
$$\Omega(x,y) = \frac{1}{2} n^2 (x^2 + y^2) + \frac{G m_1}{r_1} + \frac{G m_2}{r_2}$$

As superfícies de velocidade zero ($u^2+v^2+w^2 = 0$) definem curvas de nível limitantes que impedem que um asteroide cruze determinadas regiões do espaço se seu $C_J$ for suficientemente alto.

#### A Esfera de Hill (*Hill Sphere*)
A região de dominância gravitacional de um planeta de massa $m_p$ orbitando o Sol ($M_\odot$) a uma distância $a_p$ com excentricidade $e_p$ é delimitada pelo raio da Esfera de Hill:
$$r_H \approx a_p (1 - e_p) \left( \frac{m_p}{3 M_\odot} \right)^{1/3}$$

*Aplicação prática ao asteroide 2024 YR4:*
- Massa da Terra: $m_\oplus \approx 5.972 \times 10^{24}\text{ kg}$
- Massa do Sol: $M_\odot \approx 1.989 \times 10^{30}\text{ kg}$
- Raio de Hill da Terra:
  $$r_{H, \oplus} \approx 1\text{ AU} \times \left( \frac{5.972 \times 10^{24}}{3 \times 1.989 \times 10^{30}} \right)^{1/3} \approx 0.010\text{ AU} \approx 1.496.000\text{ km} \approx 3.89\text{ LD}$$
- **Aproximação de 2024 YR4:** A distância mínima nominal projetada para 22/12/2032 é de **$0.00186\text{ AU} \approx 278.330\text{ km} \approx 0.72\text{ LD}$**.  
  Isso significa que o asteroide **penetra profundamente no interior da Esfera de Hill da Terra** ($0.18\, r_H$), onde a gravidade terrestre supera a atração de maré solar, tornando a perturbação geocêntrica o efeito dominante da trajetória.

---

### 2.3 Teoria de Perturbações e Equações de Gauss (Capítulos 6 e 7)

Na presença de acelerações perturbadoras $\mathbf{F}_{pert} = R \mathbf{\hat{r}} + T \mathbf{\hat{\theta}} + N \mathbf{\hat{n}}$ (componente radial $R$, transversal $T$ e normal $N$), a variação instantânea dos elementos osculantes é regida pelas **Equações Planetárias de Gauss**:

$$\frac{da}{dt} = \frac{2 a^2 \sqrt{1-e^2}}{\mu} \left[ R e \sin f + T \frac{p}{r} \right]$$
$$\frac{de}{dt} = \frac{\sqrt{1-e^2}}{\mu} \left[ R \sin f + T \left( \cos f + \frac{e + \cos f}{1 + e \cos f} \right) \right]$$
$$\frac{di}{dt} = \frac{r \cos(\omega + f)}{\mu \sqrt{1-e^2}} N$$

Essas equações são cruciais para modelar forças não-gravitacionais (como o **Efeito Yarkovsky**, onde a emissão térmica anisotrópica gera uma força na componente transversal $T$, causando deriva secular no semieixo maior $a$).

---

### 2.4 Encontros Próximos e o Invariante de Tisserand (Capítulo 8)

Durante um encontro próximo rápido com um planeta, a integral de Jacobi no referencial heliocêntrico se traduz no **Parâmetro de Tisserand ($T_P$)**:
$$T_P = \frac{a_p}{a} + 2 \sqrt{\frac{a}{a_p} (1 - e^2)} \cos i$$

- O Parâmetro de Tisserand permanece aproximadamente constante antes e depois de um *close encounter*.
- Ele conecta o semieixo maior ($a$), a excentricidade ($e$) e a inclinação ($i$): se o encontro planetário reduzir o semieixo maior ($a$), ele obrigatoriamente forçará uma compensação na excentricidade ($e$) e/ou na inclinação ($i$).
- Este critério classifica populações de pequenos corpos:
  - Asteroides do Cinturão Principal: $T_J > 3.1$
  - Família de Cometas de Júpiter (JFCs): $2 < T_J < 3$
  - Cometasquase-parabólicos (Nuvem de Oort): $T_J < 2$

---

## 3. Por Que a Simulação Numérica de Alta Ordem é Essencial

Murray & Dermott demonstram formalmente que as soluções analíticas e perturbações seculares (desenvolvimentos em séries de potências de $e$ e $i$) **divergem e colapsam** quando ocorrem:
1. Encontros dentro da esfera de Hill ($r < r_H$);
2. Ressonâncias de movimento médio de ordem alta com cruzamento de órbitas;
3. Caos hamiltoniano por sobreposição de ressonâncias (critério de Chirikov).

Nesses regimes, a integração numérica direta através de algoritmos de alta precisão como o **IAS15** (desenvolvido por Rein & Spiegel com base em Everhart) torna-se o único instrumento científico capaz de propagar a probabilidade de impacto e a evolução a longo prazo do corpo.

---

## 4. Integração no Repositório OrbitREBOUND

- **Módulo de Leitura e Propagação:** [rebound_runner.py](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech\Projetos\2024YR4-Rebound-Benchmarking/src/rebound_runner.py)
- **Leitor de Vetores de Estado e Condições Iniciais:** [datos_parser.py](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech\Projetos\2024YR4-Rebound-Benchmarking/src/datos_parser.py)
- **Relatório Oficial do Asteroide 2024 YR4:** [NASA_JPL_2024YR4_CNEOS_Report.md](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech\Projetos\2024YR4-Rebound-Benchmarking/references/NASA_JPL_2024YR4_CNEOS_Report.md)
