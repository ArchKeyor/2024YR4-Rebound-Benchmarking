# O Algoritmo Histórico de Everhart (1985) e o Integrador RADAU15

> **Referência Bibliográfica Original:**  
> Everhart, E. (1985). *An efficient integrator that uses Gauss-Radau spacings.*  
> In: Carusi, A., Valsecchi, G. B. (eds) *Dynamics of Comets: Their Origin and Evolution.*  
> Proceedings of IAU Colloquium No. 83, Rome, Italy, June 11–15, 1984.  
> Astrophysics and Space Science Library, vol 115. Cambridge University Press / D. Reidel, Dordrecht, pp. 185–202.  
> **DOI:** [10.1017/S0252921100083913](https://doi.org/10.1017/S0252921100083913)  
> **Artigo Prévio Fundamental:** Everhart, E. (1974). *Implicit single-sequence methods for integrating orbits.* Celestial Mechanics, 10(1), 35–55.

---

## 1. Contexto Histórico e Motivação

Durante os anos 1960 e 1970, o astrônomo e físico **Edgar Everhart** (Universidade de Denver) dedicava-se ao estudo da evolução dinâmica de cometas e da captura de pequenos corpos pelos planetas gigantes do Sistema Solar (Júpiter e Saturno).

### O Dilema dos Integradores Numéricos Convencionais
As equações diferenciais da mecânica celeste para $N$ corpos gravitacionais:
$$\mathbf{\ddot{r}}_i = - G \sum_{j \neq i}^{N} m_j \frac{\mathbf{r}_i - \mathbf{r}_j}{\|\mathbf{r}_i - \mathbf{r}_j\|^3}$$
apresentam singularidades quando $\|\mathbf{r}_i - \mathbf{r}_j\| \to 0$ (encontros próximos / *close encounters*).

Na época, os métodos predominantes eram:
1. **Runge-Kutta de 4ª ordem (RK4):** Baixa ordem; requeria passos minúsculos, acumulando erros de truncamento gigantescos ao longo de milhares de órbitas.
2. **Métodos Multipasso (Adams-Bashforth-Moulton):** Exigiam histórico de passos prévios uniformes; sofriam para alterar o passo dinamicamente durante aproximações rápidas de periélio ou encontros planetários sem perda severa de estabilidade.
3. **Bulirsch-Stoer (Extrapolação de Richardson):** Excelente precisão, porém computacionalmente caro devido ao grande número de avaliações de aceleração por passo.

Everhart propôs um método **de passo único (*single-step*), implícito, de altíssima ordem (15ª ordem)**, que fosse auto-iniciante (*self-starting*) e capaz de adaptar o tamanho do passo $h$ de forma suave e eficiente através do monitoramento de ordens superiores de diferenças divididas. Esse método ficou mundialmente conhecido como **RADAU15** ou **RA15**.

---

## 2. Fundamentos Matemáticos do RADAU15

### 2.1 A Quadratura de Gauss-Radau
Na integração numérica clássica de Gauss-Legendre, nenhum dos extremos do intervalo $[0, h]$ é fixado como nó de amostragem.  
Na **quadratura de Gauss-Radau**, um dos pontos extremos (neste caso, o ponto inicial $t_0$, ou seja, $h_0 = 0$) é fixado obrigatoriamente como nó de quadratura. Os demais nós $h_1, h_2, \dots, h_k$ são escolhidos como as raízes de polinômios ortogonais específicos (derivados dos polinômios de Jacobi/Legendre).

Para um sistema de ordem 15 ($k=7$ nós intermediários além do ponto inicial):
$$0 = h_0 < h_1 < h_2 < h_3 < h_4 < h_5 < h_6 < h_7 = 1$$
*(normalizados para o intervalo de passo unitário $s \in [0, 1]$)*.

Os nós exatos de Gauss-Radau de 15ª ordem calculados por Everhart são dados pelas raízes:
- $h_1 \approx 0.0562625605369$
- $h_2 \approx 0.1802406917369$
- $h_3 \approx 0.3526249659989$
- $h_4 \approx 0.5471536138135$
- $h_5 \approx 0.7342124755499$
- $h_6 \approx 0.8853209468391$
- $h_7 = 1.0000000000000$

### 2.2 Polinômio de Interpolação de Aceleração
Everhart aproxima a aceleração vetorial $\mathbf{a}(t)$ no intervalo $[t_n, t_n + h]$ por um polinômio interpolador de grau 7 expresso em termos de coeficientes vetoriais $\mathbf{b}_1, \dots, \mathbf{b}_7$:
$$\mathbf{a}(t_n + \tau h) = \mathbf{a}_0 + \sum_{j=1}^{7} \mathbf{b}_j \, \prod_{m=0}^{j-1} (\tau - h_m)$$
onde $\tau = (t - t_n)/h \in [0, 1]$.

Integrando analiticamente duas vezes o polinômio em relação a $\tau$, obtêm-se expressões fechadas para a velocidade e posição em qualquer ponto intermediário ou no final do passo:
$$\mathbf{v}(t_n + h) = \mathbf{v}_n + h \left[ \mathbf{a}_0 + \sum_{j=1}^{7} \beta_j \mathbf{b}_j \right]$$
$$\mathbf{r}(t_n + h) = \mathbf{r}_n + h \mathbf{v}_n + h^2 \left[ \frac{1}{2}\mathbf{a}_0 + \sum_{j=1}^{7} \gamma_j \mathbf{b}_j \right]$$
onde os escalares $\beta_j$ e $\gamma_j$ são constantes pré-calculadas analiticamente a partir da integração dos polinômios de Radau.

### 2.3 O Ciclo Preditor-Corretor (Loop de Gauss-Seidel)
Como a aceleração depende da posição e da velocidade que estão sendo integradas ($\mathbf{a} = \mathbf{a}(t, \mathbf{r}, \mathbf{v})$), o sistema é resolvido de forma **implícita**.
1. **Previsão (*Predictor*):** Extrapolam-se os coeficientes $\mathbf{b}_j$ do passo anterior para o novo passo.
2. **Subpassos:** Avaliam-se as forças nos 7 nós $h_1, \dots, h_7$.
3. **Correção (*Corrector*):** Recalculam-se as diferenças divididas e ajustam-se os $\mathbf{b}_j$.
4. **Convergência:** Repetem-se 2 a 4 iterações até que a variação de $\mathbf{b}_7$ seja inferior à tolerância estipulada.

### 2.4 Controle Adaptativo de Passo
O termo de mais alta ordem, $\mathbf{b}_7$, representa a derivada de 15ª ordem da posição. Everhart utilizou esse termo para prever o erro de truncamento local $E_{local} \propto \|\mathbf{b}_7\| h^8$.  
O próximo passo $h_{novo}$ é dimensionado por:
$$h_{novo} = h_{atual} \left( \frac{\epsilon}{\|\mathbf{b}_7\|} \right)^{1/7}$$
onde $\epsilon$ é o parâmetro de tolerância de precisão fornecido pelo usuário (ex: $10^{-12}$).

---

## 3. O Legado EVORB (Fortran 77)

O simulador herdado neste projeto (`data/evorb_raw/evorb15corr.for`) é uma implementação direta do algoritmo RADAU15 de Everhart (1985), adaptado com:
- Blocos globais `COMMON` para compartilhamento de vetores de estado entre sub-rotinas (`EVORB`, `FORCE`, `RADAU`).
- Suporte a forças não-gravitacionais comutadas via *flags* inteiras.
- Formato de entrada de dados de largura fixa (`datos.dat`).

### Limitações do EVORB15 em Fortran 77
1. **Acúmulo de Erro de Ponto Flutuante (*Roundoff Drift*):** Em integrações de mais de $10^5$ passos, a soma sucessiva $x_{n+1} = x_n + \Delta x$ em precisão dupla padrão (IEEE 754 de 64 bits) sofre de cancelamento catastrófico. O erro de energia acumula-se proporcionalmente a $\sqrt{N}$ (passeio aleatório) ou $N$ (viés sistemático).
2. **Arquitetura Não-Reentrante:** O uso de `COMMON /VARS/` impossibilita paralelização multi-thread (OpenMP/MPI) ou execução de múltiplos sistemas simultâneos no mesmo processo.
3. **Falta de Interface Moderna:** A amarração a arquivos de texto posicionais com `FORMAT(3D22.15)` impedia automação e integração com pipelines modernos de dados astronômicos (JSON/AstroPy/Parquet).

---

## 4. Como o IAS15 (Rein & Spiegel, 2015) Aperfeiçoou Everhart (1985)

O integrador **IAS15** (*Integrator with Adaptive Step-size 15th-order*), criado por Hanno Rein e Daniel S. Spiegel e incorporado nativamente ao **REBOUND**, é o herdeiro direto do trabalho de Everhart.

As quatro grandes inovações do IAS15 sobre o Everhart (1985) foram:

| Característica | Everhart (1985) / EVORB15 | IAS15 (Rein & Spiegel, 2015) |
| :--- | :--- | :--- |
| **Precisão de Máquina** | Limitado por erro de arredondamento em passos longos | Mantém precisão de máquina ($10^{-16}$) mesmo após $10^9$ órbitas |
| **Tratamento de Ponto Flutuante** | Soma ingênua ($\mathbf{r} \leftarrow \mathbf{r} + \Delta \mathbf{r}$) | **Soma Compensada de Babuška-Neumaier** (elimina perda de dígitos significativos) |
| **Controle de Passo Adaptativo** | Baseado apenas em $\|\mathbf{b}_7\|$ | Duplo controle: avalia tanto $\mathbf{b}_7$ quanto a conservação infinitesimal de energia |
| **Estrutura de Software** | Fortran 77 monossessão com variáveis estáticas | C puro ISO C99 thread-safe, com structs de simulação independentes e bindings Python de alto desempenho |
| **Forças Não-Conservativas** | Requer modificação manual de sub-rotinas | Suporte unificado via *function pointers* sem overhead de compilação |

---

## 5. Referências Cruzadas no Repositório

- **Código Legado Fortran:** [evorb15corr.for](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/data/evorb_raw/evorb15corr.for)
- **Implementação Moderna (REBOUND/IAS15):** [rebound_runner.py](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/src/rebound_runner.py)
- **Artigo IAS15 em PDF:** [Rein_Spiegel_2015_IAS15.pdf](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/references/Rein_Spiegel_2015_IAS15.pdf)
- **Script de Validação:** [verify_confronto.py](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/scripts/verify_confronto.py)
