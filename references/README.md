# Biblioteca Científica e Referências Teóricas — OrbitREBOUND

Este diretório contém a base bibliográfica, artigos originais em PDF, compêndios matemáticos e relatórios astronômicos oficiais da NASA JPL / CNEOS que fundamentam o desenvolvimento e os benchmarks do projeto **OrbitREBOUND**.

---

## 1. Índice do Acervo

| Documento / Artigo | Autores / Ano | Formato / Arquivo | Descrição & Relevância |
| :--- | :--- | :---: | :--- |
| **O Artigo do IAS15** | Rein, H. & Spiegel, D. S. (2015) | [PDF Local](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/references/Rein_Spiegel_2015_IAS15.pdf) (547 KB) | *MNRAS 446, 1424–1437*. Define a formulação do integrador IAS15 de 15ª ordem, soma compensada e precisão de máquina. |
| **O Artigo Original do REBOUND** | Rein, H. & Liu, S.-F. (2012) | [PDF Local](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/references/Rein_Liu_2012_REBOUND.pdf) (291 KB) | *A&A 537, A128*. Apresenta a arquitetura de simulação $N$-corpos em C com alta eficiência computacional. |
| **O Algoritmo Histórico de Everhart (1985)** | Everhart, E. (1985) | [Dossiê Matemático](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/references/Everhart_1985_RADAU15_Algorithmic_Reference.md) | *IAU Colloq. 83, pp. 185–202*. A gênese do algoritmo RADAU15 e da rotina `evorb15corr.for` original em Fortran 77. |
| **A "Bíblia" da Dinâmica Orbital** | Murray, C. D. & Dermott, S. F. (1999) | [Compêndio Teórico](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/references/Murray_Dermott_Solar_System_Dynamics_Key_Extracts.md) | *Solar System Dynamics, Cambridge Univ. Press*. Formulação analítica: problema dos 2 e 3 corpos, esfera de Hill, equações de Gauss e parâmetro de Tisserand. |
| **Relatório Oficial Asteroide 2024 YR4** | NASA JPL / CNEOS (2025–2026) | [Relatório Técnico](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/references/NASA_JPL_2024YR4_CNEOS_Report.md) | Parâmetros astrométricos, tabela de aproximações próximas (CNEOS CAD) e histórico no sistema de risco de impacto Sentry para o encontro de 2032. |
| **Dados Brutos de Telemetria JPL** | NASA JPL Horizons / SBDB | [Pasta `data/`](file:///C:/Users/Kauê/Documents/02_Programacao_e_Tech/Projetos/2024YR4-Rebound-Benchmarking/references/data/) | Payloads JSON originais da API JPL Small-Body Database e Close Approach Data. |

---

## 2. Guia de Estudos Fundamentais

### 2.1 Por que o REBOUND é tão rápido?
1. **Núcleo em C Puro (ISO C99):** Toda a física, vetores de estado, cálculo de forças gravitacionais e passos de integração são executados em código nativo compilado, sem overhead do interpretador Python.
2. **Layout de Memória Contíguo e Cache-Friendly:** As partículas são armazenadas em estruturas contíguas na memória (`struct reb_particle*`), facilitando a vetorização por instruções SIMD (AVX2/AVX-512) da CPU e minimizando *cache misses*.
3. **Ausência de Alocações Dinâmicas no Loop Quente:** Não ocorrem alocações de memória (`malloc`) nem *garbage collection* durante o avanço temporal dos passos.
4. **Algoritmos Especializados e Paralelização Opcional:** Suporte nativo a OpenMP e MPI, além de acelerações em árvore (*octree / Barnes-Hut*) quando o número de partículas é grande.

### 2.2 O que é o IAS15 e por que ele supera o EVORB15?
- **Gênese:** O IAS15 (*Integrator with Adaptive Step-size 15th-order*) é o refinamento moderno do método de quadratura de Gauss-Radau desenvolvido por Edgar Everhart em 1985 (a base do EVORB15).
- **Ordem 15 com 7 Subpassos:** Em cada passo temporal, avalia a aceleração em 7 nós posicionados estrategicamente nas raízes dos polinômios de Radau, alcançando precisão de 15ª ordem.
- **Soma Compensada de Babuška-Neumaier:** No EVORB15 (Fortran 77), a adição contínua de pequenos incrementos $\Delta \mathbf{r}$ a posições grandes $\mathbf{r}$ causa perda gradativa dos dígitos menos significativos em precisão dupla (64 bits). O IAS15 rastreia os bits descartados em registradores de compensação, garantindo **conservação de energia e momento angular até a precisão da máquina ($10^{-16}$) mesmo após $10^9$ passos**.
- **Segurança de Memória e Modernidade:** O EVORB15 utiliza blocos estáticos `COMMON /VARS/` do Fortran 77, que impossibilitam concorrência *multi-thread* e tornam a compilação vulnerável. O IAS15 é reentrante, encapsulado e integrado a ferramentas modernas de ciência de dados.

### 2.3 A Linha do Tempo Histórica da Gravitação e Algoritmos $N$-Corpos

```
[1687] Isaac Newton: Principia Mathematica (Lei da Gravitação Universal)
         │
         ▼
[1760] Leonhard Euler: Método de Euler (1ª ordem numérica)
         │
         ▼
[1772] Joseph-Louis Lagrange: Formulação Hamiltoniana e os 5 Pontos de Equilíbrio (L1 a L5)
         │
         ▼
[1895] Carl Runge e Wilhelm Kutta: Integradores Runge-Kutta (RK4)
         │
         ▼
[1974] Edgar Everhart: Artigo precursor sobre passos implícitos de Gauss-Radau
         │
         ▼
[1985] Edgar Everhart: Algoritmo RADAU15 (origem da rotina EVORB em Fortran 77)
         │
         ▼
[1991] Jack Wisdom & Matthew Holman: Integradores Simpléticos para dinâmica solar de longo prazo
         │
         ▼
[1999] Carl Murray & Stanley Dermott: Publicação de "Solar System Dynamics"
         │
         ▼
[2012] Hanno Rein & Shang-Fei Liu: Publicação do framework REBOUND
         │
         ▼
[2015] Hanno Rein & Daniel Spiegel: Criação do IAS15 (Precisão de máquina com Babuška-Neumaier)
         │
         ▼
[2024-2026] Descoberta do asteroide 2024 YR4 e Benchmarking no projeto OrbitREBOUND
```

---

## 3. Como Atualizar os Dados e Relatórios

Para atualizar os dados astrométricos do asteroide 2024 YR4 a qualquer momento a partir das APIs da NASA JPL:
```bash
python scripts/fetch_jpl_references.py
```
Esse script atualiza os arquivos JSON em `references/data/` e regenera automaticamente o relatório `references/NASA_JPL_2024YR4_CNEOS_Report.md`.
