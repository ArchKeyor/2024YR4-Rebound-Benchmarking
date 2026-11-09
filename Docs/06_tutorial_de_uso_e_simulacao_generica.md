# 🪐 Tutorial de Uso: Simulador Orbital Universal (OrbitREBOUND)

Este documento é o guia prático e tutorial completo para utilizar o **OrbitREBOUND** (`scripts/run_simulation.py`), um simulador orbital de N-corpos em Python/C desenvolvido como substituto moderno, aberto e de altíssima performance para o software legado **EVORB15** (Prof. Tabaré Gallardo, 2009).

O simulador é **100% genérico**: ele não está restrito ao asteroide 2024 YR4. Você pode utilizá-lo para simular **qualquer asteroide, cometa ou corpo menor do Sistema Solar**, seja para uma órbita nominal única ou para milhares de clones de incerteza orbital.

---

## 🎯 1. Como Funciona a Arquitetura

O simulador é composto por 3 módulos principais:

```
[ Arquivo datos.dat ] 
         │
         ▼
[ src/parsers/evorb_parser.py ] ──> Lê parâmetros, planetas e N partículas de teste
         │
         ▼
[ src/simulation/rebound_runner.py ] ──> Integrador IAS15 (Gauss-Radau 15ª ordem em C)
         │                               Calibração da constante de Gauss (k) e LLEVA
         ▼
[ scripts/run_simulation.py ] ──> Monitoramento de encontros, colisões e energia
         │
         ├──📂 final/    (Snapshots consolidados em T = T_final)
         └──📂 temporal/ (Séries temporais completas T = 0 até T_final)
```

---

## 🚀 2. Guia de Execução Rápida

### A) Modos de Execução: Individual ou Varredura Recursiva em Lote

O OrbitREBOUND possui inteligência de descoberta de arquivos (`rglob`): você pode apontar para uma pasta individual ou para uma pasta raiz com centenas de subpastas aninhadas.

#### 1. Execução Individual (Conjunto Único)
Para rodar uma pasta específica contendo um `datos.dat`:
```bash
python scripts/run_simulation.py -i data/evorb_raw/datos1
```
* **O que acontece:** O simulador processa os 408 corpos da pasta `datos1`, integrando pelo tempo definido e salvando em `data/rebound_output/datos1/`.

#### 2. Execução em Lote por Descoberta Recursiva Automática
Se você possui múltiplos experimentos ou pastas de dados (ex: `datos1`, `datos2`, ..., `datos21`, ou subníveis como `campanhas/2026/loteA/datos.dat`), basta passar o diretório pai:
```bash
python scripts/run_simulation.py -i data/evorb_raw
# ou para qualquer pasta 'datos':
python scripts/run_simulation.py -i caminho/para/datos
```
* **Como funciona o rastreamento recursivo:**
  1. O simulador vasculha a pasta de entrada e **entra recursivamente em todas as subpastas e níveis de profundidade** procurando por arquivos `datos.dat`.
  2. Apresenta um índice consolidado com todos os conjuntos encontrados na tela.
  3. Executa as simulações sequencialmente em lote, preservando a integridade de cada conjunto.
  4. Organiza automaticamente os resultados em `data/rebound_output/<nome_da_subpasta>/`, mantendo cada série de clones isolada e catalogada.

---

### B) Parâmetros e Opções da Linha de Comando (CLI)

| Flag | Tipo | Padrão | Descrição |
| :--- | :---: | :---: | :--- |
| `-i`, `--input-dir` | `str` | `data/evorb_raw/datos1` | Arquivo `datos.dat`, pasta individual ou pasta raiz para busca recursiva. |
| `-o`, `--output-dir` | `str` | `None` (automático) | Pasta onde as saídas serão salvas (padrão: `data/rebound_output/`). |
| `-y`, `--years` | `float` | `None` (lê do `.dat`) | Tempo de simulação em anos (ex: `-y 50.0` para 50 anos). |
| `--dt-temporal` | `float` | `None` (lê do `.dat`) | Passo de amostragem na pasta `temporal/`. Por padrão lê exatamente o `tsal` do `datos.dat` (tipicamente `0.001` anos = ~8,76 horas). |

#### Exemplos de comandos:

1. **Varredura Recursiva Automática de Toda a Pasta de Dados:**
   ```bash
   python scripts/run_simulation.py -i data/evorb_raw
   ```
   *(Localiza e roda todos os `datos.dat` presentes em quaisquer subpastas de `data/evorb_raw`).*

2. **Execução Padrão com Resolução Física Total (`tsal = 0.001` anos):**
   ```bash
   python scripts/run_simulation.py -i data/evorb_raw/datos1
   ```
   *(Lê o `tsal = 0.001` anos do arquivo original, executando os 20.001 passos de amostragem temporal com máxima fidelidade científica).*

3. **Simular por 50 anos com saída personalizada:**
   ```bash
   python scripts/run_simulation.py -i data/evorb_raw/datos2 -y 50.0 -o resultados/run_50anos
   ```

4. **Simular com amostragem customizada (ex: a cada 36 dias para arquivos menores):**
   ```bash
   python scripts/run_simulation.py -i data/evorb_raw/datos1 --dt-temporal 0.1
   ```

5. **Execução Isolada via Docker (com montagem de volume para saídas locais):**
   ```bash
   docker run --rm -v "${PWD}/data:/app/data" 2024yr4-rebound python scripts/run_simulation.py -i data/evorb_raw/datos1
   ```

---

## 🔬 3. Transparência Científica: Unidades, Conversões Físicas e Ordenação

Para evitar qualquer confusão ao comparar saídas do REBOUND com integradores legados em Fortran (EVORB15, Mercury, Swift), esta seção documenta com rigor a física e a convenção de dados adotada.

### A) Sistema de Unidades e Conservação de Energia ($GM$)
* **No Fortran EVORB:** O programa adota unidades astronômicas clássicas: massa em Massas Solares ($M_\odot$), distância em $\text{AU}$ e tempo em Anos Julianos ($365,25\text{ dias}$).
  Nesse sistema, a constante de gravitação do Sol é $GM = (k \cdot 365,25)^2 \approx 39,476926\text{ AU}^3/\text{ano}^2$ (onde $k = 0,01720209895$ é a constante de Gauss).
  A energia mecânica do sistema ($E = T + V$) é calculada em unidades de:
  $$[E] = M_\odot \cdot \frac{\text{AU}^2}{\text{ano}^2}$$
* **No REBOUND:** O motor opera internamente em unidades adimensionais onde $G = 1$ e o tempo é medido em radianos ($1/2\pi\text{ anos}$).
* **A Conversão do OrbitREBOUND:** Para garantir transparência e equivalência com o que a comunidade científica espera, o OrbitREBOUND converte automaticamente a energia interna do REBOUND multiplicando pelo fator $GM \approx 39,476926$.
  * Resultado: o `energia.csv` cospe diretamente $E \approx -0,0044324\ M_\odot\text{AU}^2/\text{ano}^2$, batendo com o `energia.sal` do Fortran ($-0,0044325$) na 5ª casa decimal.

---

### B) Estrutura das Pastas e Critérios de Ordenação dos Dados

A saída é dividida em duas pastas especializadas para atender a propósitos analíticos distintos:

#### 1. Pasta `temporal/` (Série Temporal Cronológica Contínua)
* **Objetivo:** Registro passo a passo da evolução dinâmica desde $T = 0$ até $T = T_{\text{max}}$.
* **Resolução Temporal:** Lê rigorosamente o `tsal` do `datos.dat` (`0,001\text{ anos} \approx 8,76\text{ horas}`). Essa alta densidade é indispensável para resolver encontros hiperbólicos na esfera de Hill sem efeito de *aliasing*.
* **Critérios de Ordenação:**
  * `encuent.csv`: Ordenado estritamente por **`[t_anos, clone_id]`** (ordem cronológica primária; em caso de passagens no mesmo instante, ordena pelo número do clone: 9, 10, 11...). Espelha com exatidão a ordem de emissão do `encuent.sal`.
  * `orbitas.csv` e `xyzv.csv`: Ordenados por **`[t_anos, N]`** (série temporal completa de 20.001 passos para cada corpo).
  * `energia.csv`: Ordenado por **`t_anos`** (20.001 amostras da conservação de energia com drift $\Delta E / E \sim 10^{-16}$).

#### 2. Pasta `final/` (Consolidação no Instante Final $T = T_{\text{max}}$)
* **Objetivo:** Snapshot consolidado do estado do sistema ao final da simulação para re-arranque e análise de perigo de impacto.
* **Critérios de Ordenação:**
  * `encuent.csv`: Ordenado por **`dist_km` crescente** (as menores distâncias aparecem primeiro no topo). Isso permite ao astrofísico identificar instantaneamente os clones mais rasantes e perigosos (ex: Clone 28/22 em 2032).
  * `orbitas.csv`: Ordenado por **`N`** (corpos de 1 a 408 no instante $T=20$).
  * `condin.csv`: Ordenado por **`N`** (substituto direto do binário `condin.xxx` para re-arranque da simulação).
  * `colisio.csv` e `finales.csv`: Relatórios consolidados de integridade (0 colisões e 0 ejeções).

---

## 📁 4. Tabela Resumo dos 7 Arquivos Gerados

As saídas são geradas em formato **CSV limpo, padronizado e com cabeçalho**, divididas em duas pastas:

```
data/rebound_output/<pasta>/
├── 📂 final/       -> Recorte estático no instante final (T = T_max)
└── 📂 temporal/    -> Série histórica completa desde o início (T = 0 até T = T_max)
```

### Detalhamento dos 7 Arquivos:

| # | Arquivo CSV | Equivalente EVORB | O que contém na pasta `final/` | O que contém na pasta `temporal/` |
| :-: | :--- | :--- | :--- | :--- |
| 1 | **`orbitas.csv`** | `orbitas.sal` | Elementos keplerianos osculadores ($a, e, i, \Omega, \omega, M$) no tempo final. | Evolução temporal completa de cada corpo desde $T=0$. |
| 2 | **`condin.csv`** | `condin.xxx` | Coordenadas cartesianas ($X, Y, Z, Vx, Vy, Vz$) no final (para re-arranque). | Histórico de vetores de estado cartesianos heliocêntricos. |
| 3 | **`encuent.csv`** | `encuent.sal` | Encontros mais próximos com a Terra ordenados por menor distância. | Cronograma de passagens rasantes na ordem em que ocorreram. |
| 4 | **`colisio.csv`** | `colisio.sal` | Balanço consolidado de impactos físicos com planetas. | Linha do tempo de colisões físicas com data e velocidade. |
| 5 | **`finales.csv`** | `finales.sal` | Balanço de partículas eliminadas (queda no Sol ou ejeção além de 100 AU). | Linha do tempo de cada eliminação no momento exato. |
| 6 | **`energia.csv`** | `energia.sal` | Resumo final do erro relativo de energia $\Delta E / E$. | Curva temporal de energia cinética, potencial e total. |
| 7 | **`xyzv.csv`** | `xyzv.sal` | Snapshot cartesiano final de todos os corpos. | Série temporal completa de trajetórias em 3D. |

---

## 🌌 4. Como Simular Qualquer Outro Asteroide

Para simular qualquer asteroide ou cometa, basta criar ou editar o arquivo `datos.dat` na pasta de entrada desejada.

### A Estrutura do Arquivo `datos.dat`:
O arquivo é composto por parâmetros gerais, seguidos pela tabela dos planetas e pela tabela das partículas:

```text
*INDICE CORRECCION RELATIVISTA (0 NO, 1 Gallardo-Venturini, 2 standard)
1
*INDICE PARA GUARDAR (1) O NO (0) ARCHIVO CON (X,Y,Z,Vx,Vy,Vz)
1
*LIMITES INFERIOR Y SUPERIOR EN UAs PARA ELIMINACION DE CUERPOS 
0.1  100.0
*MASA CENTRAL EN MASAS SOLARES Y RADIO EN kms
1.0  696000.0
*TIEMPO MAXIMO DE LA SIMULACION EN AÑOS:
20.0
*INTERVALO DE ESCRITURA EN AÑOS (multiplo del PASO):
0.01
*PASO DE INTEGRACION en años:
0.0001
*INDICE PARA RE-ARRANQUE (0 ES UNA CORRIDA NUEVA, 1 ES CONTINUACION)
0
*NUMERO DE PLANETAS (minimo uno)
8
*JD epoca de planetas
2451400.5
*elementos J2000 de los planetas: N a e i nodo argper anomedia masa radio(km)
1  0.38709897  0.205622188  7.0050053  48.3307302  29.1239378 303.4539666 1.66013679D-07  2440.0 
2  0.72333793  0.006771606  3.3945712  76.6821802  54.8201213 178.9704486 2.44783833D-06  6052.0 
3  0.99998836  0.016736826  0.0001213  27.6464539  75.3942093 215.0034486 3.04043264D-06  6371.0
4  1.52367503  0.093290001  1.8498920  49.5616397 286.5090192 303.6562023 3.22715144D-07  3390.0
5  5.20399661  0.048719684  1.3046775 100.4829904 275.0927233   6.8013147 9.54791938D-04 69911.0
6  9.58195033  0.055272778  2.4852156 113.6436321 335.5040835 315.9800041 2.85885980D-04 58232.0
7 19.24744286  0.043713914  0.7727772  74.0128993  97.4115217 140.3579592 4.36624404D-05 25362.0
8 30.14125054  0.011166945  1.7676331 131.7957546 259.2549938 273.2914807 5.15138902D-05 24624.0
*NUMERO DE PARTICULAS (minimo una)
1
*JD epoca de particulas
2460600.5
*elementos J2000 de ESAS particulas: N a e i nodo argper anomedia
9  2.539027  0.664131  3.4528  271.4123  134.6422  351.0758
```

---

### Exemplo 1: Simulando o Asteroide 99942 Apophis
O asteroide **Apophis** terá um encontro próximo histórico em **13 de abril de 2029** a apenas ~31.000 km da superfície da Terra.

Para simular o Apophis:
1. Obtenha os 6 elementos keplerianos no [NASA JPL Small-Body Database](https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=99942):
   * Época JD: `2460600.5`
   * $a = 0.922438\text{ AU}$
   * $e = 0.191503$
   * $i = 3.3314^\circ$
   * $\Omega = 204.0416^\circ$
   * $\omega = 126.3986^\circ$
   * $M = 230.8271^\circ$
2. Edite o bloco de partículas no `datos.dat`:
   ```text
   *NUMERO DE PARTICULAS (minimo una)
   1
   *JD epoca de particulas
   2460600.5
   *elementos J2000 de ESAS particulas: N a e i nodo argper anomedia
   9  0.922438  0.191503  3.3314  204.0416  126.3986  230.8271
   ```
3. Defina o tempo de simulação para `10.0` anos e execute:
   ```bash
   python scripts/run_simulation.py -i pasta_apophis -y 10.0
   ```
4. Abra o arquivo `pasta_apophis/final/encuent.csv`: ele registrará a aproximação de Abril de 2029 a menos de $6 R_\oplus$ com velocidade relativa de ~7.4 km/s!

---

### Exemplo 2: Simulando o Asteroide 101955 Bennu
Para simular o asteroide **Bennu** (alvo da missão OSIRIS-REx):
1. Elementos J2000:
   * $a = 1.126391\text{ AU}$
   * $e = 0.203745$
   * $i = 6.0349^\circ$
   * $\Omega = 2.0609^\circ$
   * $\omega = 66.2231^\circ$
   * $M = 101.7039^\circ$
2. Coloque esses valores na linha da partícula e execute por 50 ou 100 anos.

---

## 📊 5. Análise de Dados e Visualização

Como todos os outputs são arquivos CSV regulares com cabeçalhos expressivos, você pode abri-los diretamente em qualquer ferramenta:

### No Python (Pandas & Matplotlib)
```python
import pandas as pd
import matplotlib.pyplot as plt

# 1. Carrega a evolução orbital no tempo
df_temp = pd.read_csv("data/rebound_output/datos1/temporal/orbitas.csv")

# 2. Filtra o Clone 9
clone9 = df_temp[df_temp["nome"] == "clone_9"]

# 3. Plota a transição do semi-eixo maior provocada pelo encontro de 2032
plt.figure(figsize=(10, 5))
plt.plot(clone9["t_anos"], clone9["a"], label="Semi-eixo maior (a)", color="crimson")
plt.axvline(x=8.18, color="navy", linestyle="--", label="Encontro com a Terra (Dez/2032)")
plt.xlabel("Tempo (anos)")
plt.ylabel("Semi-eixo a (AU)")
plt.title("Evolução Orbital do Asteroide 2024 YR4")
plt.legend()
plt.grid(True)
plt.show()
```

### No Excel ou Power BI
Basta clicar duas vezes em qualquer arquivo da pasta `final/` ou `temporal/`. Eles abrem nativamente com colunas separadas e valores formatados.
