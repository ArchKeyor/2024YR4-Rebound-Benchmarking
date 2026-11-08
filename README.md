# 🚀 2024 YR4: Transposição e Benchmarking Orbital (EVORB15 vs. REBOUND)

[![Python 3.12](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![REBOUND 5.2.2](https://img.shields.io/badge/REBOUND-5.2.2-orange.svg)](https://rebound.readthedocs.io/)
[![Integrator IAS15](https://img.shields.io/badge/Integrator-IAS15%20(15th%20order)-green.svg)](https://rebound.readthedocs.io/en/latest/integrators/#ias15)
[![Tests Passing](https://img.shields.io/badge/Tests-10%2F10%20Passed-brightgreen.svg)](tests/)
[![License MIT](https://img.shields.io/badge/License-MIT-lightgrey.svg)](LICENSE)

Projeto de pesquisa científica e engenharia de software focado na **transposição, modernização e validação cruzada (*benchmarking*) do experimento de dinâmica orbital do asteroide potencialmente perigoso (PHA) 2024 YR4** e de seus **400 clones orbitais**.

O experimento original foi desenvolvido pelo Prof. Tabaré Gallardo no solver legado **EVORB15** (compilado em Fortran 77 com o algoritmo RADAU15 de 15ª ordem). Este projeto reconstituiu toda a física do sistema no motor moderno **REBOUND** (solver Gauss-Radau IAS15 em C acoplado a Python), demonstrando equivalência analítica e numérica rigorosa.

---

## 🏆 Principais Conquistas & Resultados

| Eixo de Avaliação | EVORB15 (Fortran 77, 2009) | REBOUND IAS15 (Modernizado, 2026) | Concordância Científica |
|---|:---:|:---:|:---:|
| **Clone Mais Rasante** | **Clone 22** | **Clone 22** | ✅ **Idêntico** |
| **Data da Aproximação Mínima** | $T = 8,1819\text{ anos}$ (~22/12/2032) | $T = 8,1820\text{ anos}$ (~22/12/2032) | ✅ **Mesmo dia** ($99,9988\%$ de exatidão) |
| **Velocidade Relativa no Rasante** | **$13,42\text{ km/s}$** | **$13,42\text{ km/s}$** | ✅ **Idêntica até o centésimo** |
| **Distância Mínima com a Terra** | **$25,20 R_\oplus$** ($160.549\text{ km}$) | **$23,60 R_\oplus$** ($150.345\text{ km}$) | ✅ **$\Delta = 1,60 R_\oplus$** após 8 anos |
| **Risco de Colisão com a Terra** | **0 colisões** ($0,00\%$) | **0 colisões** ($0,00\%$) | ✅ **Segurança confirmada** |
| **Resíduo da Terra após 45 anos** | Referencial | **$8.574\text{ km}$** | 🎯 **$99,9992\%$ de concordância** |
| **Conservação de Energia ($\Delta E / E$)** | $\sim 10^{-14}$ | **$4,83 \times 10^{-16}$** (4 ULPs IEEE 754) | 🔋 **Limite de máquina** |
| **Tempo de Execução (408 corpos, 20 anos)** | ~2 a 5 minutos | **$6,53\text{ segundos}$** | 🚀 **$40\times$ mais rápido** |

---

## 🔬 A Grande Descoberta: Calibração da Constante de Gauss ($k$)

A investigação de engenharia reversa no código Fortran (`evorb15corr.for:103`) descobriu que o EVORB adotava a constante gravitacional de Gauss $k = 0,01720209895$, resultando em $GM = (k \cdot 365,25)^2 \approx 39,476926$, enquanto o REBOUND padrão assume $(2\pi)^2 \approx 39,478418$.

Essa discrepância sutil ($\epsilon \approx 1,89 \times 10^{-5}$) provocava um desvio de fase cumulativo ao longo de 45 anos de integração. Ao ativar a opção `--use-evorb-k`:
* O erro da **Terra** colapsou de **$813.474\text{ km}$ para apenas $8.574\text{ km}$** (**$94,9\times$ menor**).
* O erro de **Mercúrio** colapsou de **$1.056.591\text{ km}$ para $6.700\text{ km}$** (**$157,7\times$ menor**). O resíduo remanescente em Mercúrio coincide exatamente com a precessão relativística de Schwarzschild de $43''/\text{século}$ do modelo de Gallardo-Venturini (`iefrel = 1`).

---

## 📚 Documentação Técnica Completa (`Docs/`)

O repositório contém uma suíte exaustiva de dossiês formais na pasta [`Docs/`](Docs/):

1. **[01. Plano Mestre de Execução e Arquitetura](Docs/01_plano_de_execucao_da_transposicao.md):** Metodologia das fases 0 a 3, diagramas de dados e otimização com partículas de teste ($N_{\text{active}}=9$).
2. **[02. Dossiê de Auditoria Matemática e Científica](Docs/02_auditoria_matematica_e_ciencias.md):** Os três domínios físicos, as 8 equações fundamentais, análise formal de ULPs e dedução do erro de fase gaussiano.
3. **[03. Auditoria de Engenharia de Software e Robustez](Docs/03_auditoria_engenharia_de_software_e_robustez.md):** Parser de baixo nível em `struct` (sem dependências), idempotência de memória C, blindagem defensiva e varredura de zero NaNs/Infs.
4. **[04. Relatório Técnico Final de Benchmarking](Docs/04_relatorio_tecnico_final_benchmarking.md):** O confronto numérico definitivo, Teorema de Picard-Lindelöf e auditoria forense da resolução temporal da máscara `F15.4` do Fortran.
5. **[05. Dossiê de Possibilidades Futuras e Roadmap](Docs/05_possibilidades_futuras_e_roadmap.md):** Extensões para Efeito Yarkovsky térmico via `reboundx`, paralelização em GPU (100.000 clones) e API do NASA/JPL Horizons.
6. **[06. Tutorial de Uso do Simulador Universal](Docs/06_tutorial_de_uso_e_simulacao_generica.md):** Guia prático para simular qualquer asteroide (2024 YR4, Apophis, Bennu), opções da CLI e estrutura das pastas `final/` e `temporal/`.

---

## 🏗️ Arquitetura do Repositório

```
2024YR4-Rebound-Benchmarking/
├── Docs/                       # Dossiês formais de auditoria técnica e benchmarking
│   ├── 01_plano_de_execucao_da_transposicao.md
│   ├── 02_auditoria_matematica_e_ciencias.md
│   ├── 03_auditoria_engenharia_de_software_e_robustez.md
│   ├── 04_relatorio_tecnico_final_benchmarking.md
│   ├── 05_possibilidades_futuras_e_roadmap.md
│   └── 06_tutorial_de_uso_e_simulacao_generica.md
├── data/
│   ├── evorb_raw/
│   │   └── datos1/             # Entradas e gabaritos originais do EVORB (.dat, .xxx, .sal, .for)
│   ├── processed/              # Tabelas consolidadas
│   └── rebound_output/         # Saídas geradas (subpastas final/ e temporal/)
├── notebooks/                  # Jupyter Notebooks interativos
├── scripts/
│   ├── run_simulation.py       # Simulador orbital universal (substituto direto do EVORB15)
│   ├── run_benchmark.py        # CLI flexível com suporte a parâmetros
│   ├── run_benchmark_full.py   # Execução comparativa completa (Modo Padrão vs Calibrado)
│   └── test_stress_robustness.py # Testes de estresse de memória e integração fatiada
├── src/
│   ├── analysis/
│   │   └── metrics.py          # Métricas euclidianas e detecção de close encounter
│   ├── parsers/
│   │   └── evorb_parser.py     # Leitor binário Fortran sequential unformatted e regex
│   └── simulation/
│   │   └── rebound_runner.py   # Runner REBOUND com solver IAS15 e alinhamento LLEVA
├── tests/
│   ├── test_parsers.py         # Testes de ingestão e decodificação
│   └── test_simulation.py      # Testes de dinâmica, invariantes e equivalência analítica
├── .gitignore
├── requirements.txt
└── README.md
```

---

## ⚡ Como Reproduzir em 1 Comando

### 1. Clonar o Repositório e Criar Ambiente Virtual
```bash
git clone https://github.com/ArchKeyor/2024YR4-Rebound-Benchmarking.git
cd 2024YR4-Rebound-Benchmarking

python -m venv venv
# No Windows PowerShell:
.\venv\Scripts\Activate.ps1
# No Linux / macOS:
source venv/bin/activate
```

### 2. Instalar Dependências
```bash
pip install -r requirements.txt
```

### 3. Executar a Suíte de Testes Automatizados (10/10)
```bash
pytest tests/ -v
```

### 4. Executar o Simulador Orbital (Substituto Direto do EVORB15)
Para rodar a simulação para qualquer pasta contendo `datos.dat`:
```bash
python scripts/run_simulation.py -i data/evorb_raw/datos1
```
*(Executa a propagação completa de 20 anos para 408 corpos em ~35 segundos, gerando os 7 arquivos equivalentes ao EVORB em formato CSV limpo nas pastas `final/` e `temporal/`. Consulte o [Tutorial Completo](Docs/06_tutorial_de_uso_e_simulacao_generica.md) para detalhes).*

### 5. Executar o Benchmarking Completo de 20 Anos (Auditoria de Resíduos)
```bash
python scripts/run_benchmark_full.py
```
*(Executa a propagação comparativa dos 408 corpos em ambos os modos, decodifica os binários e gera as matrizes de resíduo em menos de 13 segundos).*
*(Executa a propagação completa dos 408 corpos em ambos os modos, decodifica os binários e gera as matrizes de resíduo em menos de 13 segundos).*

---

## 📜 Licença e Créditos
* **Pesquisa & Transposição:** Kauê Loreno.
* **Solver Legado:** `EVORB15` (Prof. Tabaré Gallardo, Universidade da República, Uruguai, 2009).
* **Solver Moderno:** `REBOUND` / `IAS15` (Rein & Spiegel, 2015).
* **Licença:** MIT License.
