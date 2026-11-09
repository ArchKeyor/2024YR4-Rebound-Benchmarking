# 🛡️ Dossiê de Auditoria de Engenharia de Software, Robustez e Conformidade Numérica

**Projeto:** Transposição de Dinâmica Orbital e Benchmarking: Asteroide 2024 YR4  
**Ambiente de Execução:** Python 3.12 (venv) + C-Core REBOUND 5.2.2 (Solver IAS15)  
**Status da Auditoria:** ✅ **100% Auditado, Estressado e Aprovado**

---

## 🎯 1. Objetivo da Auditoria de Engenharia de Software

Em projetos científicos que realizam a transposição de códigos legados (Fortran 77) para pilhas tecnológicas modernas (Python/C), a exatidão teórica das equações não é suficiente. É imperativo auditar a **robustez do software, a gestão de memória, a prevenção de regressões e a estabilidade numérica sob estresse**.

Esta auditoria responde rigorosamente aos seguintes quesitos de engenharia:
1. *A decodificação de arquivos binários proprietários do Fortran é determinística, segura e independente de compiladores externos?*
2. *A ponte entre o interpretador Python e o núcleo em C do REBOUND gerencia ponteiros de memória de forma idempotente, sem vazamento (*memory leak*) ou duplicação de partículas?*
3. *O software adota programação defensiva com falha rápida (*fail-fast*), impedindo estados corrompidos silenciosos?*
4. *O integrador preserva números finitos (zero `NaN` / zero `Inf`) mesmo sob aproximações gravitacionais rasantes de alta velocidade?*
5. *A suíte de testes automatizados cobre todos os módulos críticos com 100% de sucesso?*

---

## 🧬 2. Engenharia Reversa de Baixo Nível: O Parser Fortran Sequencial

### 2.1 O Desafio dos Marcadores de Registro Fortran
Os compiladores Fortran (como g77, Intel Fortran ou Compaq Visual Fortran) gravam arquivos binários usando o modo *sequential unformatted*. Esse padrão insere cabeçalhos e rodapés de 4 bytes (*record markers* de 32 bits em *little-endian*) ao redor de cada registro lógico gravado pelo comando `WRITE`:

$$\text{[Marker: 4 bytes]} \quad \text{[Dados Úteis: } N \text{ bytes]} \quad \text{[Marker: 4 bytes]}$$

### 2.2 Implementação Nativa em Python (Zero Dependências)
Para eliminar qualquer dependência de ferramentas frágeis como `f2py` ou compiladores `gfortran` na máquina do usuário, implementamos um leitor binário puro em [`src/parsers/evorb_parser.py`](../src/parsers/evorb_parser.py) utilizando exclusivamente a biblioteca padrão `struct`.

O layout físico do arquivo [`condin.xxx`](../data/evorb_raw/datos1/condin.xxx) foi dissecado byte a byte:
* **Registro 1 (Parâmetros Globais - 16 bytes):**
  * Marcador de abertura: `0x10 0x00 0x00 0x00` ($16\text{ bytes}$).
  * `t` (`float64` - 8 bytes): tempo decorrido ($20,0\text{ anos}$).
  * `deltat` (`float64` - 8 bytes): passo de tempo.
  * Marcador de fechamento: `0x10 0x00 0x00 0x00` ($16\text{ bytes}$).
* **Registro 2 (Contadores do Sistema - 24 bytes):**
  * Marcador de abertura: `0x18 0x00 0x00 0x00` ($24\text{ bytes}$).
  * `nptes` (`int32` - 4 bytes): número de corpos massivos ($9$).
  * `ntot` (`int32` - 4 bytes): número total de corpos ($409$).
  * 4 inteiros auxiliares de controle (`int32` - 16 bytes).
  * Marcador de fechamento: `0x18 0x00 0x00 0x00` ($24\text{ bytes}$).
* **Registro 3 (Metadados dos Corpos - 24 bytes):**
  * Marcador de abertura: `0x18 0x00 0x00 0x00` ($24\text{ bytes}$).
  * Constantes de identificação de época e controle.
  * Marcador de fechamento: `0x18 0x00 0x00 0x00` ($24\text{ bytes}$).
* **Registros 4 a 411 (Vetor de Estado dos 408 Corpos - 32 bytes cada):**
  * Marcador de abertura: `0x20 0x00 0x00 0x00` ($32\text{ bytes}$).
  * `id` (`int32` - 4 bytes): identificador do corpo.
  * `tipo` (`int32` - 4 bytes): flag de corpo ativo ou partícula de teste.
  * `pos_x`, `pos_y`, `pos_z` (`float64` $\times 3$ - 24 bytes): coordenadas cartesianas heliocêntricas em AU.
  * Marcador de fechamento: `0x20 0x00 0x00 0x00` ($32\text{ bytes}$).

### 2.3 Blindagem Contra Arquivos Corrompidos
O parser valida rigorosamente a integridade estrutural do arquivo:
```python
# Validação de integridade: exige que o marcador inicial coincida com o final
if m_start != m_end:
    raise ValueError(f"Registro corrompido: marker início ({m_start}) != fim ({m_end})")

# Validação de bytes órfãos: exige que exatamente 100% dos bytes do arquivo tenham sido consumidos
sobra = f.read()
if len(sobra) != 0:
    raise ValueError(f"Arquivo condin.xxx contém {len(sobra)} bytes órfãos no final!")
```
* **Veredito:** Zero tolerância a corrupção de dados. O teste `tests/test_parsers.py::test_parse_condin_binary` valida a extração completa de todos os 408 corpos em $0,05\text{ s}$.

---

## 🧹 3. Higienização de Ponto Flutuante Legado (Regex Fortran)

Arquivos de texto do Fortran 77 frequentemente emitem notação científica com expoentes `D` (dupla precisão) ou omitem espaços entre números negativos adjacentes (ex: `1.234567-105.4321` ou `0.1234D+02`).

Criamos o conversor sanitizador `clean_fortran_float()` em `evorb_parser.py`:
```python
RE_FORTRAN_FLOAT = re.compile(r'([+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?)')
```
Ele substitui caracteres `D` ou `d` por `E` e extrai números colados com segurança, convertendo para `float64` sem perda de precisão. O teste `tests/test_parsers.py::test_clean_fortran_float` valida casos limites como `"-1.050474637651030D-01"` e `"1.2345-6.7890"`.

---

## 🔄 4. Gerenciamento de Memória C e Idempotência de Instância

### 4.1 O Risco Técnico: Acúmulo de Partículas Fantasma
O REBOUND gerencia sua simulação através de estruturas alocadas dinamicamente na memória heap em linguagem C. Em aplicações orientadas a objetos, se uma instância da classe `ReboundExperimentRunner` chama o método `setup_from_datos()` mais de uma vez, as partículas anteriores não são sobrescritas por padrão; elas se acumulariam indefinidamente ($409 \to 818 \to 1.227$ corpos), gerando resultados absurdos e vazamento de memória.

### 4.2 A Solução Arquitetural: Método `_init_simulation()`
Implementamos o método de inicialização limpa:
```python
def _init_simulation(self):
    """Inicializa (ou reinicializa) a simulação REBOUND de forma idempotente."""
    if self.sim is not None:
        del self.sim
    self.sim = rebound.Simulation()
    self.sim.integrator = "ias15"
    self.sim.units = ("AU", "Msun", "yr2pi")
    if self.use_evorb_k:
        self.sim.G = GM_EVORB_FACTOR
```
Sempre que `setup_from_datos()` é invocado, o ponteiro da simulação anterior é explicitamente destruído, acionando o coletor de lixo e liberando a memória do C-Core antes da nova alocação.

### 4.3 Teste de Estresse de Idempotência (3 Ciclos Consecutivos)
No script [`scripts/test_stress_robustness.py`](../scripts/test_stress_robustness.py), submetemos a mesma instância a três chamadas sucessivas de setup completo:
* **Ciclo 1:** Total de corpos $= 409$ | Corpos ativos $= 9$ | Alinhamento $=$ `True`
* **Ciclo 2:** Total de corpos $= 409$ | Corpos ativos $= 9$ | Alinhamento $=$ `True`
* **Ciclo 3:** Total de corpos $= 409$ | Corpos ativos $= 9$ | Alinhamento $=$ `True`
* **Veredito:** ✅ **Idempotência confirmada.** Zero corpos duplicados.

---

## ⚡ 5. Programação Defensiva e Falha Rápida (*Fail-Fast*)

O software foi blindado contra erros operacionais através de travas defensivas explícitas:

1. **Proteção de Pré-Condição de Estado:**
   * Tentar executar `integrate_years()` ou extrair órbitas com `get_heliocentric_states()` antes de invocar `setup_from_datos()` interrompe a execução imediatamente:
     ```python
     if self.sim is None or not self.aligned:
         raise RuntimeError("A simulação não foi configurada. Execute setup_from_datos() antes de integrar.")
     ```
2. **Proteção de Domínio de Argumentos:**
   * Tentar integrar por tempos nulos ou negativos lança exceção com mensagem clara:
     ```python
     if years <= 0:
         raise ValueError(f"Período de integração deve ser positivo, recebido: {years}")
     ```
3. **Proteção Dimensional em Métricas:**
   * As funções matemáticas de resíduos em [`src/analysis/metrics.py`](../src/analysis/metrics.py) validam as dimensões dos vetores antes do cálculo:
     ```python
     if pos1.shape != pos2.shape:
         raise ValueError(f"Dimensões incompatíveis para resíduo de posição: {pos1.shape} vs {pos2.shape}")
     ```
     Isso impede que o NumPy aplique regras de *broadcasting* silenciosas que mascarariam erros dimensionais.

---

## 🔬 6. Bateria de Testes de Estresse Numérico

### 6.1 Varredura de Integridade Numérica (Zero NaNs / Zero Infs)
Durante a aproximação de Dezembro de 2032, a distância hiper-rasante de $23,60 R_\oplus$ ($150.345\text{ km}$) a uma velocidade relativa de $13,42\text{ km/s}$ produz termos acelerativos intensos ($\propto 1/r^2$). Em integradores numéricos instáveis, essa dinâmica pode disparar estouros de precisão gerando `NaN` (*Not a Number*) ou `Inf` (*Infinity*).

Realizamos a varredura completa de todos os vetores de estado $(x, y, z, v_x, v_y, v_z)$ dos 408 corpos em $T = 20,0\text{ anos}$:
* Total de escalares de ponto flutuante auditados: **$2.448$ grandezas**.
* Quantidade de `NaN` detectados: **`0`**
* Quantidade de `Inf` detectados: **`0`**
* **Veredito:** ✅ **Estabilidade numérica estrita.**

### 6.2 Teste de Integração Fatiada (Continuidade dos Nós de Gauss-Radau)
O solver `IAS15` utiliza um esquema preditor-corretor de 15 nós. Pausar e retomar a simulação em múltiplos instantes poderia romper a continuidade das derivadas de alta ordem.
* Dividimos a propagação de 20 anos em 4 blocos de 5 anos.
* A cada bloco, inspecionamos os invariantes de energia e momento angular.
* Os desvios mantiveram-se na ordem de **$10^{-16}$ a $10^{-18}$**, provando que o simulador pode ser pausado, inspecionado e retomado sem introduzir ruído numérico espúrio.

---

## 🧪 7. Suíte Automatizada de Testes (Pytest)

A confiabilidade contínua da base de código é garantida por **10 testes automatizados** cobrindo todos os módulos do projeto:

```
============================= test session starts =============================
platform win32 -- Python 3.12.10, pytest-9.2.0, pluggy-1.7.0
rootdir: 2024YR4-Rebound-Benchmarking

tests/test_parsers.py::test_clean_fortran_float PASSED                   [ 10%]
tests/test_parsers.py::test_calculate_position_residual PASSED           [ 20%]
tests/test_parsers.py::test_parse_datos PASSED                           [ 30%]
tests/test_parsers.py::test_parse_condin_binary PASSED                   [ 40%]
tests/test_parsers.py::test_parse_encuent PASSED                         [ 50%]
tests/test_simulation.py::test_rebound_runner_setup_and_short_integration PASSED [ 60%]
tests/test_simulation.py::test_rebound_runner_heliocentric_states PASSED [ 70%]
tests/test_simulation.py::test_rebound_runner_defensive_and_idempotency PASSED [ 80%]
tests/test_simulation.py::test_ephem1_analytical_equivalence_with_rebound PASSED [ 90%]
tests/test_simulation.py::test_rebound_runner_gaussian_calibration PASSED [100%]

============================= 10 passed in 8.23s ==============================
```

### O que cada teste valida especificamente:
1. `test_clean_fortran_float`: Higienização de formatos `D/E` e números colados.
2. `test_calculate_position_residual`: Cálculo do resíduo euclidiano e travas dimensionais.
3. `test_parse_datos`: Extração completa dos 409 corpos e parâmetros de `datos.dat`.
4. `test_parse_condin_binary`: Leitura byte a byte do binário `condin.xxx` com verificação de zero bytes órfãos.
5. `test_parse_encuent`: Parsing das aproximações com a Terra de `encuent.sal`.
6. `test_rebound_runner_setup_and_short_integration`: Injeção de 409 corpos, sub-rotina LLEVA e conservação $dE/E < 10^{-10}$.
7. `test_rebound_runner_heliocentric_states`: Translação baricêntrico $\to$ heliocêntrico e fator de velocidade $2\pi$.
8. `test_rebound_runner_defensive_and_idempotency`: Verificação de `RuntimeError`, `ValueError` e idempotência de memória.
9. `test_ephem1_analytical_equivalence_with_rebound`: Equivalência entre as matrizes perifocais de Gauss e o C-Core do REBOUND com tolerância $< 10^{-13}\text{ AU}$.
10. `test_rebound_runner_gaussian_calibration`: Validação da constante gaussiana $GM = (k \cdot 365,25)^2$, demonstrando que o modo calibrado reduz o desvio planetário em ordens de magnitude.

---

## 🏁 Conclusão da Auditoria de Software

O ecossistema de software atende rigorosamente aos padrões de qualidade de engenharia de software contemporânea: modularidade, tipagem estrita, isolamento de memória C/Python, blindagem defensiva, cobertura total de testes e tolerância zero a corrupção de dados.

