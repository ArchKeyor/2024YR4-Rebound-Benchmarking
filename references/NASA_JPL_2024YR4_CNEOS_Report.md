# Relatório Técnico Oficial: Asteroide 2024 YR4 (NASA JPL / CNEOS)

> **Fonte Oficial:** NASA Jet Propulsion Laboratory (JPL) Solar System Dynamics (SSD) Group & Center for Near Earth Object Studies (CNEOS).
> **Status de Dados:** Consultado via JPL SBDB API & CNEOS CAD API.

---

## 1. Identificação e Parâmetros Astrométricos

- **Designação Oficial:** `(2024 YR4)`
- **Classe Orbital:** Apollo (Near Earth Asteroid)
- **Solução Orbital (Orbit ID):** `428`
- **Época da Órbita:** 2461200.5 JD
- **Arco de Observação (Data Arc):** 428 dias
- **Primeira Observação:** 2024-12-25
- **Última Observação:** 2026-02-26
- **Condição de Incerteza (Condition Code / U):** `1` (0 = órbita muito bem determinada, 9 = altamente incerta)
- **RMS de Ajuste:** `0.45` segundos de arco

---

## 2. Elementos Orbitais Keplerianos Osculantes

| Elemento | Símbolo | Valor | Incerteza (1σ) | Unidade | Descrição |
| :--- | :---: | :---: | :---: | :---: | :--- |
| eccentricity | `e` | `0.661` | `3.5e-08` | None | eccentricity |
| semi-major axis | `a` | `2.52` | `2.7e-07` | au | semi-major axis |
| perihelion distance | `q` | `0.852` | `2.7e-08` | au | perihelion distance |
| inclination; angle with respect to x-y ecliptic plane | `i` | `3.41` | `9.9e-07` | deg | inclination; angle with respect to x-y ecliptic plane |
| longitude of the ascending node | `om` | `271` | `8.1e-06` | deg | longitude of the ascending node |
| argument of perihelion | `w` | `134` | `8.1e-06` | deg | argument of perihelion |
| mean anomaly | `ma` | `139` | `2.2e-05` | deg | mean anomaly |
| time of perihelion passage | `tp` | `2460636.817` | `3.2e-06` | TDB | time of perihelion passage |
| sidereal orbital period | `per` | `1460` | `0.00023` | d | sidereal orbital period |
| mean motion | `n` | `0.247` | `3.9e-08` | deg/d | mean motion |
| aphelion distance | `ad` | `4.18` | `4.4e-07` | au | aphelion distance |

---

## 3. Propriedades Físicas Estimadas

| Propriedade | Valor | Unidade | Descrição |
| :--- | :---: | :---: | :--- |
| absolute magnitude | `23.93` | None | absolute magnitude (magnitude at 1 au from Sun and observer) |

---

## 4. Tabela de Aproximações Próximas (Close Approaches - CNEOS CAD)

| Data / Hora (TDB) | Distância Nominal (AU) | Distância Mínima (AU) | Distância Máxima (AU) | Distância Nominal (LD / km) | Velocidade Relativa ($v_{rel}$) | Velocidade Infinita ($v_{\infty}$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **2024-Dec-25 04:46** | `0.0055400` | `0.00553999118089998` | `0.00554003939302755` | **2.16 LD** (~828,774 km) | `13.7514611451871 km/s` | `13.7164419724507 km/s` |
| **2032-Dec-22 08:22** | `0.0018605` | `0.00184602884773144` | `0.00187505441666802` | **0.72 LD** (~278,333 km) | `13.3071686728549 km/s` | `13.1991113554716 km/s` |

> *Nota:* 1 LD (Distância Lunar média) $\approx 384.400\text{ km} \approx 0.00257\text{ AU}$.

---

## 5. Histórico no Sistema Sentry de Risco de Impacto

### A Situação de 2032 e a Remoção Oficial
- **Descoberta:** Final de dezembro de 2024.
- **Inclusão no Sentry:** Inicialmente, devido ao arco curto de observação, a elipse de incerteza para o encontro de **22 de dezembro de 2032** intersectava a seção de choque gravitacional da Terra, registrando probabilidade de impacto superior a 1% e nível na Escala de Palermo / Turim.
- **Campanha Astrométrica Global:** Observações subsequentes de telescópios ópticos e de radar refinaram os parâmetros orbitais, estreitando drasticamente o semieixo maior da elipse de covariância.
- **Remoção do Sentry:** 2026-03-05 16:02:50.
- **Resultado Final:** A distância nominal calculada para 2032-12-22 é de **0.00186 AU (~278.330 km / 0.72 LD)**, confirmando passagem segura sem risco de colisão.

---

## 6. Relevância para o Projeto OrbitREBOUND

1. **Validação Orbital:** As efemérides geradas pelo simulador `orbit-rebound` (com integrador IAS15 e o legado EVORB15) utilizam os vetores de estado de 2024 YR4 para reproduzir a trajetória de aproximação de 2024 e projetar a dinâmica de 2032.
2. **Sensibilidade a Perturbações Planetárias:** Por passar dentro da esfera de influência do sistema Terra-Lua (~0.72 LD), a integração exige integrador de 15ª ordem com passo adaptativo (IAS15) para resolver a aceleração diferencial lunar e terrestre sem degenerar a conservação de energia e momento angular.
