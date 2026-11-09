"""
Módulo de Simulação N-Corpos com REBOUND (IAS15).
Transposição e execução do experimento do asteroide 2024 YR4.
"""

from typing import List, Dict, Any, Optional
import numpy as np
import rebound

# Constante astronômica de conversão: 1 AU em km
AU_IN_KM = 149597870.7

# Constantes de calibração gaussiana do EVORB15 (evorb15corr.for:103)
GAUSS_K = 0.01720209895
JULIAN_YEAR_DAYS = 365.25
GM_EVORB_FACTOR = (GAUSS_K * JULIAN_YEAR_DAYS / (2.0 * np.pi)) ** 2


class ReboundExperimentRunner:
    """
    Gerenciador da simulação N-corpos em REBOUND espelhando rigorosamente
    o modelo físico do EVORB15 (integradora IAS15 de 15ª ordem).
    
    Características de Engenharia:
    - Idempotência: chamadas repetidas de setup limpam a memória C sem vazamentos.
    - Validação Defensiva: impede integração ou consultas em estado não-inicializado.
    - Calibração de G: suporta modo padrão (G=1.0) e modo gaussiano EVORB (GM = (k*365.25)^2).
    - Tipagem Estrita e Blindagem Numérica.
    """

    def __init__(self, use_evorb_k: bool = False):
        """Inicializa o simulador em estado limpo."""
        self.use_evorb_k: bool = use_evorb_k
        self.sim: rebound.Simulation = rebound.Simulation()
        self.t_alinhamento: float = 0.0
        self.alinhado: bool = False
        self._init_simulation()

    def _init_simulation(self):
        """Configura a instância do REBOUND com unidades e integrador padrão."""
        self.sim = rebound.Simulation()
        self.sim.integrator = "ias15"
        self.sim.units = ("AU", "Msun", "yr2pi")
        if self.use_evorb_k:
            self.sim.G = GM_EVORB_FACTOR
        self.alinhado = False

    def setup_from_datos(
        self,
        params_or_datos: Dict[str, Any],
        planets: Optional[List[Dict[str, Any]]] = None,
        particles: Optional[List[Dict[str, Any]]] = None,
        use_evorb_k: Optional[bool] = None
    ):
        """
        Carrega as condições iniciais a partir do parser do datos.dat:
        1. Reinicializa a simulação para garantir estado limpo (idempotência).
        2. Adiciona o Sol no baricentro inicial com seu raio físico.
        3. Adiciona os 8 planetas com elementos osculadores na época J2000.
        4. Executa o alinhamento temporal (LLEVA) avançando os planetas até out/2024.
        5. Injeta os clones de partículas de teste (m=0) na época sincronizada.
        6. Define N_active = 9 (otimização de partículas sem massa).
        """
        if use_evorb_k is not None:
            self.use_evorb_k = use_evorb_k

        # Garante estado imutável/limpo mesmo se chamado repetidas vezes
        self._init_simulation()

        if planets is None and particles is None and "params" in params_or_datos:
            params = params_or_datos["params"]
            planets = params_or_datos["planets"]
            particles = params_or_datos["particles"]
        else:
            params = params_or_datos

        if not planets or not particles:
            raise ValueError("Planetas ou partículas ausentes para o setup da simulação.")

        # 1. Adicionar Sol central com raio físico (~696.340 km)
        cmass = params.get("cmass_msun", 1.0)
        rstar_km = params.get("rstar_km", 696340.0)
        self.sim.add(m=cmass, r=rstar_km / AU_IN_KM, name="Sol")

        # 2. Adicionar 8 planetas com seus elementos orbitais e raios físicos
        nomes_planetas = ["Mercurio", "Venus", "Terra", "Marte", "Jupiter", "Saturno", "Urano", "Netuno"]
        for i, p in enumerate(planets):
            nome = p.get("name", nomes_planetas[i] if i < len(nomes_planetas) else f"Planeta_{p['id']}")
            m = p.get("mass_msun", p.get("mass", 0.0))
            raio_km = p.get("radius_km", 0.0)
            self.sim.add(
                m=m,
                r=raio_km / AU_IN_KM if raio_km > 0 else 0.0,
                a=p["a"],
                e=p["e"],
                inc=np.radians(p["inc_deg"]),
                Omega=np.radians(p["node_deg"]),
                omega=np.radians(p["argper_deg"]),
                M=np.radians(p["M_deg"]),
                primary=self.sim.particles[0],
                name=nome
            )

        # Ajusta sistema para o Centro de Massa Inercial
        self.sim.move_to_com()

        # 3. Alinhamento temporal (LLEVA): avança planetas até out/2024
        delta_dias = params["t0_particles_jd"] - params["t0_planets_jd"]
        if delta_dias < 0:
            raise ValueError(f"Diferença temporal negativa entre épocas: {delta_dias} dias.")
        
        delta_anos = delta_dias / 365.25
        self.t_alinhamento = delta_anos * (2.0 * np.pi)
        self.sim.integrate(self.t_alinhamento)
        self.alinhado = True

        # 4. Injeção dos clones de teste na época sincronizada (m = 0)
        for part in particles:
            self.sim.add(
                m=0.0,
                a=part["a"],
                e=part["e"],
                inc=np.radians(part["inc_deg"]),
                Omega=np.radians(part["node_deg"]),
                omega=np.radians(part["argper_deg"]),
                M=np.radians(part["M_deg"]),
                primary=self.sim.particles[0],
                name=f"clone_{part['id']}"
            )

        # 5. Configurar partículas ativas (Sol + 8 planetas = 9)
        self.sim.N_active = 9

    def integrate_years(self, years: float):
        """
        Integra a simulação por um período adicional especificado em anos (Julianos).
        Inclui validação defensiva contra chamadas antes do setup.
        """
        if not self.alinhado or self.sim is None or self.sim.N == 0:
            raise RuntimeError("A simulação não foi configurada. Execute setup_from_datos() antes de integrar.")
        if years <= 0:
            raise ValueError(f"Período de integração deve ser positivo, recebido: {years}")

        target_t = self.sim.t + years * (2.0 * np.pi)
        self.sim.integrate(target_t)

    def energy(self) -> float:
        """Retorna a energia mecânica total do sistema ativo."""
        if self.sim is None or self.sim.N == 0:
            raise RuntimeError("Simulação não inicializada para cálculo de energia.")
        return float(self.sim.energy())

    def angular_momentum(self) -> np.ndarray:
        """
        Retorna o vetor de momento angular total (Lx, Ly, Lz) do sistema ativo.
        Definição: L = sum(m_i * (r_i x v_i))
        """
        if self.sim is None or self.sim.N == 0:
            raise RuntimeError("Simulação não inicializada para cálculo de momento angular.")

        L = np.zeros(3, dtype=np.float64)
        for i in range(self.sim.N_active):
            p = self.sim.particles[i]
            r = np.array([p.x, p.y, p.z], dtype=np.float64)
            v = np.array([p.vx, p.vy, p.vz], dtype=np.float64)
            L += p.m * np.cross(r, v)
        return L

    def get_body_state(self, index_or_name) -> Dict[str, float]:
        """Retorna posição e velocidade de um corpo."""
        if self.sim is None or self.sim.N == 0:
            raise RuntimeError("Simulação não inicializada.")
        p = self.sim.particles[index_or_name]
        return {
            "m": float(p.m),
            "x": float(p.x), "y": float(p.y), "z": float(p.z),
            "vx": float(p.vx), "vy": float(p.vy), "vz": float(p.vz)
        }

    def get_heliocentric_states(self) -> List[Dict[str, Any]]:
        """
        Retorna os vetores de estado de todos os 408 corpos (8 planetas + 400 clones)
        em coordenadas HELIOCÊNTRICAS e unidades idênticas ao EVORB (AU e AU/ano).
        
        Conversões aplicadas:
        - Posição: r_helio = r_corpo - r_sol
        - Velocidade: v_evorb = (v_corpo - v_sol) * 2pi
        """
        if self.sim is None or self.sim.N < 2:
            raise RuntimeError("Simulação requer no mínimo o Sol e mais corpos celestes.")

        sol = self.sim.particles[0]
        states = []
        for i in range(1, self.sim.N):
            p = self.sim.particles[i]
            # No REBOUND 5.x, acessa p.name diretamente
            nome = p.name if hasattr(p, "name") and p.name else f"body_{i}"
            states.append({
                "index": i - 1,  # 0 a 407 (espelhando os 408 corpos do EVORB)
                "name": nome,
                "m": float(p.m),
                "x": float(p.x - sol.x),
                "y": float(p.y - sol.y),
                "z": float(p.z - sol.z),
                "vx": float((p.vx - sol.vx) * (2.0 * np.pi)),
                "vy": float((p.vy - sol.vy) * (2.0 * np.pi)),
                "vz": float((p.vz - sol.vz) * (2.0 * np.pi)),
            })
        return states
