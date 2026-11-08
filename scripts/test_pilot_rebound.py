"""Pilot simulation script for REBOUND Phase 1."""

import time
from pathlib import Path
import sys
import numpy as np

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import rebound
from src.parsers.evorb_parser import parse_datos


def run_pilot():
    print("=" * 65)
    print("🚀 FASE 1: PROVA DE CONCEITO - SIMULAÇÃO PILOTO NO REBOUND (IAS15)")
    print("=" * 65)

    datos_path = PROJECT_ROOT / "data" / "evorb_raw" / "datos1" / "datos.dat"
    if not datos_path.exists():
        print(f"[-] Erro: Arquivo {datos_path} não encontrado!")
        return

    # 1. Carregar parâmetros do datos.dat
    print(f"[1/5] Carregando condições do arquivo: {datos_path.name}...")
    data = parse_datos(datos_path)
    params = data["params"]
    planets = data["planets"]
    particles = data["particles"]

    print(f"      -> {len(planets)} planetas principais decodificados.")
    print(f"      -> {len(particles)} clones do asteroide 2024 YR4 decodificados.")

    # 2. Configurar Simulação REBOUND
    print("\n[2/5] Configurando simulador REBOUND com motor IAS15...")
    sim = rebound.Simulation()
    sim.units = ("AU", "Msun", "yr2pi")
    sim.integrator = "ias15"

    # Adicionar Sol
    sim.add(m=params["cmass_msun"], name="Sol")

    # Adicionar os 8 Planetas em J2000
    nomes_planetas = ["Mercurio", "Venus", "Terra", "Marte", "Jupiter", "Saturno", "Urano", "Netuno"]
    for p, nome in zip(planets, nomes_planetas):
        sim.add(
            m=p["mass_msun"],
            a=p["a"],
            e=p["e"],
            inc=np.radians(p["inc_deg"]),
            Omega=np.radians(p["node_deg"]),
            omega=np.radians(p["argper_deg"]),
            M=np.radians(p["M_deg"]),
            primary=sim.particles[0],
            name=nome
        )
    sim.move_to_com()
    print(f"      -> Sol e 8 planetas posicionados na época J2000 (JD {params['t0_planets_jd']:.1f}).")

    # 3. Alinhamento Temporal de Épocas (Sub-rotina LLEVA)
    delta_dias = params["t0_particles_jd"] - params["t0_planets_jd"]
    delta_anos = delta_dias / 365.25
    t_alinhamento = delta_anos * 2 * np.pi

    print(f"\n[3/5] Executando alinhamento de épocas (LLEVA):")
    print(f"      -> Avançando planetas em {delta_dias:.1f} dias ({delta_anos:.4f} anos) até out/2024...")
    t0_start = time.perf_counter()
    sim.integrate(t_alinhamento)
    dt_alinhamento = time.perf_counter() - t0_start
    print(f"      -> Alinhamento concluído em {dt_alinhamento:.3f} segundos.")

    # 4. Injeção dos 400 Clones na Época Sincronizada
    print(f"\n[4/5] Injetando {len(particles)} clones do asteroide como partículas sem massa (m = 0)...")
    for p in particles:
        sim.add(
            m=0.0,
            a=p["a"],
            e=p["e"],
            inc=np.radians(p["inc_deg"]),
            Omega=np.radians(p["node_deg"]),
            omega=np.radians(p["argper_deg"]),
            M=np.radians(p["M_deg"]),
            primary=sim.particles[0],
            name=f"clone_{p['id']}"
        )

    print(f"      -> Sistema completo montado com {sim.N} corpos celestes.")
    sim.N_active = 9  # Sol + 8 planetas ativos (clones são partículas de teste m=0)
    E_ini = sim.energy()

    # 5. Execução do Piloto de 1 Ano
    t_piloto = t_alinhamento + 1.0 * (2 * np.pi)
    print(f"\n[5/5] Executando simulação piloto de 1 ano completo...")
    t1_start = time.perf_counter()
    sim.integrate(t_piloto)
    dt_piloto = time.perf_counter() - t1_start
    E_fim = sim.energy()

    dE = abs((E_fim - E_ini) / E_ini)

    # Resultados e Diagnóstico
    terra = sim.particles["Terra"]
    clone9 = sim.particles["clone_9"]
    r_terra = np.sqrt(terra.x**2 + terra.y**2 + terra.z**2)
    r_clone9 = np.sqrt(clone9.x**2 + clone9.y**2 + clone9.z**2)
    dist_terra_clone9 = np.sqrt((clone9.x - terra.x)**2 + (clone9.y - terra.y)**2 + (clone9.z - terra.z)**2)

    print("\n" + "=" * 65)
    print("📊 RESULTADOS E DIAGNÓSTICO DO PILOTO (FASE 1)")
    print("=" * 65)
    print(f"⚡ Tempo de execução para 1 ano (409 corpos): {dt_piloto:.3f} segundos")
    print(f"🔋 Conservação de Energia (dE/E):             {dE:.4e}")
    print(f"🌍 Distância Terra - Sol:                     {r_terra:.4f} AU")
    print(f"☄️  Distância Asteroide (Clone 9) - Sol:        {r_clone9:.4f} AU")
    print(f"📍 Distância Terra - Asteroide (Clone 9):     {dist_terra_clone9:.4f} AU ({dist_terra_clone9 * 149597870.7 / 1e6:.2f} milhões de km)")

    # Critérios de Aceite
    print("\n📋 VERIFICAÇÃO DOS CRITÉRIOS DE ACEITE:")
    status_energia = "✅ APROVADO" if dE < 1e-10 else "❌ REPROVADO"
    status_tempo = "✅ APROVADO" if dt_piloto < 3.0 else "⚠️ ALERTA"
    status_corpos = "✅ APROVADO" if sim.N == 409 else "❌ REPROVADO"

    print(f"   [1] Conservação de Energia (< 1e-10):   {status_energia} ({dE:.2e})")
    print(f"   [2] Tempo de Execução (< 3.0s):         {status_tempo} ({dt_piloto:.3f}s)")
    print(f"   [3] Contagem de Corpos (409 ativos):    {status_corpos} ({sim.N} corpos)")

    if dE < 1e-10 and dt_piloto < 3.0 and sim.N == 409:
        print("\n🎉 FASE 1 CONCLUÍDA COM SUCESSO! O motor está perfeitamente calibrado.")


if __name__ == "__main__":
    run_pilot()
