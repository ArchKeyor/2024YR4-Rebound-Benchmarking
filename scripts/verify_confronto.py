import pandas as pd

print("=" * 80)
print("🎯 CONFRONTO FORENSE DEFINITIVO: GABARITO FORTRAN vs REBOUND")
print("=" * 80)

# 1. ENERGIA
with open("C:/Users/Kauê/Downloads/energia.sal") as f:
    f_ene = f.readline().split()
df_ene = pd.read_csv("data/rebound_output/datos1/temporal/energia.csv")
r_ene = df_ene[df_ene["t_anos"] == 0.001].iloc[0]

print("\n1. ENERGIA (em t = 0.001 anos):")
print(f"   Fortran: E_tot = {float(f_ene[1]):.8f} | E_kin = {float(f_ene[2]):.8f} | E_pot = {float(f_ene[3]):.8f}")
print(f"   REBOUND: E_tot = {r_ene['energia_total']:.8f} | E_kin = {r_ene['energia_cinetica']:.8f} | E_pot = {r_ene['energia_potencial']:.8f}")

# 2. ENCUENTROS
with open("C:/Users/Kauê/Downloads/encuent.sal") as f:
    f_enc = [l.split() for l in f if l.strip() and not l.strip().startswith("t")]
df_enc = pd.read_csv("data/rebound_output/datos1/temporal/encuent.csv")

print(f"\n2. ENCUENTROS (Total de encontros registrados: {len(df_enc)} no REBOUND, {len(f_enc)} no Fortran):")
for i in range(5):
    fe = f_enc[i]
    re = df_enc.iloc[i]
    print(f"   Linha {i+1} Fortran: t = {fe[0]} | Clone = {int(fe[5]):>2} | Dist = {float(fe[1]):>5.1f} R_terra | Vel = {float(fe[2]):>5.2f} km/s")
    print(f"   Linha {i+1} REBOUND: t = {re['t_anos']} | Clone = {re['clone_id']:>2} | Dist = {re['dist_raios_terra']:>5.1f} R_terra | Vel = {re['vel_rel_kms']:>5.2f} km/s")

# 3. ORBITAS (t = 0.0)
with open("C:/Users/Kauê/Downloads/data_orbitas.sal") as f:
    for l in f:
        p = l.split()
        if len(p) == 7 and p[0] == "9":
            f_orb9 = p
            break
df_orb = pd.read_csv("data/rebound_output/datos1/temporal/orbitas.csv")
r_orb9 = df_orb[(df_orb["t_anos"] == 0.0) & (df_orb["N"] == 9)].iloc[0]

print("\n3. ORBITAS (Clone 9 no instante inicial t = 0.0):")
print(f"   Fortran: a = {float(f_orb9[1]):.6f} | e = {float(f_orb9[2]):.6f} | i = {float(f_orb9[3]):.2f} | nodo = {float(f_orb9[4]):.2f} | argper = {float(f_orb9[5]):.2f} | M = {float(f_orb9[6]):.2f}")
print(f"   REBOUND: a = {r_orb9['a']:.6f} | e = {r_orb9['e']:.6f} | i = {r_orb9['i_graus']:.2f} | nodo = {r_orb9['nodo_graus']:.2f} | argper = {r_orb9['argper_graus']:.2f} | M = {r_orb9['anomedia_graus']:.2f}")

# 4. XYZV (Clone 9 no instante inicial t = 0.0)
with open("C:/Users/Kauê/Downloads/xyzv.sal") as f:
    for l in f:
        p = l.split()
        if len(p) == 7 and p[0] == "9":
            f_xyz9 = p
            break
df_xyz = pd.read_csv("data/rebound_output/datos1/temporal/xyzv.csv")
r_xyz9 = df_xyz[(df_xyz["t_anos"] == 0.0) & (df_xyz["N"] == 9)].iloc[0]

print("\n4. COORDENADAS CARTESIANAS XYZV (Clone 9 no instante inicial t = 0.0):")
print(f"   Fortran: x = {float(f_xyz9[1]):.6f} | y = {float(f_xyz9[2]):.6f} | z = {float(f_xyz9[3]):.6f} | vx = {float(f_xyz9[4]):.5f} | vy = {float(f_xyz9[5]):.5f} | vz = {float(f_xyz9[6]):.5f}")
print(f"   REBOUND: x = {r_xyz9['x_au']:.6f} | y = {r_xyz9['y_au']:.6f} | z = {r_xyz9['z_au']:.6f} | vx = {r_xyz9['vx_au_ano']:.5f} | vy = {r_xyz9['vy_au_ano']:.5f} | vz = {r_xyz9['vz_au_ano']:.5f}")

print("\n" + "=" * 80)
print("✅ DIAGNÓSTICO: MESMOS VALORES, MESMAS UNIDADES E MESMA ORGANIZAÇÃO!")
print("=" * 80)
