import numpy as np
import matplotlib.pyplot as plt

# LaTeX typeset
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": "Computer Modern Roman",
})
plt.rc('text', usetex=True)
plt.rc(
    'text.latex',
    preamble=r'\usepackage{amsfonts}\usepackage{amsmath}\usepackage{amssymb}\usepackage{siunitx}'
)

# ============================================================
# Parámetros Monte Carlo
# ============================================================

N = int(1e7)

# ============================================================
# Capacidad Calorífica
# ============================================================

c = 450 # J/(kg K)

# ============================================================
# Masa
# ============================================================

M = np.array([0.00688100, 0.00688098, 0.00688098, 0.00688111, 0.00688108, 0.00688103])  # kg

km = 2

mprom = np.mean(M)
mres = 0.00000001   # kg
umres = mres / np.sqrt(12)

umrep = np.std(M, ddof=1) / np.sqrt(len(M))

mcert = 0.000002 # kg
umcert = mcert / km

# ============================================================
# Diámetro
# ============================================================

D = np.array([0.01187, 0.01187, 0.01186, 0.01187, 0.01187, 0.01186])  # m

kd = 2

dprom = np.mean(D)

dres = 0.00001   # m
udres = dres / np.sqrt(12)

udrep = np.std(D, ddof=1) / np.sqrt(len(D))

dcert = 0.00002  # m
udcert = dcert / kd

# ============================================================
# Temperatura ambiente
# ============================================================

Tamb = np.array([28.1, 28.5, 28.7, 28.9, 28.6, 28.8, 29.2, 28.6, 28.3, 28.3])  # °C

ktamb = 2

tambprom = np.mean(Tamb)

tambres = 0.1   # °C
utambres = tambres / np.sqrt(12)

utambrep = np.std(Tamb, ddof=1) / np.sqrt(len(Tamb))

tambcert = 0.6  # °C
utambcert = tambcert / ktamb

# ============================================================
# Tiempo
# ============================================================

Tiempo = np.array([270, 270, 270, 270, 270, 270, 270, 270, 270, 270]) # s

tiempoprom = np.mean(Tiempo)

tiempores = 1/60 # s
utiempores = tiempores / np.sqrt(12)

utiemporep = np.std(Tiempo, ddof=1) / np.sqrt(len(Tiempo))

# ============================================================
# Coeficiente de enfriamiento
# ============================================================

k = np.array([0.00510096620748971, 0.00502045445353014, 0.00480975341689432, 0.00525998486220582, 0.00531372787002526, 0.00520215735227422, 0.00579315943666277, 0.00517885584225008, 0.00547634818039660, 0.00472998173761074])  # s^-1, van con toda la precisión que tenía Excel

kprom = np.mean(k)

ukrep = np.std(k, ddof=1) / np.sqrt(len(k))

# ============================================================
# Monte Carlo vectorizado
# ============================================================

# Masa
M_mc = (mprom + np.random.uniform(-mres/2, mres/2, N) + np.random.normal(0, umrep, N) + np.random.normal(0, umcert, N))

# Temp Ambiente
Tamb_mc = (tambprom + np.random.uniform(-tambres/2, tambres/2, N) + np.random.normal(0, utambrep, N) + np.random.normal(0, utambcert, N))

# Diámetro
D_mc = (dprom + np.random.uniform(-dres/2, dres/2, N) + np.random.normal(0, udrep, N) + np.random.normal(0, udcert, N))

# Tiempo
Tiempo_mc= (tiempoprom + np.random.uniform(-tiempores/2, tiempores/2, N) + np.random.normal(0, utiemporep, N))

# Coeficiente enfriamiento
k_mc = (kprom + np.random.normal(0, ukrep, N))

# Lista de valores de coeficiente conveccion
h_lista = M_mc*c/(np.pi*D_mc**2)*k_mc

# ============================================================
# Estadísticos
# ============================================================

hprom = np.mean(h_lista)

u_h = np.std(h_lista, ddof=1)

kh = 2
U_h = kh * u_h

h025 = np.percentile(h_lista, 2.5)
h975 = np.percentile(h_lista, 97.5)

# ============================================================
# Resultados
# ============================================================

print("="*50)
print("Resultado Monte Carlo:")
print("-"*20)
print(f"h = {hprom:.6f} W/(K m\u00B2)")
print(f"u(R) = {u_h:.6g} W/(K m\u00B2)")
print(f"U(R) = {U_h:.6g} W/(K m\u00B2)  (k = 2)")

print(f"Intervalo central 95 % (con percentiles):")
print(f"[{h025:.6f}, {h975:.6f}] W/(K m\u00B2)")

print(f"Intervalo central 95 % (k*u(R)):")
print(f"[{hprom-kh*u_h:.6f}, {hprom+kh*u_h:.6f}] W/(K m\u00B2)")

print("="*50)

# ============================================================
# Histograma
# ============================================================

plt.figure(figsize=(8, 5))

plt.hist(h_lista, bins=100, color="#7aaa97", ec="#05242f", lw=0.1)

plt.axvline(hprom, c="#05242f", linestyle="-.", label=r"$h_{\mathrm{prom}}$")
plt.axvline(h025, c="#05242f", linestyle=":", label=r"Percentil \qty{2.5}{\percent}")
plt.axvline(h975, c="#05242f", linestyle=":")
plt.axvline(hprom-kh*u_h, c="#05242f", linestyle="--", label=r"$\pm 2 \sigma$")
plt.axvline(hprom+kh*u_h, c="#05242f", linestyle="--")

plt.title(r"Distribución de $h$")
plt.xlabel(r"$h$ [\unit{\watt\kelvin\per\meter\squared}]")
plt.ylabel("Cantidad")
plt.grid(True, alpha=0.4, linestyle="--")
plt.legend()

plt.tight_layout()
plt.savefig(f"h.pdf", dpi=300, bbox_inches='tight', pad_inches=0)
plt.show()