import numpy as np
import sympy as sp
import pandas as pd
import matplotlib.pyplot as plt

#Variables
h, m, c, D, tau, T_0, T_inf, T_tau = sp.symbols(r"h m c D \tau T_0 T_\infty T(\tau)", real = True)

h = -m*c/(sp.pi*D**2*tau)*sp.ln((T_tau-T_inf)/(T_0-T_inf))

vars_dep = [m, c, D, tau, T_0, T_inf, T_tau]
coefs = []
unidades = [r"\watt\kelvin\per\kg\per\m\squared", r"\kg\per\s\per\m\squared", r"\watt\kelvin\per\m\cubed", r"\watt\kelvin\per\s\per\m\squared", r"\watt\per\m\squared", r"\watt\per\m\squared", r"\watt\per\m\squared"]

for var in vars_dep:
    der = sp.Derivative(h, var)
    coefs.append(sp.simplify(sp.factor(der)))

with open("coeficientes_sensibilidad_h.tex", "w") as f:
    f.write(r"\begin{align}"+"\n")
    for var, coef, uni in zip(vars_dep, coefs, unidades):
        f.write(rf"c_{{{sp.latex(var)}}} &= \frac{{\partial h}}{{\partial {sp.latex(var)}}} = {sp.latex(coef)} \left[\unit{{{uni}}}\right]")
        if coef != coefs[-1]:
            f.write(r"\\" + "\n")
    f.write("\n" + r"\end{align}")
    f.close()