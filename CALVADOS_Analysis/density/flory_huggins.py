import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 26,
    "axes.titlesize": 26,
    "axes.labelsize": 26,
    "xtick.labelsize": 24,
    "ytick.labelsize": 24,
    "legend.fontsize": 24,
    "lines.linewidth": 3,
    "lines.markersize": 12
})

def analyze_flory_huggins(temp, rho_H, rho_L, rho_protein=900):
    """
    Analyze phase behavior using the Flory-Huggins theory.
    
    Parameters:
    - temp: list or array of temperatures (K)
    - rho_H: list or array of high density values (mg/cm^3)
    - rho_L: list or array of low density values (mg/cm^3)
    - rho_protein: protein density in mg/cm^3 (default: 900)

    The function plots:
    - Phase diagram (T vs φ)
    - χ vs 1/T
    - Free energy landscapes
    - Theoretical coexistence curve
    """

    temp = np.array(temp)
    rho_H = np.array(rho_H)
    rho_L = np.array(rho_L)

    phi_H = rho_H / rho_protein
    phi_L = rho_L / rho_protein

    # ---- A. Experimental Phase Diagram ----
    plt.figure(figsize=(8, 6))
    plt.scatter(phi_L, temp, label=r"$\phi'$", color='#5ecef1', marker='o')
    plt.scatter(phi_H, temp, label=r"$\phi''$", color='#ff7979', marker='s')
    plt.xlabel(r'$\phi$')
    plt.ylabel('T (K)')
    plt.tight_layout()
    plt.show()

    # ---- B. χ(T) Fit ----
    N1, N2 = 163, 1

    def calc_chi(phi_H, phi_L):
        num = (1/N1)*np.log(phi_H/phi_L) + (1/N2)*np.log((1-phi_L)/(1-phi_H))
        den = 2 * (phi_H - phi_L)
        return num / den

    chi_vals = calc_chi(phi_H, phi_L)

    def chi_model(T, A, B):
        return A + B / T

    chi_params, _ = curve_fit(chi_model, temp, chi_vals)
    A_fit, B_fit = chi_params
    T_fit = np.linspace(min(temp)-5, max(temp)+5, 300)
    chi_fit_vals = chi_model(T_fit, A_fit, B_fit)

    plt.figure(figsize=(8, 6))
    plt.plot(1/T_fit, chi_fit_vals, '-', color='orange', label='Fit: A + B/T')
    plt.scatter(1/temp, chi_vals, color='#456cee', label='Calculated χ')
    plt.xlabel(r'$1/T$ (K$^{-1}$)')
    plt.ylabel(r'$\chi$')
    plt.tight_layout()
    plt.show()

    # ---- C. Free Energy Landscape ----
    def F_mix(phi, chi, N1=163, N2=1):
        return (phi/N1)*np.log(phi) + ((1-phi)/N2)*np.log(1 - phi) + chi * phi * (1 - phi)

    T_list = np.arange(200, 320, 0.5)
    phi_range = np.linspace(0.001, 0.999, 500)

    plt.figure(figsize=(8, 6))
    cmap = plt.get_cmap('coolwarm')
    norm = plt.Normalize(vmin=T_list.min(), vmax=T_list.max())

    for T in T_list:
        chi = chi_model(T, A_fit, B_fit)
        F_vals = F_mix(phi_range, chi)

        if T == T_list.min():
            color = '#5ecef1'
        elif T == T_list.max():
            color = '#ff7979'
        else:
            color = cmap(norm(T))

        plt.plot(phi_range, F_vals, color=color, linewidth=0.5)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=plt.gca())
    cbar.set_label('Temperature (K)')

    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\Delta F_{\mathrm{mix}} $/ kT')
    plt.tight_layout()
    plt.show()

    # ---- D. Coexistence Curve ----
    coexisting_points = []

    for T in T_list:
        chi = chi_model(T, A_fit, B_fit)
        F_vals = F_mix(phi_range, chi)
        minima_indices = argrelextrema(F_vals, np.less)[0]

        if len(minima_indices) >= 2:
            phis = phi_range[minima_indices]
            phis_sorted = np.sort(phis)
            coexisting_points.append([T, phis_sorted[0], phis_sorted[-1]])

    coexisting_points = np.array(coexisting_points)

    plt.figure(figsize=(8, 6))
    plt.plot(coexisting_points[:,1], coexisting_points[:,0], color='blue', linewidth=2, label="Low φ")
    plt.plot(coexisting_points[:,2], coexisting_points[:,0], color='red', linewidth=2, label="High φ")
    plt.xlabel(r'$\phi$')
    plt.ylabel('T (K)')
    plt.title('Theoretical Coexistence Curve')
    plt.tight_layout()
    plt.show()
