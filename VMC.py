import numpy as np

# ==========================================
# BAGIAN 1: Kode dari eloc.py (Local Energy)
# ==========================================

def psi(al, oms, X):
    """
    Calculate trial wave function.
    """
    rexp = np.sum(oms * X**2)
    return np.exp(-0.5 * al * rexp)

def ekin(al, oms, X, h=0.01, hom=1.):
    """
    Calculate kinetic energy using finite difference.
    """
    npart, ndim = X.shape
    psiold = psi(al, oms, X)
    kin = 0.
    
    for (j, el), r in np.ndenumerate(X):
        X[j, el] = r + h
        psip = psi(al, oms, X)
        X[j, el] = r - h
        psim = psi(al, oms, X)
        X[j, el] = r 
        lapl = (psip + psim - 2. * psiold) / h**2
        kin += -0.5 * hom * lapl
        
    return kin / psiold

def epot(oms, X, g=3.):
    """
    Calculate potential energy.
    """
    npart, ndim = X.shape
    
    v1 = 0.5 * np.sum(oms**2 * X**2)
    v2 = 0.
    for k in range(1, npart):
        for j in range(k):
            r = np.sqrt(np.sum((X[j, :] - X[k, :])**2))
            v2 += np.exp(-r**2)
    return v1 + g * v2

# ==============================================
# BAGIAN 2: Kode dari vmc.py (Metropolis Algorithm)
# ==============================================

def vmc(n_particles, nd, al, oms, seed=314159):
    """
    Perform Variational Monte Carlo simulation.
    """
    y = 10**4    
    nm = 10      
    th = 0.5     
    
    np.random.seed(seed)
    
    X = np.random.uniform(-1., 1., (n_particles, nd))
    
    elocs = []
    
    # Loop Metropolis
    for i in range(y):
        Xnew = X + np.random.uniform(-th, th, (n_particles, nd))        
        ratio = (psi(al, oms, Xnew) / psi(al, oms, X))**2
        if ratio > np.random.uniform():
            X = Xnew
        if i % nm == 0:
            el = ekin(al, oms, X) + epot(oms, X)
            elocs.append(el)
    return np.mean(elocs), np.std(elocs) / np.sqrt(len(elocs))

# ==========================================
# BAGIAN 3: Eksekusi Utama (Main Block)
# ==========================================

import matplotlib.pyplot as plt

def run_simulation(alpha, n_particles, n_dimensions, omegas, return_trace=False):
    """
    Runs VMC simulation for a given alpha.
    Returns (mean_energy, std_error) or (mean_energy, std_error, energy_trace)
    """
    n_samples = 2000 
    nm = 100
    n_steps = n_samples * nm
    th = 0.8
    
    np.random.seed(42) 
    X = np.random.uniform(-1., 1., (n_particles, n_dimensions))
    elocs = []
    
    for i in range(n_steps):
        Xnew = X + np.random.uniform(-th, th, (n_particles, n_dimensions))
        ratio = (psi(alpha, omegas, Xnew) / psi(alpha, omegas, X))**2
        
        if ratio > np.random.uniform():
            X = Xnew
            
        if i % nm == 0:
            el = ekin(alpha, omegas, X) + epot(omegas, X)
            elocs.append(el)
            
    mean_E = np.mean(elocs)
    std_err = np.std(elocs) / np.sqrt(len(elocs))
    
    if return_trace:
        return mean_E, std_err, elocs
    return mean_E, std_err

def golden_section_search(f, a, b, tol=1e-4):
    """
    Finds the minimum of function f in the interval [a, b] 
    using the Golden Section Search algorithm.
    """
    gr = (np.sqrt(5) + 1) / 2
    
    c = b - (b - a) / gr
    d = a + (b - a) / gr
    
    while abs(b - a) > tol:
        if f(c) < f(d):
            b = d
        else:
            a = c
        c = b - (b - a) / gr
        d = a + (b - a) / gr
        
    return (b + a) / 2, f((b + a) / 2)

def run_optimization():
    """
    1. Scans alpha for initial plot and logging.
    2. Uses custom Golden Section Search for precise optimum.
    """
    n_particles = 4
    n_dimensions = 3
    omegas = np.arange(1, n_dimensions + 1, dtype=float)
    
    print("Melakukan scanning alpha...")
    alphas = np.arange(0.5, 1.15, 0.05)
    energies = []
    errors = []
    
    for al in alphas:
        mean_E, std_err = run_simulation(al, n_particles, n_dimensions, omegas)
        energies.append(mean_E)
        errors.append(std_err)
        print(f"Alpha: {al:.2f}, Energi: {mean_E:.5f} +/- {std_err:.5f}")
        
    # Plot Energy vs Alpha
    plt.figure(figsize=(8, 5))
    plt.errorbar(alphas, energies, yerr=errors, fmt='o-', capsize=3, label='VMC Energy')
    plt.xlabel(r'Variational Parameter ($\alpha$)',  fontsize=16)
    plt.ylabel(r'Energy ($\hbar\omega$)', fontsize=16)
    # plt.title('VMC Optimization: Energy vs Alpha', fontsize=18)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig('energy_vs_alpha.pdf')
    print("Grafik disimpan: energy_vs_alpha.pdf")
    
    
    # --- Part 2: Global Optimization ---
    print("\nMencari global minimum presisi untuk alpha (Golden Section Search)...")
    
    def objective(alpha):
        e, _ = run_simulation(alpha, n_particles, n_dimensions, omegas)
        return e

    opt_alpha, min_energy = golden_section_search(objective, 0.5, 1.1, tol=1e-3)
    
    print(f"\nOptimal Alpha: {opt_alpha:.5f}")
    print(f"Minimum Energy: {min_energy:.5f}")
    
    print(f"Generating energy trace for optimal alpha...")
    target_alpha = opt_alpha
    _, _, energy_trace = run_simulation(target_alpha, n_particles, n_dimensions, omegas, return_trace=True)
            
    plt.figure(figsize=(8, 5))
    plt.plot(energy_trace, 'k.', markersize=2, alpha=0.5)
    plt.title(f'Local Energy Trace (Alpha={target_alpha:.4f})', fontsize=18)
    plt.xlabel('Sample Block', fontsize=16)
    plt.ylabel(r'Energy ($\hbar\omega$)', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.savefig('energy_trace.pdf')
    print("Grafik disimpan: energy_trace.pdf")

if __name__ == '__main__':
    run_optimization()