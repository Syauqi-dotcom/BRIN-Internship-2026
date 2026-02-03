# H2 Molecule Ground State Energy Simulation
### Classical Variational Monte Carlo (VMC) vs. Quantum Variational Quantum Eigensolver (VQE)

This repository contains a comparative study of classical and quantum algorithms for calculating the ground state energy of the Hydrogen molecule ($H_2$) at a bond distance of $0.735 \mathring{A}$.

## Project Overview

The goal of this project is to solve the time-independent Schrödinger equation for the $H_2$ molecule using variational methods:
1.  **Classical Approach**: A Variational Monte Carlo (VMC) simulation using a trial wavefunction and stochastic integration.
2.  **Quantum Approach**: A Variational Quantum Eigensolver (VQE) utilizing Qiskit's V2 Primitives to simulate a quantum circuit ansatz.

## Files Description

### 1. `H2Classic.ipynb` (Classical VMC)
This notebook implements the Variational Monte Carlo method from scratch using `numpy`.

*   **Physics**:
    *   Uses a trial wavefunction $\Psi_T(\vec{r}, \alpha)$ parametrized by a variational parameter $\alpha$.
    *   Calculates the **Local Energy** $E_L(\vec{r}) = \frac{\hat{H}\Psi_T}{\Psi_T}$ at sampled walker positions.
    *   Approximates the ground state energy integral $\langle E \rangle = \int P(\vec{r}) E_L(\vec{r}) d\vec{r}$ using Metropolis-Hastings sampling.

*   **Key Features**:
    *   **Custom Hamiltonian**: Implements the full real-space Hamiltonian for $H_2$ (kinetic + electron-Proton + electron-electron terms).
    *   **Optimization Strategies**:
        *   **Grid Search**: Manually scanning $\alpha$ values to find the global minimum.
        *   **Gradient Descent**: An on-the-fly optimization loop that updates $\alpha$ based on energy gradients.
    *   **Visualization**: Plots of Energy vs. Iteration and Energy vs. Alpha.

### 2. `H2VQE.ipynb` (Quantum VQE)
This notebook leverages **Qiskit** to solve the same problem using a quantum computer simulator.

*   **Physics**:
    *   **Hamiltonian**: Uses a mapped qubit Hamiltonian (Parity/Jordan-Wigner) reduced to a 1-qubit operator for the active space.
    *   **Ansatz**: A parameterized quantum circuit $U(\theta)$ designed to prepare the ground state.
        *   Structure: $R_x(\theta_0) R_z(\theta_1) R_x(\theta_2)$ (3 parameters).
    *   **Estimator**: Uses Qiskit's `StatevectorEstimator` (Primitives V2) for exact expectation value calculation (infinite shots limit).

*   **Key Features**:
    *   **Optimizer**: Uses the classical **COBYLA** optimizer via `scipy.optimize.minimize`.
    *   **Custom Logging**: A real-time callback logger (`AnsatzLogger` / `VQELogger`) to track parameters and energy convergence at every iteration.
    *   **Landscape Analysis**: Visualizes the optimization path and parameter convergence.

## Results Comparison

| Method | Description | Result (Hartree) |
| :--- | :--- | :--- |
| **Exact** | Full CI / Physical Constant | **-1.1744** |
| **VQE (Quantum)** | Qiskit Statevector (3-param ansatz) | **-1.1460** |
| **VMC (Classical)** | Simple Exponential Trial Wavefunction | **~ -1.1280** |

*   **VQE Performance**: The quantum method achieves a result very close to the exact value, limited primarily by the expressibility of the chosen ansatz and the reduction of the Hamiltonian.
*   **VMC Performance**: The classical method produces a reasonable upper bound but is limited by the lack of explicit electron-electron correlation terms (Jastrow factor) in the simple trial wavefunction used.

## Requirements

*   Python 3.8+
*   `numpy`
*   `scipy`
*   `matplotlib`
*   `qiskit>=1.0.0`

## Usage

Run the notebooks using Jupyter Lab or Jupyter Notebook:

```bash
jupyter notebook H2VQE.ipynb
```
or
```bash
jupyter notebook H2Classic.ipynb
```
