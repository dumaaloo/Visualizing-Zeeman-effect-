# Breit-Rabi Simulation: Spin-Orbit Coupling & The Zeeman Effect

This repository contains a computational physics model written in **Scilab** that simulates the transition of an atomic $p$-orbital ($l=1$) electron from the weak-field regime (dominated by internal spin-orbit coupling) to the strong-field regime (the Paschen-Back effect).

## Physical Overview
The simulation builds and diagonalizes a $6 \times 6$ Hamiltonian matrix accounting for:
1. **Spin-Orbit Interaction ($\vec{L} \cdot \vec{S}$):** Handled via $L_zS_z$ on the diagonal and raising/lowering ladder operators  $\frac{1}{2}(L_+S_- + L_-S_+)$  on the off-diagonals.
2. **External Magnetic Field ($B_z$):** Applied along the quantization axis to interact linearly with the magnetic moments.

## Key Insights & Results
By analyzing the eigenvectors at different magnetic field strengths ($B_z$), the model successfully tracks the uncoupling of the quantum states:
* **Weak Field ($B_z = 0$):** States are highly mixed, naturally yielding the exact analytic **Clebsch-Gordan coefficients** (e.g., mixtures of $\sqrt{1/3}$ and $\sqrt{2/3}$).
* **Strong Field ($B_z \to \infty$):** The external field dominates. Unique states completely uncouple into pure basis vectors (1s and 0s), while degenerate states retain their geometric combinations due to matrix symmetry.

## How to Run
1. Open the `.sce` script in Scilab.
2. Run the script to generate the eigenvalue evolution plots and output the state transition matrices to the console.
