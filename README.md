# 🔬 2D_NNF_LBM: Lattice Boltzmann Method for Non-Newtonian Fluid-Particle Interactions

>  *Advanced CFD solver for complex fluid-structure interaction problems*
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
---

## 👨‍🔬 Author Information

**Taha Rezaee**  
📧 **Email:** [rezaee.taha@gmail.com](mailto:rezaee.taha@gmail.com)  |  🔬 **ORCID:** [![ORCID](https://img.shields.io/badge/ORCID-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0008-6807-0765)

---

## 📖 About

This repository contains a Fortran code implementing a Lattice Boltzmann Method (LBM) solver, developed during my PhD research. The solver is designed for complex multi-physics simulations, including:

*   **Fluid-Solid Interaction (FSI):** Simulating the motion of rigid and deformable particles in Newtonian and non-Newtonian fluids.
*   **Non-Newtonian Fluids:** Modeling viscoplastic (Bingham-Casson) fluids with yield stress.
*   **Porous Media:** Simulating the dynamics of porous particles within a fluid domain.

The code has been rigorously validated against established benchmarks and finite element method (FEM) results, demonstrating excellent accuracy and robustness. It has been used as the core simulation engine for several peer-reviewed publications.

## 🎯 Featured Publications

This code has been used in the following scientific publications:

1.  **Rezaee, T.** (2025). Why Do Viscoplastic Fluids Adhere to Moving Surfaces? *Physics of Fluids*. [![DOI](https://img.shields.io/badge/DOI-10.1063%2F5.0265952-blue)](https://doi.org/10.1063/5.0265952)
2.  **Rezaee, T.** (2025). Dynamics of a neutrally-buoyant elliptic particle in a square cavity with two adjacent moving walls: a lattice Boltzmann simulation. *European Journal of Mechanics B/Fluids*. [![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.euromechflu.2025.20426-blue)](https://doi.org/10.1016/j.euromechflu.2025.20426)
3.  **Rezaee, T.**, et al. (2024). On the use of sinusoidal vibrations for disaggregating clusters of non-settling inertial particles immersed in yield-stress fluids. *Journal of Non-Newtonian Fluid Mechanics*.  [![DOI](https://img.shields.io/badge/DOI-10.1016/j.jnnfm.2024.105261-blue)]([https://doi.org/10.1016/j.euromechflu.2025.20426](https://doi.org/10.1016/j.jnnfm.2024.105261)) 
4.  **Rezaee, T.**, et al. (2023). Settling dynamics of circular particles in vibrating tanks filled with a yield-stress liquid. *Physics of Fluids*, 35(5). [![DOI](https://img.shields.io/badge/DOI-10.1063/5.0150359-blue)](https://doi.org/10.1063/5.0150359)
5.  **Rezaee, T.**, et al. (2019). Effect of porosity on the settling behavior of a 2D elliptic particle in a narrow vessel: a lattice Boltzmann simulation. *Physics of Fluids*, 31(12). [![DOI](https://img.shields.io/badge/DOI-10.1063/1.5125172-blue)](https://doi.org/10.1063/1.5125172)

---

## 🧪 Validation Cases

The solver's capabilities have been verified through the following benchmark problems. Each case has a dedicated folder containing the Fortran source code and necessary input files to reproduce the results.

### 1. Free Fall of a Rigid Circular Particle in a Newtonian Fluid
*   **Folder:** `1_FreeFall_Circular`
*   **Description:** A rigid circular particle is released from rest in a closed 2D channel filled with a Newtonian fluid and settles under gravity. The simulation captures the particle's trajectory, velocity, translational kinetic energy, and Reynolds number until it impacts the bottom wall.
*   **Validation:** Results for vertical position, velocity, and kinetic energy show excellent agreement with Finite Element Method (FEM) results from **Wan & Turek**.
*   **Key Highlight:** The present LBM algorithm correctly predicts a zero final velocity upon impact, accurately modeling the particle as rigid, unlike some FEM approaches that show residual velocity due to mesh deformation.

<p align="center">
  <img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/sehe1.jpg" alt="Schematic 1" width="10%">
  <br>
  <em>Figure 1.1: Geometry of a circular particle settling in a Newtonian fluid.</em>
</p>

<p align="center">
  <img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/y.jpeg" alt="Vertical Position" width="40%">
  <img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/uy.jpeg" alt="Vertical Velocity" width="40%">
  <br>
  <em>Figure 1.2: Time history of vertical position (left) and vertical velocity (right).</em>
</p>

---

### 2. Interaction of Two Settling Particles in a Newtonian Fluid
*   **Folder:** `2_DKT_Newtonian`
*   **Description:** This case tests the implemented force scheme by simulating the interaction and drafting-kissing-tumbling (DKT) dynamics of two settling circular particles in a Newtonian fluid.
*   **Validation:** The time histories of vertical/horizontal position and velocity for both particles match well with the FEM results of **Wan & Turek**.
*   **Key Highlight:** The LBM solver successfully captures the complex interaction where the trailing particle (P1) catches up and collides with the leading particle (P2) due to reduced drag in its wake.

<p align="center">
  <img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/sehe2.jpg" alt="Schematic 2" width="15%">
  <br>
  <em>Figure 2.1: Geometry of two circular particles settling in a Newtonian fluid.</em>
</p>

<p align="center">
  <img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/y_dkt.jpeg" alt="Vertical Pos 2P" width="40%">
  <img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/uy_dkt.jpeg" alt="Horizontal Pos 2P" width="40%">
  <br>
  <em>Figure 2.2: Time history of vertical position (left) and vertical velocity (right) for the two interacting particles. (Solid lines: Present LBM, Symbols: Wan & Turek).</em>
</p>

---

### 3. Sedimentation of Multiple Particles in a 2D Closed Cavity
*   **Folder:** `3_MultiParticle_Sedimentation`
*   **Description:** This case evaluates the solver's ability to handle a large number of particles. 504 circular particles are released from rest in a closed 2D cavity.
*   **Validation:** The overall sedimentation process and the formation of two vortices near the side walls are consistent with the finite difference results of **Glowinski et al.**.
*   **Key Highlight:** The solver efficiently manages multi-body interactions and collision dynamics in a dense system.

<p align="center">
  <table align="center">
    <tr>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/n0.jpg" alt="Image 1" width="90%"></td>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/n1.png" alt="Image 2" width="90%"></td>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/n2.jpg" alt="Image 3" width="90%"></td>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/n3.jpg" alt="Image 4" width="90%"></td>
    </tr>
    <tr>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/n4.jpg" alt="Image 5" width="90%"></td>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/n8.jpg" alt="Image 6" width="90%"></td>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/n12.jpg" alt="Image 7" width="90%"></td>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/n24.jpg" alt="Image 8" width="90%"></td>
    </tr>
  </table>
  <br>
</p>
<p align="center">
  <em>Figure 3.1: Instantaneous positions of particles at different times during multi-particle sedimentation.</em>
</p>

---

### 4. Oscillating Circular Particle in a Stationary Fluid
*   **Folder:** `4_Oscillating_Particle`
*   **Description:** A circular particle oscillates horizontally with a prescribed velocity in a stationary Newtonian fluid at a high Reynolds number (Re=100) and Keulegan-Carpenter number (KC=5).
*   **Validation:** The computed drag coefficient history and velocity profiles in the domain show excellent agreement with numerical and experimental data from **Dütsch et al. [75]**.
*   **Key Highlight:** This case underscores the critical importance of including the **internal mass (or added mass)** effect in the momentum exchange calculation for unsteady flows. Neglecting it leads to significant errors in the predicted force on the body.

<p align="center">
  <img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/sehe5.jpg" alt="Schematic Oscillating" width="30%">
  <img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/Cd.jpeg" alt="Drag Coefficient Osc" width="30%">
  <br>
  <em>Figure 4.1: (left) Geometry of an oscillating particle in a stationary fluid. (right) Time history of the drag coefficient (Re=100, KC=5). </em>
</p>

<p align="center">
  <table align="center">
    <tr>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/ux_a.jpeg" alt="Image 1" width="90%"></td>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/ux_b.jpeg" alt="Image 2" width="90%"></td>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/ux_c.jpeg" alt="Image 3" width="90%"></td>
    </tr>
    <tr>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/uy_a.jpeg" alt="Image 5" width="90%"></td>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/uy_b.jpeg" alt="Image 6" width="90%"></td>
      <td align="center"><img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/uy_c.jpeg" alt="Image 7" width="90%"></td>
    </tr>
  </table>
  <br>
  <em>Figure 4.2: Comparison of the velocity components at four cross-sections with constant x-value. Phase position (left) 180◦, (middle) 210◦ and (right) 330◦.</em>
</p>

---

### 5. Drag on a Particle in Creeping Flow of a Yield-Stress Fluid
*   **Folder:** `5_YieldStress_Drag`
*   **Description:** The drag force on a stationary circular particle in a creeping flow of a Bingham-Papanastasiou fluid is computed for a wide range of Bingham numbers (Bn) and channel confinement ratios (H/D).
*   **Validation:** The computed drag coefficients show remarkable agreement with Finite Element Method results from **Mitsoulis**.
*   **Key Highlight:** The solver reliably computes forces in yield-stress fluids. Results show that drag becomes independent of channel confinement for Bn ≥ 10 and increases sharply with Bn.

<p align="center">
  <img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/sehe4.jpg" alt="Schematic Yield Stress Drag" width="30%">
  <br>
  <em>Figure 5.1: Geometry for drag on a particle in creeping flow of a yield-stress fluid.</em>
</p>

<p align="center">
  <img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/bcc2.jpeg" alt="Drag Coefficient Yield" width="40%">
  <br>
  <em>Figure 5.2: Drag coefficient on a particle in Bingham plastic creeping flow. (Lines: Present LBM, Symbols: FEM).</em>
</p>

---

### 6. Free Fall of a Rigid Particle in a Yield-Stress Fluid
*   **Folder:** `6_FreeFall_YieldStress`
*   **Description:** A rigid circular particle settles under gravity in a 2D channel filled with a Bingham fluid. The terminal velocity and settling behavior are investigated for different yield numbers (Y).
*   **Validation:** The time history of the settling velocity for a specific yield number (Y=0.033) agrees well with the results of **Wachs & Frigaard**, who used an Augmented Lagrangian (AL) with Distributed Lagrange Multiplier/Fictitious Domain (DLM/FD) method.
*   **Key Highlight:** Confirms the solver's capability to handle fully coupled fluid-particle dynamics in complex non-Newtonian fluids.

<p align="center">
  <img src="https://github.com/taha-re/2D_NNF_LBM/blob/main/media/uy_bingham.jpeg" alt="Settling Velocity Yield" width="40%">
  <br>
  <em>Figure 6.1: Time history of settling velocity in a yield-stress fluid for Y=0.033.</em>
</p>

---

## 🚀 Getting Started

1.  **Navigate** to the folder of the validation case you are interested in.
2.  **Compile** the Fortran source file using a standard Fortran compiler (e.g., `gfortran`).
    ```bash
    gfortran -O3 -o lbm_solver main.f90
    ```
3.  **Run** the executable.
    ```bash
    ./lbm_solver
    ```
4.  The code will generate output files (which should be post-processesed by python or Tecplot) in the folder.

## 📄 License

This project is licensed under the MIT License.

