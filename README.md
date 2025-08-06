# 3D Finite Element Simulation of Heterogeneous Domains in MATLAB

This repository provides a customizable MATLAB implementation for simulating 3D solid structures using the Finite Element Method (FEM). It supports domains composed of elements with different materials and sizes, enabling the modeling of heterogeneous media.

---

## Motivation

This project was developed to study and prototype advanced structures composed of various materials, using FEM with 8-node isoparametric hexahedral elements. It is particularly suited for modeling structures with spatially varying properties — including materials with frequency-dependent or anisotropic behavior — and offers flexibility to simulate periodic and resonant systems.

While a separate repository contains the theoretical formulation of the finite elements used ([see here](https://github.com/Arincon94/Assembly-of-Elemental-Matrices-in-the-FEM)), this implementation focuses on the **assembly of complex domains**, allowing regions to be defined with distinct materials, including the possibility of voids (zero-stiffness regions).

---

## Test Example: Unit Cell of a Sandwich Metastructure

As a demonstrative case, this repository simulates a **unit cell of a sandwich-type metastructure**, composed of:

- **PolyJet DM 8510**: used in the walls of the unit cell.
- **An aproximation of PolyJet DM 8430**: used in the internal resonator.
- **AS4/3501-6 composite**: a carbon fiber/epoxy material with anisotropic properties, used in the top and bottom constraining layers.
- An **empty (void) region** that represents zones with no material contribution (modeled with zero stiffness).

This anisotropic material is defined by explicitly specifying its element-wise stiffness matrix, which can be adapted by the user for other composite laminates or fiber orientations.

---

##  Repository Structure

```bash
FEM-3D-Heterogeneous
├── main.m                     # Main script to run the simulation
├── DM_Mid.mat                 # DM 8430 properties in the frequency domain
├── Res_top.mat                # Different resonators topologies
├── ElementStiffness.m         # Defines the elementary stiffness matrix for the isotropic materials
├── ElementMass.m              # Defines the elementary mass matrix for the isotropic materials
├── orthocl.m                  # Defines the constitutive matrix, D, for the anisotropic materials (by the Professor Marcelo Trindade)
├── ElementStiffness_comp.m    # Defines the elementary stiffness matrix for the anisotropic materials
