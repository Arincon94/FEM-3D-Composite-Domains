# 3D Finite Element Simulation of Heterogeneous Domains in MATLAB

This repository provides a customizable MATLAB implementation for simulating 3D solid structures using the Finite Element Method (FEM). It supports domains composed of elements with different materials and sizes, enabling the modeling of heterogeneous media such as optimized metamaterials, periodic unit cells, and structures with embedded resonating regions.

---

## Motivation

This project was developed to study and prototype advanced structures composed of various materials, using FEM with 8-node isoparametric hexahedral elements. It is particularly suited for modeling structures with spatially varying properties — including materials with frequency-dependent behavior — and offers flexibility to simulate periodic and resonant systems.

While a separate repository contains the theoretical formulation of the finite elements used, this implementation focuses on the **assembly of complex domains**, allowing regions to be defined as different materials, including the possibility of voids (zero-stiffness regions).
