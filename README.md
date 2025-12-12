# HEMEX-HP

This repository contains a MATLAB implementation of the **HEmodynamic and MEtabolic eXchange (HEMEX) Model** for fitting of dynamic Hyperpolarized (HP) <sup>13</sup>C MR images of metabolites. Currently only two simultaneous metabolites are supported (e.g. Pyruvate and Lactate). A 4-metabolite (Pyruvate, Lactate, Bicarbonate, Alanine) implementation is currently being worked on (see code/HEMEX_main_4met.m).

A biological illustration of the model parameters can be seen here:

<p align="center">
  <img src="docs/Model_illustration.png" width="80%">
</p>

An overview of the implementation:

<p align="center">
  <img src="docs/Implementation_overview.png" width="60%">
</p>

A dataset example of Pig Kidneys can be found under *data/*. The supplied MATLAB GUI under *code/* gives a quick and easy way to test the model fitting on <sup>13</sup>C data, as exemplified here:

<img src="docs/GUI_example.png" width=100% height=100%>

## Usage

A userguide for the GUI is available under *docs/*.

A simulation example of the model can be run in the **example_usage_Simulation.m** script.

## Requirements

```
- Optimization Toolbox  
- Image Processing Toolbox  

Code tested on MATLAB version 24.2 (R2024b)
```

## Citing this Work

If you use this repository in your research or publication, please cite our paper:

**N. V. Christensen, M. Redda, N. Bøgh, E. S. S. Hansen, S. Jespersen, and C. Laustsen, "A Combined Hemodynamic and Metabolic Exchange (HEMEX) Model for In Vivo Hyperpolarized 13C MRI," Magnetic Resonance in Medicine (2025): 1–12, https://doi.org/10.1002/mrm.70182.**