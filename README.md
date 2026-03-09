# HEMEX-HP

This repository contains a MATLAB implementation of the **HEmodynamic and MEtabolic eXchange (HEMEX) Model** for fitting of dynamic Hyperpolarized (HP) <sup>13</sup>C MR images of metabolites. Currently only two simultaneous metabolites are supported (e.g. Pyruvate and Lactate). A 4-metabolite (Pyruvate, Lactate, Bicarbonate, Alanine) implementation is currently being worked on (see *code/HEMEX_main_4met.m*).

A biological illustration of the model parameters can be seen here:

<p align="center">
  <img src="docs/Model_illustration.png" width="80%">
</p>

An overview of the implementation:

<p align="center">
  <img src="docs/Implementation_overview.png" width="60%">
</p>

## Usage

A dataset example of Pig Kidneys can be found: *data/Pig_kidney_demo.mat* 

If you have your data in your MATLAB environment, simply call the HEMEX_main_2met function (*code/HEMEX_main_2met.m*) to fit the HEMEX model:
```
model_fit_results = HEMEX_main_2met(AIF_raw, pyr_norm, lac_norm, TR, flip_P, flip_L)
```
**AIF_raw**: mean pyruvate signal of AIF over time flip-angle corrected (Mz)

**pyr_norm**: mean pyruvate signal in ROI or voxel over time flip-angle corrected (Mz)

**lac_norm**: mean lactate signal in ROI or voxel over time flip-angle corrected (Mz)

**TR**: repetition time (s)

**flip_P**: flip-angle on pyruvate (°)

**flip_L**: flip-angle on lactate (°)

Visualize results:
```
subplot(1,2,1)
plot(model_fit_results.t,model_fit_results.AIF_raw,'.k','markersize',15)
hold on
plot(model_fit_results.t,model_fit_results.AIF_fit,'-k','linewidth',2)
xlabel('Time (s)')
ylabel('Signal (A.U.)')
legend('AIF','AIF fit')

subplot(1,2,2)
plot(model_fit_results.t,model_fit_results.pyr_norm,'.r','markersize',15);
hold on
plot(model_fit_results.t,model_fit_results.lac_norm,'.b','markersize',15);
plot(model_fit_results.t,model_fit_results.pyr_fits,'-r','linewidth',2);
plot(model_fit_results.t,model_fit_results.lac_fits,'-b','linewidth',2);
xlabel('Time (s)')
ylabel('Signal (A.U.)')
legend('Pyr','Lac','Pyr fit','Lac fit')
s = model_fit_results;
title_str = sprintf(['$k_{PL}=%.3f$ $s^{-1}$, $k_{LP}=%.3f$ $s^{-1}$, $R_{P}=%.3f$ s, $R_{L}=%.3f$ s,\n' ...
                     '$k=%.3f$ $s^{-1}$, $t_0=%.3f$ s, $\\mu=%.3f$ s, $\\sigma=%.3f$ s, $R^2=%.3f$\n'], ...
                     s.kpl_fits, s.klp_fits, s.r1p_fits, s.r1l_fits, ...
                     s.k_fits, s.t0_delay_fits, s.mu_fits, s.sigma_fits, s.R2);
title(title_str,'Interpreter','latex')
```

A simulation example of the model can be found in the *code/example_usage_Simulation.m* script.

The supplied MATLAB GUI *code/HEMEX_GUI_main.m* gives a quick and easy way to test the model fitting on <sup>13</sup>C data. A userguide for the GUI is available: *docs/GUI_usage_guide.pdf*. Example usage (*code/example_usage_GUI.m*):

<img src="docs/GUI_example.png" width=100% height=100%>


## Requirements

```
- Optimization Toolbox  
- Image Processing Toolbox  

Code tested on MATLAB version 24.2 (R2024b)
```

## Citing this Work

If you use this repository in your research or publication, please cite our paper:

**N. V. Christensen, M. Redda, N. Bøgh, E. S. S. Hansen, S. Jespersen, and C. Laustsen, "A Combined Hemodynamic and Metabolic Exchange (HEMEX) Model for In Vivo Hyperpolarized 13C MRI," Magnetic Resonance in Medicine (2025): 1–12, https://doi.org/10.1002/mrm.70182.**