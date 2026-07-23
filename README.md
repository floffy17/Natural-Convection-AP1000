# Analysis and Design of Passive Decay Heat Removal Systems (DHRS)

![Politecnico di Torino](https://img.shields.io/badge/University-Politecnico_di_Torino-blue.svg)
![Field](https://img.shields.io/badge/Field-Nuclear_Engineering-green.svg)
![Standard](https://img.shields.io/badge/Focus-Passive_Safety_Systems-red.svg)

---

## Overview
This repository contains the computational model, analytical formulations, and calculation scripts developed for **Assignments 1.1 & 1.2** of the **Nuclear Fission Plants** course at **Politecnico di Torino** (A.Y. 2025/2026).

The project focuses on the thermal-hydraulic analysis and design verification of a multi-loop Passive Decay Heat Removal System (DHRS) for nuclear reactors operating under post-shutdown accidental conditions. The study evaluates buoyancy-driven natural circulation mechanism to ensure residual decay power is safely transferred to an ultimate heat sink without requiring active pumps or external electrical power.

## Key Objectives & Analysis Features
- **Natural Circulation Loop Dimensioning:** Sizes a basic single-phase/two-phase rectangular loop by balancing buoyancy driving forces against hydraulic pressure losses.
- **Pipe Roughness & Friction Effects:** Evaluates the influence of standard commercial pipe sizing, wall roughness, and fittings on localized and distributed head losses.
- **Multi-Circuit DHRS Performance:** Verifies the coupled thermal-hydraulic response across the Primary System Circuit, Intermediate System Circuit, and Tertiary System Circuit.
- **Heat Exchanger Design Verification:** Models heat transfer rates and overall thermal resistance across shell-and-tube and submerged in-pool heat exchangers.
- **Thermal Safety Margin Assessment:** Determines steady-state mass flow rates and temperature distributions to ensure coolants remain strictly in the single-phase liquid regime, avoiding boiling[cite: 3].

## Methodology & Standards
1. **Steady-State Thermal-Hydraulic Modeling:** Solves coupled energy and momentum balance equations iteratively using fluid thermophysical properties[cite: 3].
2. **Hydraulic Loss Quantification:** Evaluates friction factors and minor loss coefficients across bends, valves, core geometry, and heat exchanger tube bundles[cite: 3].
3. **Passive System Safety Verification:** Assesses system performance against reactor safety margins to confirm reliable decay heat removal during station blackout or pump failure scenarios[cite: 3].

## Repository Structure
```text
.
├── src/                # Python calculation routines and fluid property wrappers
├── data/               # Input geometry, component loss coefficients, and plant data
├── docs/               # Assignment documentation and technical reference report
├── results/            # Computed flow rates, pressure drop breakdowns, and temperature profiles
└── README.md           # Project summary and documentation
