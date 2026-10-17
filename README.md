# SMELT / RECLAIM - Microwave Aluminum Melt-and-Pour System

Thermal and drain modeling, test data, CAD, and documentation for **SMELT**, the metals-processing subsystem of **RECLAIM**, a dual metal and plastic recycling system built at Penn State for the NASA LunaRecycle Challenge. RECLAIM won a total of $125,000 in prize awards.

I led a 4-person lab team that designed, built, and tested the physical aluminum recycling and casting system, and I co-wrote the phase-2 and final competition proposals.

## What's here

| Folder | Contents |
|---|---|
| `matlab/` | MATLAB `ode45` models of the melt and pour process |
| `data/` | Melt test logs (raw timestamped, 8-point averaged, full data stream), energy fit, and pump sizing |
| `cad/` | Phase 3 chamber and lid (SolidWorks + Fusion + drawings), assembled chamber STEP/STL, crucible, ingot and mold geometry, bayonet-lock guide |
| `docs/` | SSPL SMELT final report, competition proposal, metals flowdown, P&ID (draw.io source), SOP, pump sizing |
| `figures/` | P&ID render, thermal crucible results, chamber renders |

## Models

**`smelt_batch_melt_and_drain.m`**: a two-node lumped model with a SiC crucible node and an aluminum charge node.
- Microwave power in: `eta*P`, with a provisional coupling efficiency of 0.45.
- Heat moves from the crucible to the Al through a wall conduction resistance `R_int = t/(k*A)`. Heat leaves through insulation conductance `UA_loss`.
- Latent heat of fusion is spread over a ±2.5 °C window as an effective cp.
- Integration stops at the 700 °C pour temperature (event function).
- A Torricelli/orifice drain model gives the drain time through an 8 mm tap hole (Cd = 0.80).

**`smelt_continuous_pour.m`**: a four-state coupled model with solid mass, liquid mass, crucible temperature, and pool temperature.
- Melt-first energy split: heat from the crucible melts the solid charge first, and only what's left over superheats the pool.
- The pool drains while it is molten.

**`smelt_energy_fit.m`**: a least-squares fit of measured run energy against charge mass and 1/power (`E = E0 + e*M + k/P`), with R².

Temperatures are in °C throughout. Specific-heat fits use intercepts shifted from Kelvin fits.

## Tools
MATLAB, SolidWorks, Autodesk Fusion, draw.io, Excel

## Credits
Team project, Penn State SSPL. The final report and proposal were written with my teammates.
