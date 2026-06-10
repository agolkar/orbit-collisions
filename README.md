# OrbitCollisions: LEO mega-constellation collision risk

Live demo: https://orbit-collisions.pages.dev (interactive model, plain-language guide,
technical report). Model and report prepared by Claude Fable 5 (Anthropic) under the
supervision of Alessandro Golkar (TUM); AI-generated and currently undergoing
verification. License: MIT. Feedback: golkar@tum.de.

Stage 1 deliverables (analytic models, convergence study, validation). Stage 2
(interactive 3D web app) follows on explicit request and will reuse the engine as-is.

## Contents

- `report.pdf` / `report.tex`: full analytic treatment in English. Both models,
  convergence proof in refined form, correction factors, Monte Carlo validation,
  sensitivity and cascade analysis.
- `collision_models.py`: the computation engine. `Inputs` dataclass (all free
  parameters) strictly separated from derived quantities. Implements the kinetic
  baseline (Kessler particle-in-a-box), the Keplerian kinetic model (latitude
  density field, heading-based velocity structure, rate integral), distributions
  for plotting, shell-thickness inversion, and the cascade branching criterion.
  `python3 collision_models.py` prints the anchor checks (2608 collisions/yr baseline,
  1941/yr Keplerian for the reference case).
- `montecarlo_validation.py`: Cube estimator (Liou et al.) and direct conjunction
  counting with exact two-body propagation. Direct counting validated the analytic
  rate at 1.003 +- 0.032 (954 events over 7 simulated days, N = 1000, R_c = 5 km).
- `analysis_figures.py`: regenerates `figures/*.png` and `results_summary.json`;
  every number in the report comes from these outputs.
- `mc_cube_summary.json`, `mc_direct_summary.json`, `mc_results.json`: raw Monte
  Carlo results.

## Headline numbers (reference case: N = 80,000, A = 120 m^2, sigma = 4A, 500-800 km)

| Quantity | Kinetic baseline | Keplerian |
|---|---|---|
| Natural collisions/yr | 2608 | 1941 |
| Per-satellite rate /yr | 0.065 | 0.049 (0.042-0.062 by inclination) |
| Rate-effective v_rel | 10 km/s (assumed) | 6.1 km/s (derived) |
| Mean impact speed | 10 km/s | 10.2 km/s |

Ratio 0.744 = F_spatial (1.22) x F_velocity (0.61) at 0.5 deg inclination dispersion;
the product is robust, the split depends logarithmically on the dispersion.
Cascade branching number kappa ~ 400 (supercritical); kappa = 1 near N ~ 200 at this
cross-section. Geometry cannot fix the problem (E_acc = 1/yr needs a shell beyond GEO).

## Web app (orbit_collisions_app.html)

Single-file interactive explorer (Phase 3): live parameters, both models, 3D
constellation with risk-gradient coloring (light blue nominal, flashing red at
high local collision-rate density, collision markers at the natural rate in
simulated time), parameter-sweep and distribution charts, and a "Satellite
scale" view comparing the AI1 radiator footprint to a city bus, a tennis
court, an ISS solar wing, and Hubble. Orbit-regime switch: free inclination
mix, or all-SSO (inclination locked to altitude by the J2 sun-synchronous
condition, satellites spread evenly across the band). SSO-only raises the
natural rate ~39% over the Starlink-like mix (2689 vs 1941 per yr): near-polar
planes cross at large angles. JS model layer is a verified 1:1 port of
collision_models.py.

## Requirements

Python 3 with numpy, scipy, matplotlib. Report builds with `pdflatex report.tex` (run twice).
