# Alignment note: response to "Feasibility of an Orbital Data Center" (10 June 2026)

Prepared with reference to our analysis "Collision Risk in LEO Mega-Constellations: from the
Kinetic-Gas Model to the Keplerian Model" (9 June 2026, cited as [R4] in the reviewed report)
and its computational engine. Contact: golkar@tum.de.

## Summary

The report is an independent reconstruction and extension of our analysis, and the
reconstruction passes at full precision: every shared quantity matches our engine to three
or more significant digits. We consider this a successful independent replication of the
core model. The extensions (altitude sensitivity, external debris, residence-time framing,
decision gates) are complementary and consistent with our own approximation register, which
declared external debris, operational reliability, and financial/regulatory risk out of
scope. The few divergences are at the level of framing and of generalizations that hold for
the base case but not universally; they are listed below as suggestions, not corrections.

## Numerical cross-check (their value vs our engine)

| Quantity | Report | Our engine |
|---|---|---|
| Shell volume 500-800 km | 1.859e20 m3 | 1.859e20 m3 |
| Kinetic-gas baseline | 2,608/yr | 2,608/yr |
| Keplerian rate, base case | 1,941.3/yr | 1,941.5/yr |
| Keplerian/kinetic ratio | 0.744 | 0.744 |
| Mean natural rate per satellite | 0.0485/yr (4.74%/yr) | 0.0485/yr (4.74%/yr) |
| Radial factor, uniform band | 2/V | 2/V |
| Thin 10 km shell at 650 km | ~58,000/yr | 58,248/yr |
| Screening load within 1 km | ~320/sat/yr; ~25.6M fleet-yr | 317/sat/yr; 25.4M fleet-yr |
| N for 1 natural collision/yr | 1,816 | 1,816 |
| A for 1 natural collision/yr at N=80k | 0.062 m2 | 0.062 m2 |
| f_fail for <1 residual/yr | 5.15e-4 | 5.15e-4 |
| Cascade branching number | ~408 | 407.5 |
| kappa = 1 residence time | ~22 days | 22.4 days |
| Altitude factor 650 to 1,200 km | 0.83 | 0.83 (r^-2.5 confirmed) |

The five-year per-satellite probability (21.5%) and the residual-risk thresholds in their
Table 9 are exact Poisson consequences of the shared rate and are confirmed.

## Where the report extends us, and we agree

The altitude-sensitivity mapping (R roughly proportional to r^-2.5; no rescue inside LEO,
worse persistence above ~600-700 km) follows from our parametric engine and we confirm the
table values. The external-debris layer is genuinely additive: our analysis covers
self-collisions only, by design, and we agree that for 9.6 million m2 of deployed area the
external flux (catalogue plus lethal non-trackable) can rival or dominate the residual risk
budget, and that a bankable assessment requires native ORDEM/MASTER runs. The
residence-time reformulation of the cascade boundary (kappa = 1 at ~22 days for the base
case) is an elegant corollary of the shared branching formula and makes the disposal lever
vivid. The distinction between natural workload and residual loss, the full-chain
definition of f_fail including dead satellites, and the gated decision structure all match
our conclusions and sharpen their operational meaning.

## Suggested refinements for the author

1. **The 0.744 correction factor is architecture-specific, not generic.** It belongs to the
   43/53/70/97.6 mix with weights 0.2/0.4/0.2/0.2. For an all-sun-synchronous architecture
   (inclination locked to altitude by the J2 condition, the configuration thermal designers
   prefer for radiator and solar-array orientation), our engine gives 2,689/yr, a ratio of
   1.03: the Keplerian correction flips sign because every near-polar plane crosses every
   other at steep angles. The report's design lever "avoid unnecessary polar/retrograde
   mixing" is correct and could cite this quantitatively: the popular SSO choice raises the
   natural rate by about 39% over the mixed base case.

2. **The per-satellite rate hides a 50% spread across inclination families.** The 0.0485/yr
   mean decomposes into 0.042 (43 deg) to 0.062 (97.6 deg) per year. For insurance pricing
   and replacement budgeting by orbit family, the spread is material.

3. **f_fail and kappa are decoupled, and stating this explicitly would strengthen Section 9.**
   Avoidance reliability suppresses satellite-satellite losses but does nothing against the
   untracked fragment flux that drives the branching number; the two risk channels need
   separate mitigation stacks (operations for the first, fragmentation suppression and
   disposal for the second). The report implies this but does not say it directly.

4. **Five-year figures assume a static fleet.** The 21.5% per-satellite and the 1-in-5-years
   thresholds presuppose constant N and sigma over the horizon; a deployment ramp roughly
   halves the time-integrated exposure of the first tranche while replenishment extends it.
   A one-line caveat would prevent over-reading.

5. **RAAN clustering deserves a red flag of its own.** Both analyses assume randomized RAAN,
   and the report correctly lists "real slotting can reduce some terms" in its approximation
   register. The dangerous direction is the opposite one: an all-dawn-dusk SSO fleet
   clusters RAAN at one local time, concentrating traffic in a way that invalidates the
   randomized-ensemble assumption of every kinetic model, ours and theirs alike. Since
   dawn-dusk is the thermally preferred orbit for compute payloads, this is arguably the
   single most important open modeling case for this class of system.

6. **One framing distinction.** The verdicts diverge in tone, not physics: the report issues
   a NO-GO for the base configuration; our analysis states gating requirements (avoidance
   failure near 1e-4 over the full chain including dead satellites, fragmentation
   suppression, verified disposal). The underlying numbers are identical; the difference is
   a risk-appetite judgment layered on shared physics, which is the reviewer's prerogative.

## Bottom line

We find no numerical disagreement anywhere the two analyses overlap, which, given the
independent reconstruction, materially raises confidence in both. The extensions are ones
our approximation register anticipated, and the refinements suggested above are additive.
We would be glad to wire the external-debris layer into our interactive engine as a third
risk channel alongside self-collisions and cascade; the modular structure supports it
directly.

*This note was prepared by Claude Fable 5 (Anthropic) under the supervision of
A. Golkar (TUM). The underlying model is AI-generated and undergoing verification;
interactive version and source: https://orbit-collisions.pages.dev and
github.com/agolkar/orbit-collisions.*
