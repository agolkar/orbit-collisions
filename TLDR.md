# TLDR: Can 80,000 compute satellites coexist without crashing into each other?

**The question.** Orbital data centers need big radiators to dump heat from their chips. Big radiators make big targets. We modeled how often a fleet of 80,000 satellites, each with a 120 m² radiator, flying between 500 and 800 km altitude, would collide if nobody steered.

**The headline numbers.**

- Left unmanaged, the fleet would suffer roughly **2,000 collisions per year**. Each satellite has about a **5% chance per year** of being hit.
- This is a workload figure, a measure of how hard the autopilot system must work. With collision avoidance that succeeds 99.9% of the time, residual losses drop to about 2 per year.
- You cannot fix this with geometry. Spreading the fleet thin enough to get to one collision per year would require filling space out beyond the geostationary belt, 40,000+ km up. The only real levers are active avoidance, smaller cross-sections, and fleet size.
- The stakes compound: at this density, one breakup creates debris likely to cause hundreds of follow-on hits over time (a Kessler cascade). Prevention and reliable disposal of dead satellites are everything; you cannot dodge debris you cannot track.
- Collision burden scales with the **square** of fleet size and linearly with radiator area. Every extra megawatt of compute carries a proportional collision tax.

**How much do we trust this?** We built two independent models: a simple "gas of satellites" estimate and a full orbital-mechanics model, then verified the latter by brute-force simulation of close approaches (within 3% agreement). The refined model says the simple estimate was about 25% pessimistic on the rate, and tells us where collisions happen (clustered at high latitudes where orbital planes cross) and how hard they hit (10 to 15 km/s, almost always shattering).

**Bottom line for an investor.** Mega-scale orbital compute is not blocked by collision physics, but it is gated by operations: the business case must carry an autopilot and traffic-management system of extreme reliability, radiator designs that minimize frontal area, and disposal discipline near 100%. Ask any orbital data center venture three questions: what is your collision-avoidance failure rate including dead satellites, how does your radiator area scale with compute power, and how far is your architecture from the debris-cascade threshold.

*Basis: full analytic treatment, Monte Carlo validation, and parametric engine in report.pdf and collision_models.py (June 2026). Reference scenario predates and roughly matches the scale of constellations now proposed publicly (e.g., an 88,000-satellite orbital data center filing with the FCC in February 2026).*
