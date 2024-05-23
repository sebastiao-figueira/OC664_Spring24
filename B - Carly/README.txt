entrainedSandLoad.py

Four functions:
1. dimensionless_grainsize: calculates dimensionless grain size D* from Soulsby 1997 given: grain diameter, gravity, sediment density, water density, and water kinematic viscosity. Returns constant D*
2. critical_shields: calculates critical Shields number from Soulsby 1997 given D* (dimensionless grain size). Returns constant theta_cr (critical Shields number).
3. sandload_crest: Calculates sand load for crest half-cycle
4. sandload_trough: Calculates sand load for trough half-cycle