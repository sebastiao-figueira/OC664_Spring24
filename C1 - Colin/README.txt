The working_model_maggie_emily_colin.py combines the ripple, current related friction factor, and wave related friction factor codes, as they all involve an iterative calculation.
The combined wave friction and ripples part requires:
  T_cu: partial crest period
  T_c: crest half cycle period
  T_tu: partial trough period
  T_t: trough half cycle period
  a_hat: 
  u_hat: a representative orbital velocity
  u_delta: current velocity at the boundary
  delta: boundary layer thickness
  d50: median grain size
  d90: 90th percentile grain size
  s: specific gravity [VDA13 defines this differently than some authors]
  g: gravitational constant

The function outputs:
  shields_aa: the time averaged absolute value of the Shields parameter
  f_wc: the wave related friction factor for the crest half cycle of a wave period
  f_wt: the wave related friction factor for the trough half cycle of a wave period
  ksdelta: current related bed roughness
  ksw: wave related bed roughness
  fdelta: also a product of Maggie's code, which is README'd separately in her folder
  f_w: general wave friction factor from the well-known Swart equation
