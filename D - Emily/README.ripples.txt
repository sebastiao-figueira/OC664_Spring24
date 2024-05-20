
# ripple model equations: VDA13 appendix B
# required defined factors: need to be solved for to get desired outputs 
    # ripple length (lambda)
    # ripple height (eta)
    # time-averaged absolute shields stress (named shields_aa here)
        # defined in category C1 (Colin)
    # d50, d90: grain size (grain size, static values assumed)
    # maximum mobility number (psimax)
# INPUTS: 
    # CONSTANTS:
        # d50
        # d90
        # s
        # g
    # Defined elsewhere
        # ahat (category A eqn. 9)
        # shields_aa (category C1 eqn A.3, can rename to fit)
    # Defined here (for outputs)
        # eta (ripple height, eqn. B.1 VDA13)
        # lambda (ripple length, eqn. B.2) 
        # psimax (maximum mobility number for ripple dimension eqns)
            # currently using regular flow- slightly different for irregular
 # OUTPUTS
    # mu (fine factor adjustment for sheet flow)
    # ksdelta (current-related bed roughness)
    # ksw (wave-related bed roughness, iterative with shields_aa from C because
            # mean absolute shields depends on bed roughness)