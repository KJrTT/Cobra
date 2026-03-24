"""
Core module: defines model loading and simple data containers.
"""

import cobra
from cobra.io import load_model

# Attempt to use the built-in E. coli core model from cobra.test
try:
    from cobra.test import create_test_model
    DEMO_MODEL = create_test_model("textbook")
except ImportError:
    # Fallback: load from the local SBML file if available
    import os
    if os.path.exists("e_coli_core.xml"):
        DEMO_MODEL = load_model("e_coli_core.xml")
    else:
        # Minimal synthetic model if nothing else works
        DEMO_MODEL = cobra.Model("minimal")
        # Add a few reactions to make it useful
        atp = cobra.Metabolite("atp_c", name="ATP")
        adp = cobra.Metabolite("adp_c", name="ADP")
        pi = cobra.Metabolite("pi_c", name="Pi")
        h2o = cobra.Metabolite("h2o_c", name="H2O")

        r1 = cobra.Reaction("ATPhydrolysis")
        r1.add_metabolites({atp: -1, h2o: -1, adp: 1, pi: 1})
        r1.bounds = (0, 1000)
        DEMO_MODEL.add_reactions([r1])

        # Biomass reaction placeholder
        biomass = cobra.Reaction("BIOMASS")
        biomass.add_metabolites({atp: -1})
        biomass.bounds = (0, 1000)
        DEMO_MODEL.add_reactions([biomass])
        DEMO_MODEL.objective = "BIOMASS"

def get_demo_model():
    """
    Return a preloaded COBRA model suitable for demonstrations.
    """
    return DEMO_MODEL.copy()
