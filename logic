"""
Demonstration logic for COBRApy: FBA, FVA, knockouts, and visualization.
"""

import cobra
from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion
from cobra.flux_analysis import pfba
from cobra.sampling import sample
import matplotlib.pyplot as plt
import numpy as np

from core import get_demo_model
from utils import print_model_summary, format_fluxes

def demo_fba(model):
    """
    Perform Flux Balance Analysis (FBA) and print the result.
    """
    print("\n--- Flux Balance Analysis (FBA) ---")
    solution = model.optimize()
    if solution.status == "optimal":
        print(f"Optimal objective value: {solution.objective_value:.4f}")
        print("Top fluxes:")
        print(format_fluxes(solution))
    else:
        print("No optimal solution found.")

def demo_fva(model):
    """
    Perform Flux Variability Analysis (FVA) and show ranges for reactions.
    """
    print("\n--- Flux Variability Analysis (FVA) ---")
    fva_result = flux_variability_analysis(model, model.reactions[:10])  # first 10 reactions
    print(fva_result)

def demo_gene_knockout(model):
    """
    Simulate single gene deletions and report growth impact.
    """
    print("\n--- Single Gene Knockout Simulation ---")
    deletions = single_gene_deletion(model, model.genes)
    # Show top 5 most affected
    affected = deletions.sort_values(by="growth", ascending=True).head(5)
    print(affected)

def demo_visualize(model):
    """
    Visualize the flux distribution on a network (simple bar plot of top fluxes).
    """
    print("\n--- Visualization of Fluxes (Top 15) ---")
    solution = model.optimize()
    if solution.status != "optimal":
        print("Cannot visualize: no optimal solution.")
        return
    fluxes = solution.fluxes
    # Take top 15 by absolute value
    top = fluxes.abs().sort_values(ascending=False).head(15)
    plt.figure(figsize=(10, 5))
    plt.barh(top.index, top.values, color='skyblue')
    plt.xlabel('Absolute Flux')
    plt.title('Top 15 Fluxes')
    plt.tight_layout()
    plt.show()

def demo_parsimonious_fba(model):
    """
    Perform pFBA (parsimonious FBA) to minimise total flux while achieving optimum.
    """
    print("\n--- Parsimonious FBA (pFBA) ---")
    solution = pfba(model)
    if solution:
        print(f"Optimal objective value: {solution.objective_value:.4f}")
        print("Total flux sum (pFBA):", sum(abs(solution.fluxes)))
    else:
        print("pFBA failed.")

def demo_sampling(model):
    """
    Sample the solution space (if model is feasible).
    """
    print("\n--- Sampling the Solution Space ---")
    try:
        samples = sample(model, 100)
        print(f"Sampled {len(samples)} points.")
        print("First few rows:")
        print(samples.head())
    except Exception as e:
        print(f"Sampling failed: {e}")
