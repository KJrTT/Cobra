"""
Utility functions for COBRA models: formatting, I/O, summaries.
"""

import cobra
from cobra.io import save_json_model, load_json_model

def print_model_summary(model):
    """
    Print basic statistics of a COBRA model.
    """
    print(f"Model: {model.id}")
    print(f"Number of metabolites: {len(model.metabolites)}")
    print(f"Number of reactions: {len(model.reactions)}")
    print(f"Number of genes: {len(model.genes)}")
    print(f"Objective: {model.objective.expression}")
    print("Reaction bounds:")
    for rxn in model.reactions[:5]:  # first five reactions
        print(f"  {rxn.id}: {rxn.lower_bound} - {rxn.upper_bound}")

def format_fluxes(solution, n=10):
    """
    Return a string with the top n fluxes from a solution.
    """
    fluxes = solution.fluxes
    sorted_fluxes = fluxes.abs().sort_values(ascending=False).head(n)
    lines = []
    for rxn_id, flux in sorted_fluxes.items():
        lines.append(f"{rxn_id}: {flux:.4f}")
    return "\n".join(lines)

def save_model(model, filename):
    """
    Save model to JSON file (portable and human-readable).
    """
    save_json_model(model, filename)
    print(f"Model saved to {filename}")

def load_model(filename):
    """
    Load model from JSON file.
    """
    return load_json_model(filename)
