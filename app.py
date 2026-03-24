"""
Main application with a menu to run COBRApy demonstrations.
"""

import cobra
from core import get_demo_model
import logic as demo

def main():
    print("COBRApy Demonstration App")
    print("========================")

    # Load the default model
    model = get_demo_model()
    print(f"Loaded model: {model.id}")

    while True:
        print("\nMenu:")
        print("1. Model Summary")
        print("2. Flux Balance Analysis (FBA)")
        print("3. Flux Variability Analysis (FVA)")
        print("4. Gene Knockout Simulation")
        print("5. Visualize Fluxes")
        print("6. Parsimonious FBA (pFBA)")
        print("7. Sampling")
        print("0. Exit")

        choice = input("Select an option: ").strip()

        if choice == "1":
            demo.print_model_summary(model)
        elif choice == "2":
            demo.demo_fba(model)
        elif choice == "3":
            demo.demo_fva(model)
        elif choice == "4":
            demo.demo_gene_knockout(model)
        elif choice == "5":
            demo.demo_visualize(model)
        elif choice == "6":
            demo.demo_parsimonious_fba(model)
        elif choice == "7":
            demo.demo_sampling(model)
        elif choice == "0":
            print("Exiting.")
            break
        else:
            print("Invalid choice. Please try again.")

        input("\nPress Enter to continue...")

if __name__ == "__main__":
    main()
