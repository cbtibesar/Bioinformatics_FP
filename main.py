from validate_input import *
from evolutionary_models import *
from plot import *

def main():
    """
    This is the main function of the program.
    """
    ## Get the amino acid sequence
    seq = validate_input()
    ## Get the list of distances of each generation using the Jukes-Cantor Evolutionary model
    JC_distance = simulate_JC(seq)
    ## Get the list of distances of each generation using the Kimura 2 Parameter Evolutionary model
    K2P_distance = simulate_K2P(seq)
    ## Get the list of distances of each generation using the HKY85 Evolutionary model
    HKY85_distance = simulate_HKY85(seq)
    ## Get the list of distances of each generation using the General Time Reversible Evolutionary model
    # GTR_distance = simulate_GTR(seq)
    
    plot_data(JC_distance, "Jukes-Cantor simulation for HIV1", "blue")
    plot_data(K2P_distance, "K2P simulation for HIV1", "green")
    plot_data(JC_distance, "HKY85 simulation for HIV1", "red")
    # print("GTR Distance: ", GTR_distance)

main()