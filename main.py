from validate_input import *
from evolutionary_models import *
from plot import *

def main():
    """
    This is the main function of the program.
    """

    ## Get the nucleotide sequence
    nucleotide_sequence = validate_input()

    length = len(nucleotide_sequence)

    pi_A = nucleotide_sequence.count('A') / length
    pi_C = nucleotide_sequence.count('C') / length
    pi_G = nucleotide_sequence.count('G') / length
    pi_T = nucleotide_sequence.count('T') / length

    print("Nucleotide distribution: ", {
          "A": pi_A, "C": pi_C, "G": pi_G, "T": pi_T})

    JC_distances =[]
    K2P_distances = []
    HKY85_distances = []
    GTR_distances = []


    for i in range(20):
        ## Get the list of distances of each generation using the Jukes-Cantor Evolutionary model
        JC_distance = simulate_JC(nucleotide_sequence)

        plot_data(JC_distance, "Jukes-Cantor simulation for HIV Gag (" +
                str(len(JC_distance)) + " generations)", "blue", 1)

        ## Get the list of distances of each generation using the Kimura 2 Parameter Evolutionary model
        K2P_distance = simulate_K2P(nucleotide_sequence)
        plot_data(K2P_distance, "K2P simulation for HIV Gag (" +
                str(len(K2P_distance)) + " generations)", "green", 1)

        ## Get the list of distances of each generation using the HKY85 Evolutionary model
        HKY85_distance = simulate_HKY85(nucleotide_sequence)
        plot_data(HKY85_distance, "HKY85 simulation for HIV Gag(" +
                str(len(HKY85_distance)) + " generations)", "red", 1)

        ## Get the list of distances of each generation using the General Time Reversible Evolutionary model
        GTR_distance = simulate_GTR(nucleotide_sequence)
        plot_data(GTR_distance, "GTR simulation for HIV Gag(" +
                str(len(GTR_distance)) + " generations)", "orange", 1)

        plot_all_data(JC_distance, K2P_distance, HKY85_distance, GTR_distance, "Overlay of evolutionary model simulations", 1)

        JC_distances.append(len(JC_distance))
        K2P_distances.append(len(K2P_distance))
        HKY85_distances.append(len(HKY85_distance))
        GTR_distances.append(len(GTR_distance))

    print(JC_distances)
    print(K2P_distances)
    print(HKY85_distances)
    print(GTR_distances)



main()