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

    ## calculate and print the nucleotide distribution
    pi_A = nucleotide_sequence.count('A') / length
    pi_C = nucleotide_sequence.count('C') / length
    pi_G = nucleotide_sequence.count('G') / length
    pi_T = nucleotide_sequence.count('T') / length

    print("Nucleotide distribution: ", {
          "A": pi_A, "C": pi_C, "G": pi_G, "T": pi_T})

    ## create lists to store the number of generations it took each simulation to reach the genetic distance threshold
    JC_generation_lengths =[]
    K2P_generation_lengths = []
    HKY85_generation_lengths = []
    GTR_generation_lengths = []


    for i in range(20):
        ## Get the list of genetic distances by generation of each generation using the Jukes-Cantor Evolutionary model
        JC_distance = simulate_JC(nucleotide_sequence)
        plot_data(JC_distance, "Jukes-Cantor simulation for YML093W (" +
                str(len(JC_distance)) + " generations)", "blue", 10000000)

        # Get the list of genetic distances by generation of each generation using the Kimura 2 Parameter Evolutionary model
        K2P_distance = simulate_K2P(nucleotide_sequence)
        # Plot the data
        plot_data(K2P_distance, "K2P simulation for YML093W (" +
                str(len(K2P_distance)) + " generations)", "green", 10000000)

        # Get the list of genetic distances by generation of each generation using the HKY85 Evolutionary model
        HKY85_distance = simulate_HKY85(nucleotide_sequence)
        # Plot the data
        plot_data(HKY85_distance, "HKY85 simulation for YML093W (" +
                str(len(HKY85_distance)) + " generations)", "red", 10000000)

        # Get the list of genetic distances by generation of each generation using the General Time Reversible Evolutionary model
        GTR_distance = simulate_GTR(nucleotide_sequence)
        # Plot the data
        plot_data(GTR_distance, "GTR simulation for YML093W (" +
                str(len(GTR_distance)) + " generations)", "orange", 10000000)

        ## Create the overlay plot for each of the models for this simulation
        plot_all_data(JC_distance, K2P_distance, HKY85_distance, GTR_distance, "Overlay of evolutionary model simulations", 1)

        ## append the number of generations to reach the genetic distance threshold for each model for this simulation
        ## to the corresponding list
        JC_generation_lengths.append(len(JC_distance))
        K2P_generation_lengths.append(len(K2P_distance))
        HKY85_generation_lengths.append(len(HKY85_distance))
        GTR_generation_lengths.append(len(GTR_distance))

    ## create a dictionary to store the number of generations it took to reach the genetic distance threshold for each simulation for each model
    generations_by_model = {"JC":JC_generation_lengths, "K2P":K2P_generation_lengths, "HKY85":HKY85_generation_lengths, "GTR":GTR_generation_lengths}
    print("Number of generations for each simulation by model: ", "\n", generations_by_model)

main()