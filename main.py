"""
Conrad, Khoa, Justin

The purpose of this program is to simulate each of the different evolutionary models in the study (Jukes-Cantor, Kimura 2-Parameter,
HKY85, and General Time-Reversible) until they reach maximum genetic distance (0.75) from the original sequence. Each
simulation will be ran 20 times, with the number of generations it took to reach maximum genetic distance recorded for
each simulation for each model.

IMPORTANT: you must set the model parameters (mutation rates for transitions, transversions, etc.) in the
evolutionary_models.py file! These values are formatted and found in the mutation rates txt file labeled by the name of
the gene.

Once the file has been given by the user, the program will run a simulation for each model 20 times, each time producing
a scatter plot showing the genetic distance as a function of generation, or normalized generation, and the number of
normalized generations required to reach the maximum genetic distance. It will also produce a graph of the scatter plots
of each model overlaid to provide easy comparison for a simulation between models

Finally the program will output the nucleotide frequency distribution for the gene, followed by a dictionary displaying
a list of the number of generations it took each of the 20 simulations to reach the maximum genetic distance by model
type. This information was then used to calculate the mean, standard deviation, and confidence intervals of each of the
models for that particular gene.

Please make sure that you have all the required packages available on you computer when running the program.
I used Spyder3 IDE, as it made it easy to produce the graphs right within the IDE as the program was running. 
"""
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
        JC_distance = simulate_JC(nucleotide_sequence.copy())
        plot_data(JC_distance, "Jukes-Cantor simulation for HIV Gag (" +
                  str(len(JC_distance)) + " generations)", "blue", 1)

        # Get the list of genetic distances by generation of each generation using the Kimura 2 Parameter Evolutionary model
        K2P_distance = simulate_K2P(nucleotide_sequence.copy())
        # Plot the data
        plot_data(K2P_distance, "K2P simulation for HIV Gag (" +
                  str(len(K2P_distance)) + " generations)", "green", 1)

        # Get the list of genetic distances by generation of each generation using the HKY85 Evolutionary model
        HKY85_distance = simulate_HKY85(nucleotide_sequence.copy())
        # Plot the data
        plot_data(HKY85_distance, "HKY85 simulation for HIV Gag (" +
                  str(len(HKY85_distance)) + " generations)", "red", 1)

        # Get the list of genetic distances by generation of each generation using the General Time Reversible Evolutionary model
        GTR_distance = simulate_GTR(nucleotide_sequence.copy())
        # Plot the data
        plot_data(GTR_distance, "GTR simulation for HIV Gag (" +
                str(len(GTR_distance)) + " generations)", "orange", 10000000)

        ## Create the overlay plot for each of the models for this simulation
        plot_all_data(JC_distance, K2P_distance, HKY85_distance, GTR_distance,
                      "Overlay of evolutionary model simulations", 1)

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
