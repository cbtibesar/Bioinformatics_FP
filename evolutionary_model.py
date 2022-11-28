## TODO: define evolutionary models (JC, K2P, HKY85, GTR)
import random

"""
User defined constants
"""
consensus_generations = 5

"""
Function that accepts two DNA sequences and calculates the genetic distance between them,
returning the distance as a decimal of differences over length of the sequences
"""
def calculate_distance(sequence1: list, sequence2: list) -> float:
    differences = 0
    length = len(sequence1)

    for i in range (length):
        if sequence1[i] != sequence2[i]:
            differences += 1

    return differences / length

def calculate_average_distance(partial_distance_by_generation: list) -> float:
    return sum(partial_distance_by_generation) / len(partial_distance_by_generation)


"""
Function to simulate the Jukes-Cantor Evolutionary model
takes an initial nucleotide sequence and a number of generations to simulate
"""
def simulate_JC(nucleotide_sequence: list) -> list:
    ## User defined constants for the Jukes-Cantor Model:
    alpha = .05

    original_sequence = nucleotide_sequence.copy()
    distance_by_generation = [0.0]

    while (calculate_average_distance(distance_by_generation[(-1 * consensus_generations) : ]) < 0.75):

        for i in range (len(nucleotide_sequence)):

            ## get a random number between 0 and 1 inclusive
            random_float = random.random()
            nucleotide = nucleotide_sequence[i]

            ## if the random number is less than or equal the mutation rate, mutate the nucleotide with a random nucleotide other than the original (since each is equally likely in JC Model)
            if random_float <= alpha:
                dna_nucleotides = ['A', 'C', 'G', 'T']
                ## remove the the original nucleotide to so that it does not mutate to itself
                dna_nucleotides.remove(nucleotide)
                ## pick a random nucleotide of the remaining to mutate to
                nucleotide_sequence[i] = dna_nucleotides[random.randint(0, 2)]

        ## add the distance to the end of the distance_by_generation list
        distance_by_generation.append(calculate_distance(original_sequence, nucleotide_sequence))

    return distance_by_generation





print(simulate_JC(['A','A', 'A', 'A','A', 'A','A', 'A','A', 'A','A', 'A','A', 'A','A', 'A']))



# def simulate_K2P(nucleotide_sequence: list) -> list:
#     return

# def simulate_HKY85(nucleotide_sequence: list) -> list:
#     return

# def simulate_GTR(nucleotide_sequence: list) -> list:
#     return
