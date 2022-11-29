## TODO: define evolutionary models (JC, K2P, HKY85, GTR)
import random

"""
User defined constants
"""
## Number of consecutive generations for which the average distance must be above or equal to 0.75
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

"""
Simple function to calculate the average distance given a list of distances by generation
"""
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

# print(simulate_JC(['A','A', 'A', 'A','A', 'A','A', 'A','A', 'A','A', 'A','A', 'A','A', 'A']))


def simulate_K2P(nucleotide_sequence: list) -> list:
     ## User defined constants for the Kimura 2-Parameter Model:
    alpha = 0.07
    beta = 0.05

    """
    This mutation table gives all the probabilities of a nucleotide (first key) mutating to the second nucleotide
    (second key) by K2P rules, but adds upon the previous probability. This makes it easy to progressively check a
    a random number generated for which mutation (or none mutation) should occur
    """

    mutation_table = {
        'A': {
            'A': (1 - alpha - (2 * beta)),
            'C': ((1 - alpha - (2 * beta)) + beta),
            'G': (((1 - alpha - (2 * beta)) + beta) + alpha),
            'T': 1
        },
        'C': {
            'A': beta,
            'C': (beta + ((1 - alpha - (2 * beta)))),
            'G': ((beta + ((1 - alpha - (2 * beta)))) + beta),
            'T': 1
        },
        'G': {
            'A': alpha,
            'C': (alpha + beta),
            'G': ((alpha + beta) + (1 - alpha - (2 * beta))),
            'T': 1
        },
        'T': {
            'A': beta,
            'C': (beta + alpha),
            'G': (beta + alpha + beta),
            'T': 1
        }
    }


    original_sequence = nucleotide_sequence.copy()
    distance_by_generation = [0.0]

    while (calculate_average_distance(distance_by_generation[(-1 * consensus_generations) : ]) < 0.75):

        for i in range (len(nucleotide_sequence)):

            ## get a random number between 0 and 1 inclusive
            random_float = random.random()
            nucleotide = nucleotide_sequence[i]

            sub_table = mutation_table[nucleotide]

            if (random_float <= sub_table['A']):
                nucleotide_sequence[i] = 'A'
            elif (random_float <= sub_table['C']):
                nucleotide_sequence[i] = 'C'
            elif (random_float <= sub_table['G']):
                nucleotide_sequence[i] = 'G'
            ## if it does not fall in the range of probabilities, then it falls in the probability of mutating to T
            else:
                nucleotide_sequence[i] = 'T'

        ## add the distance to the end of the distance_by_generation list
        distance_by_generation.append(calculate_distance(original_sequence, nucleotide_sequence))

    return distance_by_generation


print(simulate_K2P(['A','A', 'A', 'A','A', 'A','A', 'A','A', 'A','A', 'A','A', 'A','A', 'A']))

# def simulate_HKY85(nucleotide_sequence: list) -> list:
#     return

# def simulate_GTR(nucleotide_sequence: list) -> list:
#     return
