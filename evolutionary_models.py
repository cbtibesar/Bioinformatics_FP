"""
Conrad, Justin, and Khoa

This file holds all the functions that simulate the different evolutionary models: Jukes-Cantor, Kimura 2-Parameter,
HKY85, and GTR. Each simulation accepts the initial nucleotide sequence and will run until the average distance of
a user difined number of consensus genrations is greater than or equal to 0.75 (purely random sequence).
"""

import random

"""
User defined constants
"""
# User defined stopping point for purely random genetic distance
threshold_genetic_distance = 0.74

## The number of consensus generations is the number of consecutive generations which must average at or above the
## genetic distance threshold. This lets the program run a bit past the first instance of exceeding the threshold, which
## is important as the genetic distance should fluctuate around the threshold
consensus_generations = 5


"""
In order to simulate the evolutionary models within reasonable computation time, we must make some assumptions about
the number of mutations over a period of generations. For some of the genes, the mutation rate between individual
generations is too small, so we must speed up the process by adding in some determinism. To do this, we will assume
a mutation rate that is greater over a greater number of generations, so each generation in the produced list will
represent a certain number of generations that have past.
"""
number_of_generations_per_item = 10000000

# User defined constants for the Jukes-Cantor Model for YML093W:
# General mutation rate for the gene
JC_alpha = 6.7 * 0.0001

# User defined constants for the Kimura 2-Parameter and HKY85 models for YML093W:
# Mutation rate caused by transitions
K2P_alpha, HKY85_alpha = (.494 * JC_alpha * 2), (.494 * JC_alpha * 2)
# Mutation rate caused by transversions
K2P_beta, HKY85_beta = (.506 * JC_alpha * 2), (.506 * JC_alpha * 2)

# User defined constants for the General Time Reversible Model for YML093W:
# Mutation rate caused by A <--> G
alpha_AG = (.247 * JC_alpha * 4)
# Mutation rate caused by C <--> T
alpha_CT = (.247 * JC_alpha * 4)
# Mutation rate caused by A <--> C
beta_AC = (.146 * JC_alpha * 8)
# Mutation rate caused by A <--> T
beta_AT = (.063 * JC_alpha * 8)
# Mutation rate caused by C <--> G
beta_CG = (.152 * JC_alpha * 8)
# Mutation rate caused by T <--> G
beta_GT = (.146 * JC_alpha * 8)



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
Function to simulate any of the genetic evolutionary models given the table of mutation probabilities and the
nucleotide sequence
"""

def simulatate_genetic_evolution(mutation_table: dict, nucleotide_sequence: list) -> list:
    ## create a copy of the original sequence to compare the mutated sequences to
    original_sequence = nucleotide_sequence.copy()
    ## create the list to store the genetic distances by generation, starting with the first generation having 0 distance
    distance_by_generation = [0.0]

    ## while the average distance of the last n consensus generations is less than the genetic distance threshold, mutate to the next generation
    while (calculate_average_distance(distance_by_generation[(-1 * consensus_generations):]) < threshold_genetic_distance):

        ## iterate over the length of the nucleotide sequence
        for i in range(len(nucleotide_sequence)):

            ## get a random number between 0 and 1 inclusive
            random_float = random.random()
            ## get the nucleotide at the ith position
            nucleotide = nucleotide_sequence[i]

            ## get the mutation probabilities of that nucleotide for the givben mutation table
            sub_table = mutation_table[nucleotide]

            ## find out what range the random float falls within
            ## if it falls within the range of the probability of the nucleotide mutating (or "mutating" to itself (prob. of not mutating)),
            ## set the nucleotide at that position to the corresponding nucleotide
            if (random_float <= sub_table['A']):
                nucleotide_sequence[i] = 'A'
            elif (random_float <= sub_table['C']):
                nucleotide_sequence[i] = 'C'
            elif (random_float <= sub_table['G']):
                nucleotide_sequence[i] = 'G'
            ## if it does not fall in the range of probabilities, then it falls in the probability of mutating to T
            else:
                nucleotide_sequence[i] = 'T'

        ## calculate the genetic distance of the newly mutated sequence, and append it to the end of the distance_by_generation list
        distance_by_generation.append(calculate_distance(
            original_sequence, nucleotide_sequence))

    ## once the consensus sequences average at or above the threshold, return the list of genetic distances by generation
    return distance_by_generation


"""
Function to simulate the Jukes-Cantor Evolutionary model
Takes an initial nucleotide sequence and a number of generations to simulate
"""
def simulate_JC(nucleotide_sequence: list) -> list:

    """
    This mutation table gives all the probabilities of a nucleotide (first key) mutating to the second nucleotide
    (second key) by JC rules, but adds upon the previous probability. This makes it easy to progressively check a
    a random number generated for which mutation (or none mutation) should occur
    """

    mutation_table = {
        'A': {
            'A': 1 - (3 * JC_alpha),
            'C': 1 - (2 * JC_alpha),
            'G': 1 - JC_alpha,
            'T': 1
        },
        'C': {
            'A': JC_alpha,
            'C': 1 - (2 * JC_alpha),
            'G': 1 - JC_alpha,
            'T': 1
        },
        'G': {
            'A': JC_alpha,
            'C': 2 * JC_alpha,
            'G': 1 - JC_alpha,
            'T': 1
        },
        'T': {
            'A': JC_alpha,
            'C': 2 * JC_alpha,
            'G': 3 * JC_alpha,
            'T': 1
        }
    }


    ## use the generalized simulation function, passing in the JC mutation table and the nucleotide sequence
    return simulatate_genetic_evolution(mutation_table, nucleotide_sequence)


"""
Function to simulate the Kimura 2 Parameter Evolutionary model
Takes an initial nucleotide sequence and a number of generations to simulate
"""
def simulate_K2P(nucleotide_sequence: list) -> list:

    """
    This mutation table gives all the probabilities of a nucleotide (first key) mutating to the second nucleotide
    (second key) by K2P rules, but adds upon the previous probability. This makes it easy to progressively check a
    a random number generated for which mutation (or none mutation) should occur
    """

    mutation_table = {
        'A': {
            'A': 1 - K2P_alpha - (2 * K2P_beta),
            'C': 1 - K2P_alpha - K2P_beta,
            'G': 1 - K2P_beta,
            'T': 1
        },
        'C': {
            'A': K2P_beta,
            'C': 1 - K2P_alpha - K2P_beta,
            'G': 1 - K2P_alpha,
            'T': 1
        },
        'G': {
            'A': K2P_alpha,
            'C': K2P_alpha + K2P_beta,
            'G': 1 - K2P_beta,
            'T': 1
        },
        'T': {
            'A': K2P_beta,
            'C': K2P_beta + K2P_alpha,
            'G': K2P_beta + K2P_alpha + K2P_beta,
            'T': 1
        }
    }

    ## use the generalized simulation function, passing in the K2P mutation table and the nucleotide sequence
    return simulatate_genetic_evolution(mutation_table, nucleotide_sequence)


"""
Function to simulate the HKY85 Evolutionary model
Takes an initial nucleotide sequence and a number of generations to simulate
"""

def simulate_HKY85(nucleotide_sequence: list) -> list:

    # calculate the frequency ratios based on the frequency of each nucleotide in the sequence
    # normalize so that an even distribution of nucleotides correlates to 1
    length = len(nucleotide_sequence)

    pi_A = 4 * (nucleotide_sequence.count('A') / length)
    pi_C = 4 * (nucleotide_sequence.count('C') / length)
    pi_G = 4 * (nucleotide_sequence.count('G') / length)
    pi_T = 4 * (nucleotide_sequence.count('T') / length)

    """
    This mutation table gives all the probabilities of a nucleotide (first key) mutating to the second nucleotide
    (second key) by HKY85 rules, but adds upon the previous probability. This makes it easy to progressively check a
    a random number generated for which mutation (or none mutation) should occur
    """

    mutation_table = {
        'A': {
            'A': 1 - (HKY85_alpha * pi_G) - (HKY85_beta * pi_C) - (HKY85_beta * pi_T),
            'C': 1 - (HKY85_alpha * pi_G) - (HKY85_beta * pi_T),
            'G': 1 - (HKY85_beta * pi_T),
            'T': 1
        },
        'C': {
            'A': HKY85_beta * pi_A,
            'C': 1 - (HKY85_alpha * pi_T) - (HKY85_beta * pi_G),
            'G': 1 - (HKY85_alpha * pi_T),
            'T': 1
        },
        'G': {
            'A': HKY85_alpha * pi_A,
            'C': (HKY85_alpha * pi_A) + (HKY85_beta * pi_C),
            'G': 1 - (HKY85_alpha * pi_T),
            'T': 1
        },
        'T': {
            'A': HKY85_beta * pi_A,
            'C': (HKY85_beta * pi_A) + (HKY85_alpha * pi_C),
            'G': (HKY85_beta * pi_A) + (HKY85_alpha * pi_C) + (HKY85_beta * pi_G),
            'T': 1
        }
    }

    ## use the generalized simulation function, passing in the HKY85 mutation table and the nucleotide sequence
    return simulatate_genetic_evolution(mutation_table, nucleotide_sequence)


"""
Function to simulate the General Time Reversible Evolutionary model
Takes an initial nucleotide sequence and a number of generations to simulate
"""
def simulate_GTR(nucleotide_sequence: list) -> list:

    length = len(nucleotide_sequence)

    ## calculate the frequency ratios based on the frequency of each nucleotide in the sequence
    ## normalize so that an even distribution of nucleotides correlates to 1
    pi_A = 4 * (nucleotide_sequence.count('A') / length)
    pi_C = 4 * (nucleotide_sequence.count('C') / length)
    pi_G = 4 * (nucleotide_sequence.count('G') / length)
    pi_T = 4 * (nucleotide_sequence.count('T') / length)

    """
    This mutation table gives all the probabilities of a nucleotide (first key) mutating to the second nucleotide
    (second key) by GTR rules, but adds upon the previous probability. This makes it easy to progressively check a
    a random number generated for which mutation (or none mutation) should occur
    """

    mutation_table = {
        'A': {
            'A': 1 - (beta_AC * pi_C) - (alpha_AG * pi_G) - (beta_AT * pi_T),
            'C': 1 - (alpha_AG * pi_G) - (beta_AT * pi_T),
            'G': 1 - (beta_AT * pi_T),
            'T': 1
        },
        'C': {
            'A': beta_AC * pi_A,
            'C': 1 - (beta_CG * pi_G) - (alpha_CT * pi_T),
            'G': 1 - (alpha_CT * pi_T),
            'T': 1
        },
        'G': {
            'A': alpha_AG * pi_A,
            'C': (alpha_AG * pi_A) + (beta_CG * pi_C),
            'G': 1 - (beta_GT * pi_T),
            'T': 1
        },
        'T': {
            'A': beta_AT * pi_A,
            'C': (beta_AT * pi_A) + (alpha_CT * pi_C),
            'G': (beta_AT * pi_A) + (alpha_CT * pi_C) + (beta_GT * pi_G),
            'T': 1
        }
    }

    ## use the generalized simulation function, passing in the GTR mutation table and the nucleotide sequence
    return simulatate_genetic_evolution(mutation_table, nucleotide_sequence)
