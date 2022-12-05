"""
Conrad, Justin, and Khoa

This file holds all the funstions that simulate the different evolutionary models: Jukes-Cantor, Kimura 2-Parameter,
HKY85, and GTR. Each simulation accepts the initial nucleotide sequence and will run until the average distance of
a user difined number of consensus genrations is greater than or equal to 0.75 (purely random sequence).
"""

import random

"""
User defined constants
"""


## Number of consecutive generations for which the average distance must be above or equal to user-defined purely random
## genetic distance
consensus_generations = 5
## User defined stopping point for purely random genetic distance
pr_genetic_distance = 0.72
## User defined constants for the Jukes-Cantor Model for HIV1:
JC_alpha = 7 * .0001
## User defined constants for the Kimura 2-Parameter Model for HIV1:
K2P_alpha = 0.0006
K2P_beta = 0.0001
# #User defined constants for the HKY85 Model for HIV1:
HKY85_alpha = 0.0006
HKY85_beta = 0.0001
## User defined constants for the General Time Reversible Model for HIV1:
# alpha_AG
# alpha_CT
# beta_AC
# beta_AT
# beta_CG
# beta_GT


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
    original_sequence = nucleotide_sequence.copy()
    distance_by_generation = [0.0]

    while (calculate_average_distance(distance_by_generation[(-1 * consensus_generations):]) < pr_genetic_distance):

        for i in range(len(nucleotide_sequence)):

            # get a random number between 0 and 1 inclusive
            random_float = random.random()
            nucleotide = nucleotide_sequence[i]

            sub_table = mutation_table[nucleotide]

            if (random_float <= sub_table['A']):
                nucleotide_sequence[i] = 'A'
            elif (random_float <= sub_table['C']):
                nucleotide_sequence[i] = 'C'
            elif (random_float <= sub_table['G']):
                nucleotide_sequence[i] = 'G'
            # if it does not fall in the range of probabilities, then it falls in the probability of mutating to T
            else:
                nucleotide_sequence[i] = 'T'

        # add the distance to the end of the distance_by_generation list
        distance_by_generation.append(calculate_distance(
            original_sequence, nucleotide_sequence))

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
            'G': 1 - JC_alpha,
            'T': 1
        }
    }

    return simulatate_genetic_evolution(mutation_table, nucleotide_sequence)


"""
Function to simulate the Kimura 2 Parameter Evolutionary model
Takes an initial nucleotide sequence and a number of generations to simulate
"""
def simulate_K2P(nucleotide_sequence: list) -> list:
    ## User defined constants for the Kimura 2-Parameter Model:
    K2P_alpha = 0.0006
    K2P_beta = 0.0001

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

    return simulatate_genetic_evolution(mutation_table, nucleotide_sequence)


"""
Function to simulate the HKY85 Evolutionary model
Takes an initial nucleotide sequence and a number of generations to simulate
"""

def simulate_HKY85(nucleotide_sequence: list) -> list:

    # calculate the frequency ratios based on the frequency of each nucleotide in the sequence
    length = len(nucleotide_sequence)
    pi_A = nucleotide_sequence.count('A') / length
    pi_C = nucleotide_sequence.count('C') / length
    pi_G = nucleotide_sequence.count('G') / length
    pi_T = nucleotide_sequence.count('T') / length

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

    return simulatate_genetic_evolution(mutation_table, nucleotide_sequence)


"""
Function to simulate the General Time Reversible Evolutionary model
Takes an initial nucleotide sequence and a number of generations to simulate
"""
def simulate_GTR(nucleotide_sequence: list) -> list:

    length = len(nucleotide_sequence)

    # calculate the frequency ratios based on the frequency of each nucleotide in the sequence
    pi_A = nucleotide_sequence.count('A') / length
    pi_C = nucleotide_sequence.count('C') / length
    pi_G = nucleotide_sequence.count('G') / length
    pi_T = nucleotide_sequence.count('T') / length

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

    return simulatate_genetic_evolution(mutation_table, nucleotide_sequence)
