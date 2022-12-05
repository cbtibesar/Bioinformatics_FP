import validate_input
import evolutionary_models

def main():
    """
    This is the main function of the program.
    """
    ## Get the amino acid sequence
    seq = validate_input.validate_input()
    ## Get the list of distances of each generation using the Jukes-Cantor Evolutionary model
    JC_distance = evolutionary_models.simulate_JC(seq)
    ## Get the list of distances of each generation using the Kimura 2 Parameter Evolutionary model
    K2P_distance = evolutionary_models.simulate_K2P(seq)
    ## Get the list of distances of each generation using the HKY85 Evolutionary model
    HKY85_distance = evolutionary_models.simulate_HKY85(seq)
    ## Get the list of distances of each generation using the General Time Reversible Evolutionary model
    GTR_distance = evolutionary_models.simulate_GTR(seq)
    
    print("JC Distance: ", JC_distance)
    print("K2P Distance: ", K2P_distance)
    print("HKY85 Distance: ", HKY85_distance)
    print("GTR Distance: ", GTR_distance)
    
main()