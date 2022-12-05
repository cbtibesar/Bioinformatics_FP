"""
This file gets the user-input file, opens it, and validates that the data is clean. It also capitalizes the nucleotide
sequence for greater consistency.
"""

"""
Function to get the user input from a text file
"""
def read_file() -> str:
   file_name = input("Please enter the name of the file txt file to input: ")
   try:
       file = open(file_name, 'r')

   # check for common exceptions
   except FileNotFoundError:
       print("Error: Could not find file ", file_name)
       return "Error"
   except OSError:
       print("Error: OS error while attempting to open ", file_name)
       return "Error"
   except Exception:
       print("Error: Unexpected error while attempting to open ", file_name)
       return "Error"

   else:
       file_input = file.read()
       if file_input == "":
           print("Error: file is empty")
           return "Error"
       else:
           return file_input

       file.close()


"""
Function to validate the user input DNA sequence to ensure it only
contains legal nucleotides. Returns false if input is not legal, and
returns capitalized sequence if legal input.
"""
def validate_nucleotide_sequence(raw_DNA_input: str) -> str:
   # capitalize the input to allow for capitalization flexibility, and to
   # standardize the input for validation and interpretation

   capitalized_DNA = raw_DNA_input.upper()
   # remove all tabs, spaces, and newlines
   capitalized_DNA = capitalized_DNA.replace("\n", "")
   capitalized_DNA = capitalized_DNA.replace("\t", "")
   capitalized_DNA = capitalized_DNA.replace(" ", "")

   legal_nucleotides = set(['A', 'C', 'G', 'T'])

   for i in range(0, len(capitalized_DNA)):
       # check that each nucleotide is legal
       if capitalized_DNA[i] not in legal_nucleotides:
           print("Invalid DNA sequence: found character '", capitalized_DNA[i],
                 "' at index", i)
           return "Error"

   return capitalized_DNA

def validate_input():
    ## get the file input of initial DNA sequence for each haloid individual
    raw_data = read_file()
    ## validate the input adn return the validated sequence
    return validate_data(raw_data)