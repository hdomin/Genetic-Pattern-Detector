import sys
from lib import search_pattern



# parameters:
#   pattern -> file with pattern to search
#   input_path -> path to content files directory
#   output_file -> file to harvest sequences
if __name__ == "__main__":
    if len(sys.argv) != 3 and False:
        print( "Usage: python harvestSequence.py <pattern_file> <input_path_files> <output-file>")
        sys.exit(1)

    pattern_file = "pattern.txt"
    input_path = "/home/appuser/aln"
    input_file = "1 BS002929.1 and 499 other sequences.aln"
    output_file = "out.txt"

    search_pattern( pattern_file,  input_path, input_file, output_file )


