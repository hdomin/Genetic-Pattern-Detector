import re
from Bio import AlignIO, SeqIO

def search_pattern( pattern_file,  input_path, input_file, output_file ):    
    alignment = read_aln_file(input_path, input_file)    
    iterate_sequences(input_path, input_file)    

    #cs = consensus_sequence(alignment)
    #print(cs)


def read_aln_file(input_path, input_file):
    fileName = input_path + "/" + input_file
    alignment = AlignIO.read(fileName, 'fasta') 

    return alignment


def iterate_sequencesV2(input_path, input_file):
    regxInicial = r'TTTT' #.*GAACAA'
    regxFinal = r'GAACAA'
    fileName = input_path + "/" + input_file
    i = 1
    with open(fileName, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            # Each 'record' contains information about a single sequence
            sequence_id = record.id
            #sequence_description = record.description            
            matches = list(re.finditer(regxInicial, str(record.seq)))            

            smallest_substring = None
            smallest_id = None   
            sequence = str(record.seq)        

            for match in matches:                                
                temp = re.search(regxFinal, sequence[match.end():])

                if temp:
                    if smallest_substring and temp.end() < smallest_id:
                        smallest_substring = sequence[match.start(): match.start() + temp.end() +4]
                        smallest_id = temp.end()
                    else:
                        smallest_substring = sequence[match.start():match.start()+temp.end() +4]
                        smallest_id = temp.end()

            print(f"ID: {sequence_id}    Sequence: { len(smallest_substring)}    Line: {i}")
            #print( "\n\n" + record.seq )
            i += 1
            

def find_smallest_substring(sequence, regx_inicial, regx_final):
    matches = list(re.finditer(regx_inicial, sequence))
    smallest_substring = None
    smallest_id = None

    for match in matches:
        temp = re.search(regx_final, sequence[match.end():])
        if temp and (smallest_id is None or temp.end() < smallest_id):
            smallest_substring = sequence[match.start(): match.start() + temp.end() + 4]
            smallest_id = temp.end()

    return smallest_substring

def iterate_sequences(input_path, input_file):
    regx_inicial = r'TTTT'
    regx_final = r'GAACAA'
    file_name = input_path + "/" + input_file
    i = 1

    with open(file_name, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            sequence_id = record.id
            sequence = str(record.seq)
            
            smallest_substring = find_smallest_substring(sequence, regx_inicial, regx_final)
            
            print(f"ID: {sequence_id}    Sequence: {(smallest_substring)}    Line: {i}")
            i += 1