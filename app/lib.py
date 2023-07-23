import os
import re
from Bio import AlignIO, SeqIO

def search_pattern( pattern_file,  input_path, output_file ):       
    #alignment = read_aln_file(input_path, input_file)
    regxInicial, leftInicial, regxMiddle, regxFinal, rightFinal = read_pattern_file( pattern_file)
    
    iterate_files(regxInicial, leftInicial, regxMiddle, regxFinal, rightFinal, input_path, output_file)    


def read_aln_file(input_path, input_file):
    fileName = input_path + "/" + input_file
    alignment = AlignIO.read(fileName, 'fasta') 

    return alignment


def read_pattern_file( pattern_file):
    with open(pattern_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line != "":
                ret = line.split("-")
                return re.compile( ret[0]), int(ret[1])+1, re.compile(re.escape( ret[2])), re.compile( ret[3]), int(ret[4]) *-1

def iterate_files(regxInicial, leftInicial, rgxMiddle, rgxFinal, rightFinal, input_path, output_file):

    file_list = os.listdir(input_path)

    for file_name in file_list:
        if os.path.isfile(os.path.join(input_path, file_name)) and file_name.endswith('.aln'):
            file_path = os.path.join(input_path, file_name)
        
            iterate_sequences(regxInicial, leftInicial, rgxMiddle, rgxFinal, rightFinal,file_path, file_name)    

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
    #matches = list(re.finditer(regx_inicial, sequence))
    matches = list(regx_inicial.finditer(sequence))
    smallest_substring = None
    smallest_id = None
    add = len(regx_final.pattern) -1

    for match in matches:
        temp = regx_final.search(sequence[match.end():])
        if temp and (smallest_id is None or temp.end() < smallest_id):
            smallest_substring = sequence[match.start(): match.start() + temp.end() +add]
            smallest_id = temp.end()

    return smallest_substring


def iterate_sequences(regx_inicial, left_inicial, regx_middle, regx_final, right_final,file_path, file_name):    
    i = 1

    with open(file_path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            sequence_id = record.id
            sequence = str(record.seq)
            
            smallest_substring = find_smallest_substring(sequence, regx_inicial, regx_final)
            smallest_substring = smallest_substring[left_inicial: right_final] if smallest_substring else ""
            
            print(f"File:\t{file_name}\tID:\t{sequence_id}\tSequence:\t{(smallest_substring)}\tLong:\t{len(smallest_substring)}\tLine:\t{i}\t{sequence[1:10]}")
            i += 1