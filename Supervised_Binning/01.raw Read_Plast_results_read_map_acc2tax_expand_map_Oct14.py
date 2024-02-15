"""  
Read in a BLAST or PLAST result and write the wanted columns to into a new file.                                                  
"""                                       

#   Info  #
__author__ = 'Dandan Zhao'
__email__ = 'd.zhao@dal.ca'
__version__ = '1.0.0'
#   End Info   #
                                          
import sys,re, time
import argparse
import pandas as pd
from Bio import Entrez
Entrez.email = "d.zhao@dal.ca"
import subprocess  

def Read_map(file):
    print("Read in map file.")
    with open(file) as infile:
        dict = {}
        lines = [line.strip() for line in infile if not line == '']
        for line in lines:
            if line != '':
                #print(repr(line))
                line_split = line.split('\t')
                seq_accession = line_split[0]
                dict[seq_accession] = line_split[1], line_split[2]
    #print("Map dict:", len(dict))

    return dict

def writing_list(list, list_file):
    print("\nNow writing unknown sseqID to file: %s.\n" %list_file)
    with open(list_file,'w') as out:
        for ID in list:
            out.write("%s\n" %ID )

def writing_dict(file,  dic):
    print("\nNow writing results file %s.\n" %file)
    with open(file, 'w') as out:
        for key,value in dic.items():
            #print("NewLine for map_file: ", repr(value))
            out.write("%s\t%s\t%s\n" %(key, value[0],value[1]) )


def Main(result_file):
    map_file = '/scratch3/rogerlab_databases/ds_dz_shared/PlastDB_full_Jun2020/PlastDB_Jun2020_map.txt'
    map_dict = Read_map(map_file)
    suffix = result_file.split(".")[-1]
    out_file = result_file.replace(".%s"%suffix, "_parsed_results.txt")
    out_file2 = result_file.replace(".%s"%suffix, "_taxonomy_results.txt")
    df = pd.read_csv(result_file, sep='\t', header=None)

    df.columns = ['query ID', 'subject ID', 'ident', 'align_length', 'nb. misses', 'nb. gaps',
                    'query begin', 'query end', 'subject begin', 'subject end', 'e-value', 'bit score' ,'exact_align','query length', 
                    'query frame', 'query translated', 'query coverage', 'nb. gaps in query ',
                    'subject length', 'subject frame', 'subject translated', 'subject coverage','nb. gaps in subject']

    plouts = df.groupby('query ID', as_index=False, sort=False).first()
    subset = plouts.loc[:,['query ID', 'subject ID', 'ident', 'align_length', 'e-value', 'bit score' ,'exact_align','query length', 'query coverage']]
    df_filtered = subset[(subset['query coverage'] >= 50.0) & (subset['align_length'] >= 200) & (subset['ident'] >= 90.0)]
    print("Finish subset plast results from %s.\n" %result_file)
    df_filtered.to_csv(out_file, header=0, sep='\t', mode='w')
    print("Length of PLAST results:" , len(df_filtered))
    result_dict = df_filtered.to_dict('records')
    dict_taxo_info = {}      # dict  for sseqID in the map
    dict_subject_query = {}  # dict for sseqID not in the map, key = subjectID

    for item in result_dict:
        if item['subject ID'] in map_dict.keys():
             dict_taxo_info[item['query ID']] = map_dict[item['subject ID']]
        else:
             print("Not found in map file", item['subject ID'], item['query ID'])
    print("Length of dict_taxo_info results after reading map:" , len(dict_taxo_info))

    print("Length of dict_taxo_info results:" , len(dict_taxo_info))

    writing_dict(out_file2, dict_taxo_info)

    print("Finished")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Read in a PLAST result(using NewDB), and write the results + taxonomy into a new file.", 
                                     epilog='''Usage example:
                                           python3 Extract_contigs.py -i <result_file> -f <1 | 2>''',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i","--input",help="Input BLAST or PLAST result", required=True)
    args = parser.parse_args()

    Main(args.input)

