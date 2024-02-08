"""  
Read in results from step 1 & 2 and Metaxa2, and combine results into one table.                                                  
"""   


#   Info  #
__author__ = 'Dandan Zhao'
__email__ = 'd.zhao@dal.ca'
__version__ = '1.0.0'
#   End Info   #

import os, sys, re              
import argparse
from Bio import SeqIO
import subprocess
import pandas as pd
import numpy as np
from pandas import DataFrame, ExcelWriter


def parse_fasta(fasta_file):
    print("\tNow begin to parse fasta file: %s" % fasta_file)
    fasta_length_dict = {}
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

    for record in fasta_sequences:
        fasta_length_dict[record.id] =len(record.seq)
    df_contig_length = DataFrame.from_dict(fasta_length_dict, orient='index')
    df_contig_length.rename(columns= {0:'length'}, inplace=True)


    return fasta_sequences, df_contig_length


def parse_Metaxa2_taxonomy(infile):
    print("\tNow begin to parse Metaxa2_taxonomy file: %s" % infile)
    df_Metaxa2 = pd.read_csv(infile, sep='\t', header=None)
    df_Metaxa2_cols = 'qseqid classification pident length score'.strip().split(' ')
    df_Metaxa2.columns = df_Metaxa2_cols
    df_Metaxa2[['pident','length','score' ]]=df_Metaxa2[['pident','length','score' ]].apply(pd.to_numeric)
    df_Metaxa2.index = df_Metaxa2['qseqid']
    df_Metaxa2 = df_Metaxa2.drop('qseqid', axis=1)
    # filter for results with pident >= 80.0 and length >= 300
    df_filtered = df_Metaxa2[(df_Metaxa2['pident'] >= 80.0) & (df_Metaxa2['length'] >= 300)]
    return df_filtered



def combine_Metaxa2_result(metaxa2_dir):
    metaxa2_dict = {}
    for file in os.listdir(metaxa2_dir):
        if file.endswith("_LSU.taxonomy.txt"):
            df_LSU_filtered = parse_Metaxa2_taxonomy(os.path.join(metaxa2_dir,file))
            df_LSU_filtered.columns=['tax_LSU','ident_LSU','length_LSU','score_LSU']

        if file.endswith("_SSU.taxonomy.txt"):
            df_SSU_filtered = parse_Metaxa2_taxonomy(os.path.join(metaxa2_dir,file))
            df_SSU_filtered.columns=['tax_SSU','ident_SSU','length_SSU','score_SSU']

    df_combined_metaxa2 = pd.concat([df_LSU_filtered,df_SSU_filtered], axis =1, sort=False)
    print("\tFinish combining metaxa2 results.")    
    return df_combined_metaxa2


def parse_plast_output(single_plout):
    colnames = ['query ID', 'sseqID','ident','align_length', 'e-value', 'bit score']   
    print("\tNow begin to parse plast_output results: %s" % single_plout)
    plouts_filtered = pd.read_csv(single_plout, sep='\t', header=None, names=colnames)
    plouts_filtered.index = plouts_filtered['query ID']
    plouts_filtered = plouts_filtered.drop('query ID', axis=1)
    print("plouts_filtered\n:",plouts_filtered[0:3])
    return plouts_filtered

def parse_taxonomy_results(infile):
    colnames = ['query ID', 'Domain','species']   
    print("\tNow begin to parse plast_taxonomy results: %s" % infile)
    taxonomy_results = pd.read_csv(infile, sep='\t', header=None, names=colnames)
    taxonomy_results.index = taxonomy_results['query ID']
    taxonomy_results = taxonomy_results.drop('query ID', axis=1)
    return taxonomy_results

def parse_binning_result(infile):
    print("\tNow begin to parse binning results: %s" % infile)
    binning_results = pd.read_csv(infile, sep='\t', header=0, index_col="seqID")
    return binning_results

def parse_Eukfinder_output(infile):
    print("\tNow begin to parse Eukfinder results: %s" % infile)
    df_Eukfinder = pd.read_csv(infile, sep='\t', header=0)
    df_Eukfinder = df_Eukfinder.loc[:,['readID', 'Group']]
    df_Eukfinder.index = df_Eukfinder['readID']
    df_Eukfinder = df_Eukfinder.drop('readID', axis=1)
    df_Eukfinder = df_Eukfinder[df_Eukfinder['Group']=='Eukaryota']
    return df_Eukfinder

def parse_coverage(infile):
    print("\tNow begin to parse depth file: %s" % infile)
    df_coverage = pd.read_csv(infile, sep='\t', header=0)
    df_coverage = df_coverage.loc[:,[ 'contigName','contigLen', 'totalAvgDepth']]
    df_coverage.index = df_coverage['contigName']
    df_coverage = df_coverage.drop('contigName', axis=1)
    df_coverage = df_coverage.astype({"contigLen": int})
    df_coverage = df_coverage[df_coverage['contigLen'] >= 1000]
    return df_coverage	

def parse_Centrifuge_output(infile):
    print("\tNow begin to parse Centrifuge results: %s" % infile)
    df = pd.read_csv(infile, sep='\t', header=0)
    df = df.loc[:,['readID','seqID', 'taxID', 'score', 'hitLength']]
    df = df.groupby('readID', as_index=False, sort=False).first()
    df.index = df['readID']
    df = df.drop('readID', axis=1)
    df = df[(df['seqID'] != 'unclassified')]
    return df
	
def Main(infile, metaxa2_dir, plast_infile, binning_infile, eukfinder_results, depth, centrifuge_result, taxonomy_file):
    print("\nBegin to parse results for %s:\n" % infile)
    suffix = infile.split(".")[-1]
    basename=infile[:infile.index(suffix)]
    fasta_out = parse_fasta(infile)
    coverage_out = parse_coverage(depth)
    binning_out = parse_binning_result(binning_infile)
    ctfg_out = parse_Centrifuge_output(centrifuge_result)
    Metaxa2_out = combine_Metaxa2_result(metaxa2_dir)
    print('Metaxa2_out',Metaxa2_out[0:2])


    if eukfinder_results != '0':
        Eukfinder_out = parse_Eukfinder_output(eukfinder_results)
    else:
        Eukfinder_out = pd.DataFrame({'A' : [np.nan]})
    print("\nEukfinder_out", Eukfinder_out[0:2])


    if plast_infile != '0':
        plast_out = parse_plast_output(plast_infile)
    else:
        plast_out = pd.DataFrame({'B' : [np.nan]})
    print('\nplast_out',plast_out[0:2])


    if taxonomy_file != '0':
        taxonomy_out = parse_taxonomy_results(taxonomy_file)
    else:
        taxonomy_out = pd.DataFrame({'C' : [np.nan]})
    print('\ntaxonomy_out',taxonomy_out[0:2])

    print("\t Begin panda concat.\n")
    df_combined = pd.concat([fasta_out[1], coverage_out, Eukfinder_out, ctfg_out, plast_out, taxonomy_out, binning_out,  Metaxa2_out ], axis =1, sort=False) 
    df_combined = df_combined[df_combined['length'] > 1000 ]
    #df_combined = df_combined[df_combined['MyCC_4mer'] != '']

    #Write results to output file.
    output = '%s_combined_results.xlsx' % basename
    print("\tWrite results to %s." % output)
    writer = ExcelWriter(output)
    df_combined.to_excel(writer, 'results')
    writer.save()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Combine all results into one excel file.", 
                                     epilog='''Usage example:''',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i","--infile",help="Input fasta file (required)", required=True)
    parser.add_argument("-m","--metaxa2",help="Input directory containing all metaxa2 results ('.' or 'Metaxa2_results', or input folder name)", required=False)
    parser.add_argument("-p","--plast",help="Input plast_parsed_results file from 01.py (=0, if no file)", required=True)
    parser.add_argument("-b","--binning",help="Input binning result file (default: binning_results.tsv)", required=False)
    parser.add_argument("-d","--depth",help="Input depth file (default: depth.txt)", required=False)
    parser.add_argument("-e","--eukfinder",help="Input eukfinder result file (=0, if no file)", required=True)
    parser.add_argument("-c","--centrifuge",help="Input centrifuge result file (required)", required=True)
    parser.add_argument("-t","--taxonomy",help="Input taxonomy result file (=0, if no file)", required=True)
    args = parser.parse_args()
    if args.metaxa2:
        metaxa2_dir = args.metaxa2
    elif os.path.isdir("Metaxa2_results"):
        metaxa2_dir = "Metaxa2_results"
    else:
        cwd = os.getcwd()
        metaxa2_dir = cwd 

    if args.binning:
        binning_results = args.binning
    elif os.path.isfile("binning_results.tsv"):
        binning_results = "binning_results.tsv"
    else:
        sys.exit("No binning results")


    if args.depth:
        depth_results = args.depth
    elif os.path.isfile("depth.txt"):
        depth_results = "depth.txt"
    else:
        sys.exit("No depth.txt")

    Main(args.infile, metaxa2_dir, args.plast, binning_results, args.eukfinder, depth_results, args.centrifuge, args.taxonomy )
    #Main(infile, metaxa2_dir, plast_infile, binning_infile, eukfinder_results, depth, centrifuge_result, taxonomy_file)

    print("Finished.")
