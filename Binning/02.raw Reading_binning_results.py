"""  
Read in binning results from Maxbin, Metabat2 and MyCC, and combine results into one table.                                                  
"""   

#   Info  #
__author__ = 'Dandan Zhao'
__email__ = 'd.zhao@dal.ca'
__version__ = '1.0.0'
#   End Info   #

# 02.raw Reading_binning_results.py
import os, sys
import argparse
import glob
from Bio import SeqIO
import pandas as pd
import numpy as np
from pandas import DataFrame, ExcelWriter



def parse_Maxbin(Acc_number, cwd):
    print("Begin to parse Maxbin results.")
    dir_path = os.path.abspath(cwd)
    Maxbin_dict = {} # A dict containing IDs and Cluster number for Maxbin binning results 
    for file in os.listdir(cwd):
        if file.startswith(Acc_number) and "MaxBin" in file and file.endswith("fasta"):
            print("\tParse %s." % file)
            cluster = str(file.split('.')[-2])
            file_path = os.path.join(dir_path,file)
            #print("\tOpen %s in %s" % (file,file_path))
            fasta_sequences = SeqIO.parse(open(file_path),'fasta')
            for record in fasta_sequences:
                Maxbin_dict[record.id] = cluster
    print("\t Found %s contigs." %(len(Maxbin_dict)))
    df_Maxbin = DataFrame.from_dict(Maxbin_dict, orient='index')
    df_Maxbin.rename(columns= {0:'MaxBin'}, inplace=True)
    return df_Maxbin


def parse_metabat(Acc_number, cwd):
    print("Begin to parse metabat results.")
    dir_path = os.path.abspath(cwd)
    metabat_dict = {} # A dict containing IDs and Cluster number for Maxbin binning results 
    for file in os.listdir(cwd):
        if file.startswith(Acc_number) and "metabat" in file and file.endswith("fa"):
            print("\tParse %s." % file)
            cluster = str(file.split('.')[-2])
            file_path = os.path.join(dir_path,file)
            #print("\tOpen %s in %s" % (file,file_path))
            fasta_sequences = SeqIO.parse(open(file_path),'fasta')
            for record in fasta_sequences:
                metabat_dict[record.id] = cluster
    print("\t Found %s contigs." %(len(metabat_dict)))
    df_metabat = DataFrame.from_dict(metabat_dict, orient='index')
    df_metabat.rename(columns= {0:'metabat'}, inplace=True)
    return df_metabat 

def parse_dir(dir):
    ker = dir.split('_')[2]
    dir_path = os.path.abspath(dir)
    fasta_Cluster_dict = {} # A dict containing IDs and Cluster number for each contig
    for file in os.listdir(dir):
        if file.startswith("Cluster") and file.endswith("fasta"):
            cluster = str(file.split('.')[1])
            file_path = os.path.join(dir_path,file)
            #print("\tOpen %s in %s" % (file,file_path))
            fasta_sequences = SeqIO.parse(open(file_path),'fasta')
            for record in fasta_sequences:
                ker_cluster = "%s_%s" % (ker,cluster )
                fasta_Cluster_dict[record.id] = ker_cluster
    return fasta_Cluster_dict

def MyCC_results(MyCC_dir):
    print("Begin to parse MyCC results.")
    df = pd.DataFrame({'A' : [np.nan]})
    for root,dirs, files in os.walk(MyCC_dir):
        for dir in dirs:
            if "mer" in dir:
                print("Now check dir: %s" %dir)
                dict_mer = parse_dir(dir)
                df_Cluster_mer = DataFrame.from_dict(dict_mer, orient='index')
                df_Cluster_mer.rename(columns= {0:'MyCC_%s'%dir}, inplace=True)
                df = pd.concat([df, df_Cluster_mer], axis =1, sort=False)	 

    df = df.drop('A', axis=1)
    return df

def Main(fasta_file,MyCC_directory):
    Acc = fasta_file.split("_")[0]
    cwd = os.getcwd()
    df_Maxbin_out = parse_Maxbin(Acc, cwd)
    df_metabat_out = parse_metabat(Acc, cwd)	
    df_MyCC = MyCC_results(MyCC_directory)
    results = pd.concat([df_Maxbin_out,df_metabat_out,df_MyCC], axis =1, sort=False)
    out = "binning_results.tsv"
    print("Write results to %s." % out)
    results.to_csv(out, sep='\t',index_label="seqID")
    print("Finished.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read in a BLAST or PLAST result, translate seqID to Taxonomy, and write the results + taxonomy into a new file.", 
                                     epilog='''Usage example:
                                           python3 Extract_contigs.py -i <result_file> [-m <MyCC directory>]''',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i","--input",help="Input fasta file before binning", required=True)
    parser.add_argument("-m","--mycc",help="Input directory containing all MyCC results", required=False)
    args = parser.parse_args()
    if args.mycc:
        Main(args.input, args.mycc)
    else:
        cwd = os.getcwd()
        Main(args.input, cwd)

