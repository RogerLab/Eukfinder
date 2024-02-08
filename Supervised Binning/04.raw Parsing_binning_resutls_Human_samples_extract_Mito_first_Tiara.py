'''
Parsing Human gut metagenomic samples for Blastocystis genomes

'''

#   Info  #
__author__ = 'Dandan Zhao'
__email__ = 'd.zhao@dal.ca'
__version__ = '1.0.0'
#   End Info   #

import os, sys, re               
import argparse
import glob
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from pandas import DataFrame

# Read in Subtypes for each human sample
def readingfile(infile):
    with open(infile) as In:
        dict = {}
        for line in In:
            line = line.split('\t')
            dict[line[0]] = int(line[-1])
    #print(dict)
    return dict

# Read in fasta with Contig name and length
def parse_fasta(fasta_file):
    #print("Parse fasta file: %s" % fasta_file)
    fasta_length_dict = {}
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

    for record in fasta_sequences:
        fasta_length_dict[record.id] =len(record.seq)
    df_contig_length = DataFrame.from_dict(fasta_length_dict, orient='index')
    df_contig_length.rename(columns= {0:'length'}, inplace=True)
    df_contig_length['Contig'] = df_contig_length.index
    df_contig_length = df_contig_length[['Contig', 'length']]
    df_contig_length['length'] = df_contig_length['length'].astype(int)
    #print("df_contig_length", df_contig_length.shape)
    return fasta_sequences, df_contig_length


# Parse MyCC results    
def parse_binning(dirs):
    #print("Parse MyCC results:")
    df = pd.DataFrame({'A' : [np.nan]})
    kmer_dict = {}
    #print("MyCC_dirs:\n\t%s " %dirs)
    for dir in dirs:
        cluster_list = []
        parameters = dir.split("_")
        if "mer" in parameters[-3]:
            kmer = parameters[-3] 
        else:
            print("'mer' not in parameters[-3], change to [-3] if having cov in dir_name")    
            sys.exit()            
        kmer_full_name = "_".join(parameters[-4:-1])
        #print("\tNow check dir: ", dir)
        binning_fastas = glob.glob("%s/*.fasta"%(dir))
        fasta_Cluster_dict = {} # A dict containing IDs and Cluster number for each contig        
        for file in binning_fastas:
            file_short = file.split("/")[-1]
            #print(file, file_short)
            cluster = str(file_short.split('.')[1])
            ker_cluster = "%s_bin%s" % (kmer, cluster)
            cluster_list.append(ker_cluster)
            fasta_sequences = SeqIO.parse(open(file),'fasta')
            for record in fasta_sequences:
                
                fasta_Cluster_dict[record.id] = ker_cluster
        kmer_dict[kmer_full_name] = cluster_list
        df_Cluster_mer = DataFrame.from_dict(fasta_Cluster_dict, orient='index')
        df_Cluster_mer.rename(columns= {0:'MyCC_%s'%kmer_full_name}, inplace=True)
        df = pd.concat([df, df_Cluster_mer], axis =1, sort=False)     

    df = df.drop('A', axis=1)
    #print(kmer_dict)
    #print("df_binning", df.shape)
    return df, kmer_dict

 

def parse_Centrifuge_output(infile):
    #print("Parse Centrifuge results: %s" % infile)
    df = pd.read_csv(infile, sep='\t', header=0)
    df = df.loc[:,['readID','seqID', 'taxID', 'score', 'hitLength']]
    df = df.astype({'taxID': int, 'hitLength': int})
    df = df.groupby('readID', as_index=False, sort=False).first()
    df.index = df['readID']
    df = df.drop('readID', axis=1)
    df = df[(df['seqID'] != 'unclassified')]
    #df = df[(df['hitLength'] >= 40)]
    #print("df_Centrifuge", df.shape)
    return df

def parse_plast_output(single_plout):
    colnames = ['query ID', 'sseqID','ident','align_length', 'e-value', 'bit score']   
    #print("\tNow begin to parse plast_output results: %s" % single_plout)
    plouts_filtered = pd.read_csv(single_plout, sep='\t', header=None, names=colnames)
    plouts_filtered.index = plouts_filtered['query ID']
    plouts_filtered = plouts_filtered.drop('query ID', axis=1)
    #print("plouts_filtered:",plouts_filtered.shape)
    return plouts_filtered

def parse_taxonomy_results(infile):
    colnames = ['query ID', 'Domain','species']   
    #print("\tNow begin to parse plast_taxonomy results: %s" % infile)
    taxonomy_results = pd.read_csv(infile, sep='\t', header=None, names=colnames)
    taxonomy_results.index = taxonomy_results['query ID']
    taxonomy_results = taxonomy_results.drop('query ID', axis=1)
    return taxonomy_results

def parse_coverage(infile):
    #print("\tNow begin to parse depth file: %s" % infile)
    df_coverage = pd.read_csv(infile, sep='\t', header=0)
    df_coverage = df_coverage.loc[:,[ 'contigName', 'totalAvgDepth']]
    df_coverage.index = df_coverage['contigName']
    df_coverage = df_coverage.drop('contigName', axis=1)
    #print("df_coverage:",df_coverage.shape)
    return df_coverage	

def parse_Metaxa2_taxonomy(infile):
    #print("\tNow begin to parse Metaxa2_taxonomy file: %s" % infile)
    df_Metaxa2 = pd.read_csv(infile, sep='\t', header=None)
    df_Metaxa2_cols = 'qseqid classification pident length score'.strip().split(' ')
    df_Metaxa2.columns = df_Metaxa2_cols
    df_Metaxa2[['pident','length','score' ]]=df_Metaxa2[['pident','length','score' ]].apply(pd.to_numeric)
    df_Metaxa2.index = df_Metaxa2['qseqid']
    df_Metaxa2 = df_Metaxa2.drop('qseqid', axis=1)
    # filter for results with pident >= 75.0 and length >= 300
    df_filtered = df_Metaxa2[(df_Metaxa2['pident'] >= 75.0) & (df_Metaxa2['length'] >= 300)]
    #print(df_filtered)
    return df_filtered
   
#qseqid sseqid pident length mismatch evalue bitscore stitle
def parse_BLAST_Mit(infile):
    #print("\tNow begin to parse BLAST_Mit file: %s" % infile)
    df = pd.read_csv(infile, sep='\t', header=None)
    df_cols = 'qseqid sseqid_Mit pident_Mit length_Mit mismatch_Mit evalue_Mit bitscore_Mit stitle_Mit'.strip().split(' ')
    df.columns = df_cols
    df = df.loc[:,['qseqid','sseqid_Mit', 'pident_Mit', 'length_Mit', 'evalue_Mit']]

    #Group results by adding hit length
    aggregation_functions = {'sseqid_Mit':'first','pident_Mit':'mean','length_Mit':'sum', 'evalue_Mit':'first'}
    df = df.groupby('qseqid', as_index=False).aggregate(aggregation_functions).reindex(columns=df.columns)

    df[['pident_Mit','length_Mit' ]]=df[['pident_Mit','length_Mit' ]].apply(pd.to_numeric)
    df['sseqid_Mit'] = df['sseqid_Mit'].astype(str)
    df_filtered = df[(df['pident_Mit'] >= 90.0) & (df['length_Mit'] >= 500)]
    df_filtered['organelle'] = 'Mitochondrion'
    list_Mito_seqid = df_filtered['qseqid'].tolist()
    df_filtered.index = df_filtered['qseqid']
    df_filtered =  df_filtered.drop('qseqid', axis=1)

    #print("type for df_filtered",df_filtered.dtypes)
    #print("df_BLAST_Mit:\n", df_filtered)
    return list_Mito_seqid, df_filtered

    
def Main(infile, file_suffix, desired_ID, min_percentage, Domain, Species, List_SeqIDs):
    file_basename = infile[:infile.index(".%s" %file_suffix)]  
    print("\nBegin to parse results for %s:\n" % infile)

    # Parse fasta sequence file, results: fasta_sequences, df_contig_length
    fasta_out = parse_fasta(infile)

    # Parse MyCC binning files, results: df, kmer_dict

    binning_dir = glob.glob("%s_*mer*"%file_basename)
    binning_out = parse_binning(binning_dir)
    # Parse Centrifuge output, results:
    ctfg_file = "%s_centrifuge_UP" %file_basename
    ctfg_out = parse_Centrifuge_output(ctfg_file)

    # Parse PLAST output (after running 01-4.Read_Plast_results_read_map_acc2tax_nt.py)
    # Results: plouts_filtered; taxonomy_results
    plast_file1 = "%s_plastout_parsed_results.txt" %file_basename
    plast_file2 = "%s_plastout_taxonomy_results.txt" %file_basename
    Plast_out = parse_plast_output(plast_file1)    
    Plast_taxonomy  = parse_taxonomy_results(plast_file2)

    # Parse depth file, result: df_coverage
    depth_file = "%s_depth.txt" %file_basename
    depth_out = parse_coverage(depth_file)

    Mito_contigs = []
    # Parse BLAST_Mitochondrial results:
    BLAST_Mit_file = "%s_BLAST4Mit.out" %file_basename
    Mit_size = os.path.getsize(BLAST_Mit_file)
    if Mit_size > 0:
        BLAST_Mit = parse_BLAST_Mit(BLAST_Mit_file)
        #print("data type for BLAST_Mit: ", BLAST_Mit.dtypes)
        #print(BLAST_Mit[0])
        Mito_contigs = BLAST_Mit[0]
    else:
        Mito_contigs = []


    # Parse Metaxa2 results:
    LSU_file = "Metaxa2_results/%s_metaxa2_LSU.taxonomy.txt" %file_basename
    SSU_file = "Metaxa2_results/%s_metaxa2_SSU.taxonomy.txt" %file_basename
    df_LSU_filtered = parse_Metaxa2_taxonomy(LSU_file)
    df_LSU_filtered.columns=['tax_LSU','ident_LSU','length_LSU','score_LSU']
    #print(df_LSU_filtered)

    df_SSU_filtered = parse_Metaxa2_taxonomy(SSU_file)
    df_SSU_filtered.columns=['tax_SSU','ident_SSU','length_SSU','score_SSU']
    #print(df_SSU_filtered)

    metaxa2_out = pd.concat([df_LSU_filtered,df_SSU_filtered], axis =1, sort=False)
    #print(metaxa2_out)

    #List of contigs for Mitochondria: contig_Mito
    df_Mito = metaxa2_out[metaxa2_out['tax_LSU'].astype(str).str.contains('Mitochondria')] 
    #contig_Mito = list(df_Mito.head().index)[0]
    #print("Mitochondrial contigs: ", contig_Mito)

    SSU_contigs = []
    df_Blastocystis_LSU = metaxa2_out[metaxa2_out['tax_LSU'].astype(str).str.contains('Blastocystis')] 
    contig_Blastocystis_LSU = list(df_Blastocystis_LSU.head().index)[0]

    df_Blastocystis_SSU = metaxa2_out[metaxa2_out['tax_SSU'].astype(str).str.contains('Blastocystis')] 
    df_Eukaryota_SSU = metaxa2_out[metaxa2_out['tax_SSU'].astype(str).str.contains('Eukaryota')] 
    contig_Blastocystis_SSU = list(df_Blastocystis_SSU.head().index)[0]
    contig_Eukaryota_SSU = list(df_Eukaryota_SSU.head().index)[0:]
    #print(df_Eukaryota_SSU)
    #print("contig_Eukaryota_SSU: ", contig_Eukaryota_SSU)
    SSU_contigs.append(contig_Blastocystis_SSU)
    #print("SSU_contigs", SSU_contigs)

    SSU_fasta = '%s_SSU.fas' % file_basename
    fasta_sequences = SeqIO.parse(open(infile),'fasta')

    Eukaryota_SSU_dict = df_Eukaryota_SSU.to_dict('index')
    #print(Eukaryota_SSU_dict)
    #print(Eukaryota_SSU_dict.keys())

    result_Eukaryota_SSU_dict = {}
    for item in Eukaryota_SSU_dict.keys():
        #print(item, type(item))
        #print(Eukaryota_SSU_dict[item], type(Eukaryota_SSU_dict[item]))
        #print(Eukaryota_SSU_dict[item]['ident_SSU'], Eukaryota_SSU_dict[item]['length_SSU'])

        species = Eukaryota_SSU_dict[item]['tax_SSU'].split(";")[-2]
        ident = Eukaryota_SSU_dict[item]['ident_SSU']
        length = Eukaryota_SSU_dict[item]['length_SSU']
        line = "%s\tident=%s\tlength=%s" %(species, ident, length)
        result_Eukaryota_SSU_dict[item] =line
        
    f0 = open(SSU_fasta, "w") 
    for seq in fasta_sequences:
        if seq.id in contig_Eukaryota_SSU:
            seq.id = "%s\t%s" %(seq.id, result_Eukaryota_SSU_dict[seq.id])
            SeqIO.write([seq], f0, "fasta")
    f0.close()

    df_combined = pd.concat([fasta_out[1], depth_out, binning_out[0], Plast_out, Plast_taxonomy, ctfg_out, BLAST_Mit[1], metaxa2_out], axis =1, sort=False) #BLAST_Mit, 
    df_combined = df_combined[df_combined['organelle'] != 'Mitochondrion']

    #df_combined = df_combined.dropna()
    #print("df_combined", df_combined.shape)
    output0 = '%s_raw_results.tsv' % file_basename
    #print("\tWrite parsed raw results to file:\t%s\n" % output0)
    df_combined.to_csv(output0, sep='\t', index=False)



    #Select contigs based on Plast_NewDB results only

    df_Plast_selected = df_combined[ df_combined['species'].astype(str).str.contains('Blastocystis')] 
    print("df_Plast_selected: ",df_Plast_selected.shape)

    #Select contigs based on centrifuge results only
    df_centrifuge_selected = df_combined[df_combined['taxID'].isin(List_SeqIDs) ] 
    print("df_centrifuge_selected: ", df_centrifuge_selected.shape)

    df_selected = pd.concat([df_Plast_selected, df_centrifuge_selected]).drop_duplicates()
    print("df_selected: ", df_selected.shape)

    column_names = ["bin_names", "# contigs", "# selected"]
    df_bins = pd.DataFrame(columns = column_names)

    # Bins that have >50% of contigs are selected based on Plast_NewDB and centrifugeresults
    selected_bins = []

    #print("binning_out[1]", binning_out[1])
    for kmer in binning_out[1].keys():
        #print(kmer, binning_out[1][kmer])
        for bin_name in binning_out[1][kmer]:
            value1 = len(df_combined[df_combined["MyCC_%s"%kmer] == bin_name])     
            value2 = len(df_selected[df_selected["MyCC_%s"%kmer] == bin_name])   
            p = round(value2/value1, 2)
            if p > min_percentage:
                selected_bins.append(bin_name)
            #print(bin_name,value1,value2, p )      
    print("selected_bins: ", selected_bins)


    output_contigs_CP = df_selected['Contig'].tolist() 
    #print("\t# in output_contigs_CP", len(output_contigs_CP) )
    output_CP = '%s_CP_extracted_results.tsv' % file_basename
    #print("\nWrite parsing CP results to file:\t%s\n" % output_CP)

    df_extracted_CP = df_selected.loc[output_contigs_CP]

    #Select contigs based on all assessments 
    min_depth = df_selected["totalAvgDepth"].min()
    max_depth = df_selected["totalAvgDepth"].max()
    df_combined_depth_selected = df_combined[df_combined["totalAvgDepth"] <= max_depth]
    #print("Depth range: (%s, %s)" %( min_depth , max_depth))

    df_combined_depth_selected = df_combined_depth_selected.fillna(0)
    output_contigs = []
    manual_check_contigs = []

    result_dict = df_combined_depth_selected.to_dict('records')

    for item in result_dict:
        n = 0
        for kmer in binning_out[1].keys():
            if item["MyCC_%s"%kmer] in selected_bins:
                n = n + 1
        item["# present"] = n
        #print(type(item["sseqid_Mit"]),type(item["species"]))
        if n >= 2:
            output_contigs.append(item["Contig"])
        elif n == 1:
            #print("n=1: ",item["Contig"])
            manual_check_contigs.append(item["Contig"])
        elif item["Contig"] in output_contigs_CP:

            if item["hitLength"] > 300 and int(item["taxID"]) == desired_ID:
                output_contigs.append(item["Contig"])
            elif item["sseqID"]:
                if 'Blastocystis' in item["species"] and item["align_length"] > 300.0 and item["ident"] > 80.0:
                    output_contigs.append(item["Contig"])
                else:
                    #print(item["Contig"], n)
                    manual_check_contigs.append(item["Contig"])


    #Write results to output file.
    print("\t# in output_contigs", len(output_contigs) )

    df_new = DataFrame.from_dict(result_dict, orient='columns')

    output = '%s_extracted_results.tsv' % file_basename
    df_extracted = df_combined.loc[output_contigs]
    df_extracted.to_csv(output, sep='\t', index=False)
    result_fasta = '%s_extracted.fas' % file_basename
    Mito_fasta = '%s_Mito.fas' % file_basename


    df_manual_check = df_combined.loc[manual_check_contigs]
    output_manual_check = '%s_manual_check.tsv' % file_basename
    df_manual_check.to_csv(output_manual_check, sep='\t', index=False)
    fasta_manual_check = '%s_manual_check.fas' % file_basename

    print("\tWrite extracted fasta to file:\t%s\n" % result_fasta)
    fasta_sequences = SeqIO.parse(open(infile),'fasta')
    #print("seqid in Mito_contigs:\n%s\n" % Mito_contigs)
    print("\tWrite Mito_contigs to file:\t%s\n" % Mito_fasta)
    with open(result_fasta, "w") as f1,  open(fasta_manual_check, "w") as f2, \
         open(Mito_fasta, "w") as f3:
        for seq in fasta_sequences:
            #print(seq.id)
            if seq.id in output_contigs:
                SeqIO.write([seq], f1, "fasta")
            elif seq.id in manual_check_contigs:
                SeqIO.write([seq], f2, "fasta")
            elif seq.id in Mito_contigs:
                SeqIO.write([seq], f3, "fasta")

'''
    df_extracted_CP.to_csv(output_CP, sep='\t', index=False)
    result_fasta_CP = '%s_CP_extracted.fa' % file_basename
    #print("\tWrite extracted fasta to file:\t%s\n" % result_fasta_CP)
    fasta_sequences = SeqIO.parse(open(infile),'fasta')
    with open(result_fasta_CP, "w") as f :
        for seq in fasta_sequences:
            #print(seq.id)
            if seq.id in output_contigs_CP:
                SeqIO.write([seq], f, "fasta")


        if 'Blastocystis' in str(item["sseqid_Mit"]) and 'Blastocystis' in item["species"] and int(item["taxID"]) == desired_ID:
            Mito_contigs.append(item["Contig"])
            continue

'''

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Combine all results into one excel file.", 
                                     epilog='''Usage example:''',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i","--infile",help="Input fasta file", required=False)
    parser.add_argument("-s","--suffix",help="Input file suffix", default='fasta', required=False)
    parser.add_argument("-d","--id",help="Input desired seqID for centrifuge results", required=False)
    parser.add_argument("-m","--min",help="Input minimal percentage for contigs belonging to desired seqID", default=0.5, required=False)

    parser.add_argument("-t","--taxonomy",help="Input desired domain name for taxonomy results", required=False)
    parser.add_argument("-p","--species",help="Input desired species name for taxonomy results", required=False)
    parser.add_argument("-l","--list",help="Input desired list of seqIDs for centrifuge results, separated by ','", required=False)


    args = parser.parse_args()

    if args.taxonomy:   
        domain = args.taxonomy
    else:
        domain = 'Eukaryota'        

    if args.species:   
        species = args.species
    else:
        species = 'Blastocystis'    

    if args.list:   
        list_seqIDs = args.list.split(",")
    else:
        list_seqIDs = ['12967','12968','944036', '944160', '944168', '944170']  

    #print("min_percentage for contigs in a bin belonging to Blastocystis", args.min)

    if args.infile:
        file = args.infile
        suffix = file.split(".")[-1]
        Main(file, suffix, int(args.id), float(args.min), domain, species, list_seqIDs)

    elif args.suffix:
        suffix = args.suffix
        files = glob.glob("*.%s" % suffix)
        dict_ST = readingfile('/home/dzhao/Programming_scripts/Supervised_binning/Human_sample_ST.txt')

        for file in files:
            names =re.split("[_|.]", file)
            basename = names[1]
            print(file, "basename:",basename )
            print( "names:",names )

            ID = dict_ST[basename]
            Main(file, suffix, ID, float(args.min), domain, species, list_seqIDs)

    else:
        print("Please input fasta file for supervised binning.")    
        sys.exit() 


    print("Finished.")

