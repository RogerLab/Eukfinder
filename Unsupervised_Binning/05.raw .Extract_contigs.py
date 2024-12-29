"""  
05.raw .Extract_contigs.py
Extract configs from fasta file.  
"""                                       
#   Info  #
__author__ = 'Dandan Zhao'
__email__ = 'd.zhao@dal.ca'
__version__ = '1.0.0'
#   End Info # 
                                         
import sys,re                                                          
import argparse
from Bio import SeqIO

def Main(fasta_file, wanted,suffix, basename):
    #print("wanted:\n",wanted)
    result_file1 = fasta_file.replace('.%s'%suffix,'_%s_extracted.%s'%(basename, suffix))
    print("Extract %s contigs to %s" %(len(wanted), result_file1))

    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    with open(result_file1, "w") as f:
        for seq in fasta_sequences:
            #print(seq.id)
            if seq.id in wanted:
                #print(seq.id)
                SeqIO.write([seq], f, "fasta")
				
    print("Finished")

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Extract configs from fasta file. \nInput: txt file with list of sequence IDs , or sequence IDs separated by , or /.", 
                                     epilog='''Usage example:
                                           python3 Extract_contigs.py -i <fasta_file> [-l <list_Ids> | -s <list_of_seq>]''',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i","--infile",help="Input fasta file", required=True)
    parser.add_argument("-l","--list",help="Input interesting sequence IDs, one per line, ending with '.txt'", required=False)
    parser.add_argument("-s","--seq",help="Input the headers of the sequence ID wanted to retract, separated by ', 'or '/' ", required=False)


    args = parser.parse_args()
    fasta_file = args.infile
    suffix = fasta_file.split(".")[-1]
    #print(suffix)

    if args.seq:
        wanted = set()	
        wanted = re.split('[ |,|/]', args.seq)
        #print("seq list:", repr(wanted))
        print("Extract %s sequence(s) from %s." %(len(wanted),fasta_file))
        basename= 'seq'
        Main(fasta_file, wanted, suffix,basename)
    elif args.list:
        list_file = args.list
        basename = list_file.split(".")[0]
        wanted = set()
        wanted = [re.split('[ |\t|\n]', line)[0] for line in open(args.list) if line != ""]
        print("Extract %s sequences from %s." %(len(wanted), fasta_file))
        Main(fasta_file, wanted, suffix, basename)
    elif args.list == '' and args.seq == '':
        sys.exit("Input a list (-l) or a hearder (-s) for the sequence to be extracted\n")


