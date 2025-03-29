#!/usr/bin/env python
import json
import os
import re
import sys
import tarfile
from pathlib import Path

import ete3
import glob
import time
import shutil
import platform
import numpy as np
import pandas as pd
from subprocess import PIPE, run
from joblib import Parallel, delayed
import argparse
import urllib.request


#   Info  #
__author__ = 'Dayana E. Salas-Leiva'
__email__ = 'ds2000@cam.ac.uk'
__version__ = '1.2.4'
#   End Info   #

# database info
# NOTE: It contains test datasets to verify the correct installation and functionality of the Eukfinder pipeline
_all_db = ["eukfinder env databases", "72 GB", "https://perun.biochem.dal.ca/Eukfinder/eukfinder_dbs_env_v1.2.4.tar.gz", "eukfinder_dbs_env_v1.2.4.tar.gz"]

_database = {
    "1": ["centrifuge database", "70 GB", "https://perun.biochem.dal.ca/Eukfinder/compressed_db/centrifuge_db.tar.gz", "centrifuge_db.tar.gz"],
    "2": ["PLAST database", "2.1 GB", "https://perun.biochem.dal.ca/Eukfinder/compressed_db/PlastDB_db.tar.gz", "PlastDB_db.tar.gz"],
    "3": ["Human Genome for read decontamination", "0.92 GB", "https://perun.biochem.dal.ca/Eukfinder/compressed_db/GCF_000001405.39_GRCh38.p13_human_genome.fna.tar.gz", "GCF_000001405.39_GRCh38.p13_human_genome.fna.tar.gz"],
    "4": ["Read Adapters for Illumina sequencing", "2.4 KB", "https://perun.biochem.dal.ca/Eukfinder/compressed_db/TrueSeq2_NexteraSE-PE.fa.tar.gz", "TrueSeq2_NexteraSE-PE.fa.tar.gz"],
    "5": ["test set", "51 MB", "https://perun.biochem.dal.ca/Eukfinder/compressed_db/test_files.tar.gz", "test_files.tar.gz"],
    "6": ["test set with tiny dbs", "51 MB", "https://perun.biochem.dal.ca/Eukfinder/compressed_db/test_sample_with_tiny_DB.tar.gz", "test_sample_with_tiny_DB.tar.gz"]
}

_cdb = "Centrifuge_DB/Centrifuge_NewDB_Sept2020"
_pdb = "PLAST_DB/PlastDB.fasta"
_pmap = "PLAST_DB/PlastDB_map.txt"

# JSON
_json_path = f"{os.path.expanduser('~')}/.eukfinder/config.json"

# --- preparation ---
def trimming(bn, reads1, reads2, adapath, wsize, qscore, headcrop,
             mlenght, threads, leading_trim, trail_trim, qenc):

    r1_out = '%sR1PT.fq %sR1unPT.fq ' % (bn, bn)
    r2_out = '%sR2PT.fq %sR2unPT.fq ' % (bn, bn)
    cmd = 'trimmomatic PE -threads %s -trimlog %s.trim.log ' % (threads, bn)
    cmd += '%s %s %s %s ILLUMINACLIP:%s:2:30:10 HEADCROP:%s LEADING:%s ' % (
        reads1, reads2, r1_out, r2_out, adapath, headcrop, leading_trim)
    if qenc == 'auto':
        cmd += 'TRAILING:%s SLIDINGWINDOW:%s:%s MINLEN:%s ' % (trail_trim,
                                                          wsize, qscore,
                                                          mlenght)
    else:
        cmd += 'TRAILING:%s SLIDINGWINDOW:%s:%s MINLEN:%s -%s' % (trail_trim,
                                                          wsize, qscore,
                                                          mlenght, qenc)

    ms = 'trimmomatic cmd_line:\n%s\n' % cmd
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

    _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    lst = r1_out.split() + r2_out.split()
    catfile = '%s.catunPT.fq' % bn
    os.system('cat %s %s > %s' % (lst[1], lst[-1], catfile))
    R1 = os.path.abspath(lst[0])
    R2 = os.path.abspath(lst[2])
    cat = os.path.abspath(catfile)
    return R1, R2, cat


def bowtie2build(hostg):

    bn = os.path.split(hostg)[-1]
    bn = re.split(r'.fasta$|.fa$|.fas$', bn, flags=re.IGNORECASE)[0]
    cmd = 'bowtie2-build -f %s %s' % (hostg, bn)
    ms = 'bowtie2 build cmd_line:\n%s' % cmd
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    if len(glob.glob('*bt2')) != 0:
        return bn
    else:
        ms = 'bowtie indexes could not be built.\n'
        ms += 'Exiting program\n'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)


def bowtie2(bn, threads, r1, r2, catr1r2, hindex, typ):

    if typ == 'fastq-host':
        cmd = 'bowtie2 --local --phred33 -q --threads %s -x %s' % (threads,
                                                                   hindex)
        cmd += ' -1 %s -2 %s -U %s -S %s.sam ' % (r1, r2, catr1r2, bn)
        cmd += '--un-conc %s_p.fastq --un %s_un.fastq --no-unal' % (bn, bn)
        ms = 'bowtie2 fastq-host cmd_line:\n%s' % cmd
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        bt2faR1 = os.path.abspath('%s_p.1.fastq' % bn)
        bt2faR2 = os.path.abspath('%s_p.2.fastq' % bn)
        bt2faun = os.path.abspath('%s_un.fastq' % bn)
        return bt2faR1, bt2faR2, bt2faun, cmd

    elif typ == 'fastq':
        cmd = 'bowtie2 --local --phred33 -q --threads %s -x %s' % (threads,
                                                                   hindex)
        cmd += ' -1 %s -2 %s -U %s -S %s.sam ' % (r1, r2, catr1r2, bn)
        cmd += '--al-conc %s_bowtie2.fastq ' % bn
        cmd += '--al %s_bowtie2.un.fastq --no-unal' % bn
        ms = 'bowtie2 fastq cmd_line:\n%s' % cmd
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        bt2fqR1 = os.path.abspath('%s_bowtie2.1.fastq' % bn)
        bt2fqR2 = os.path.abspath('%s_bowtie2.2.fastq' % bn)
        bt2fqun = os.path.abspath('%s_bowtie2.un.fastq' % bn)
        return bt2fqR1, bt2fqR2, bt2fqun, cmd

    else:
        cmd = 'bowtie2 --local --threads %s -x %s ' % (threads, hindex)
        cmd += '-U %s -S %s.sam --al %s_bowtie2.fasta' % (catr1r2, bn, bn)
        cmd += '--no-unal'
        ms = 'bowtie2 fasta cmd_line:\n%s' % cmd
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        bt2fa = os.path.abspath('%s_bowtie2.fasta' % bn)
        return bt2fa, None, None, cmd


def centrifuge(bn, bn_tuple, threads, mhlen, dbpath, k, pair=True, fastq=True):

    if fastq:
        gline = 'centrifuge -q --phred33 --threads %s -k %s ' % (threads, k)
        gline += '--min-hitlen %s -x %s ' % (mhlen, dbpath)
    else:
        gline = 'centrifuge -f --threads %s -k %s ' % (threads, k)
        gline += '--min-hitlen %s -x %s ' % (mhlen, dbpath)

    if pair:
        bn_r1, bn_r2 = bn_tuple
        report = os.path.join(os.getcwd(), '%s_centrifuge_P' % bn)
        gline += '-1 %s -2 %s -S %s ' % (bn_r1, bn_r2, report)
        gline += '--report-file %s.tsv' % report
        ms = 'centrifuge pair cmd_line:\n%s' % gline
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        return gline, report

    else:
        bn_r1r2 = bn_tuple
        report = os.path.join(os.getcwd(), '%s_centrifuge_UP' % bn)

        gline += '-U %s -S %s --report-file %s.tsv ' % (bn_r1r2,
                                                        report, report)
        ms = 'centrifuge unpair cmd_line:\n%s' % gline
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        return gline, report


def centrifuge_output_checkup(infile, header):

    if os.path.exists(infile) and os.stat(infile).st_size != 0:

         v = head2(infile)
    else:
        m = 'There is something wrong with output file: %s' % infile


def cats(b_outname, dir_path):

    os.chdir('..')
    wd = os.path.join(dir_path, 'Classified_reads')
    # up_dir = cwd.split('/Temp')[0]
    outr1 = os.path.join(wd, '%s.EUnk.R1.fq' % b_outname)
    outr2 = os.path.join(wd, '%s.EUnk.R2.fq' % b_outname)
    outr1r2 = os.path.join(wd, '%s.EUnk.un.fq' % b_outname)
    return outr1, outr2, outr1r2


# ---  classification ---

def reading_reports(infile):
    """

    :param infile:
    :return:
    """

    ms = 'Processing %s ...' % infile
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    if infile is not None:
        sorting_list = ['score']
        chunksize = 1000 * 10
        # header must be guessed to avoid iteration issues
        infox = pd.read_csv(infile, sep='\t',
                            chunksize=chunksize, dtype=object)
        table = Parallel(n_jobs=-2)(delayed(minimal)
                                    ((chunk, sorting_list, 'report', 'None',
                                      'None',
                                      False)) for chunk in infox)
        table = pd.concat(table).sort_values(sorting_list,
                                             ascending=[False])
        table = table.groupby('readID', as_index=False).first()
        table = table.loc[:, ['readID', 'seqID', 'taxID']]
        ms = 'Report has been read'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        return table
    else:
        m = 'report file has not been declared '
        m += 'or the declared file is not in the directory. The file or '
        m += 'a symbolic link must exist in the working directory.'
        print(m, sep=' ', end='\n', file=sys.stdout, flush=True)
        return pd.DataFrame()


def minimal(tupla):
    """

    :param tupla:
    :return:
    """

    # presort... groupby will retain the relative position
    chunk, sorting_list, typ, pid, cov, lr = tupla
    if not lr:
        cov = 0
    if typ == 'plast':
        nchunk = chunk.copy(deep=True)
        small = nchunk.qlen < nchunk.slen  # query length smaller than
        # subject
        nchunk.loc[small, 'sguess_cov'] = nchunk.qcov
        large = nchunk.qlen >= nchunk.slen  # query larger than subject
        nchunk.loc[large, 'sguess_cov'] = nchunk.scov
        # select entries with identities and coverages that meet
        # pre-established criteria
        expr = (nchunk['pidentity'] >= pid) & (nchunk['sguess_cov'] >= cov)
        newchunk = nchunk.loc[expr].copy(deep=True)
        asc = [True, False, False]
        cols = ['query ID', 'subject ID', 'e-value',
                'aln_len', 'pidentity', 'sguess_cov']
        group = 'query ID'

    else:
        group = 'readID'
        asc = False
        cols = ['readID', 'seqID', 'taxID', 'score']
        newchunk = chunk.copy(deep=True)

    if not newchunk.empty:
        dfout = newchunk.sort_values(sorting_list, ascending=asc)
        # do groupby and take the top directly
        dataframes = dfout.groupby(group, as_index=False, sort=False
                                   ).first()
        dataframes = dataframes.reindex(columns=cols)
        return dataframes
    else:
        ms = 'Chunk does not fullfill identity and coverage criteria. '
        ms += 'This produced an empty read chunk'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        print('empty chunk is :\n', newchunk)
        return newchunk


def parseplastoutput(single_plout, ident, cov, lr):
    """

    :param single_plout:
    :return:
    """

    colnames = ['query ID', 'subject ID', 'pidentity', 'aln_length',
                'nb.misses', 'nb.gaps', 'qbegin', 'qend', 'sbegin',
                'send', 'e-value', 'bit score', 'unknown', 'qlen',
                'qframe', 'qtranslated', 'qcov', 'qnb.gaps', 'slen',
                'sframe', 'stranslated', 'scov', 'snb.gaps', 'sguess_cov']

    sorting_list = ['e-value', 'pidentity', 'sguess_cov']

    types = {entry: ('float64' if entry not in ['query ID', 'subject ID']
                     else 'str') for entry in colnames}

    # reading plast outputs by chunks using parallel parsing
    chunksize = 1000 * 10
    df = pd.read_csv(single_plout, sep='\t', header=None,
                     names=colnames, chunksize=chunksize, dtype=types)
    plouts = Parallel(n_jobs=-2)(delayed(minimal)
                                 ((chunk, sorting_list, 'plast', ident,
                                   cov, lr)) for chunk in df)

    # re-sorting as some queries may appear in different chunks
    plouts = pd.concat(plouts).sort_values(
        sorting_list, ascending=[True, False, False])

    # final best single plout
    plouts = plouts.groupby('query ID', as_index=False).first()
    #
    temp = plouts.loc[:, ['subject ID']]
    # acc2tax = acc2tax['subject ID'].str.
    #        split('|', 4, expand=True).loc[:, 3]
    temp = temp.loc[:, 'subject ID'].tolist()

    # database MUST be in ncbi format. an exception should be created here
    # to capture a different ending('|')
    acc2tax = [(entry, entry.rstrip('|').split('|')[-1].split('.')[0])
               if '|' in entry else (entry, entry.split()[0].split(
                '.')[0]) for entry in temp]

    # creating a dataframe with only accessions
    acc2tax = pd.DataFrame(acc2tax,
                           columns=['subject ID', 'Accession'])
    #
    acc2tax = acc2tax.groupby('Accession', as_index=False).first()

    # renaming column to match that of the main dataframe and acc2tax
    plouts = plouts.rename(columns={'query ID': 'readID'})

    # merge main plout output with only the accession numbers
    plouts = plouts.merge(acc2tax, on='subject ID',
                          how='outer')
    # return a dataframe with all columns plus an accession id

    return plouts


def cmd(CPUs, chunk, pDBpath, evalue, long=False):
    """

    :param CPUs:
    :param chunk:
    :param pDBpath:
    :param evalid:
    :return:
    """

    # evalue = float(evalid.split(':')[0])
    line = 'plast -e %s -max-hit-per-query 1 -outfmt 2 '
    line += '-a %s -p plastn -i %s -d %s -o %s.plout '
    line += '-force-query-order 1000 '
    if long:
        line += '-F F\n'
        line %= (evalue, CPUs, chunk, pDBpath, chunk)
    else:
        line += '\n'
        line %= (evalue, CPUs, chunk, pDBpath, chunk)
    ms = 'plast cmd line:\n%s ' % line
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    return line


def my_run(cmdline):

    """

    :param cmdline:
    :return:
    """

    run(cmdline, stderr=PIPE,
        stdout=PIPE, shell=True, bufsize=1,
        close_fds=True, universal_newlines=True)


def plastsearches(CPUs, pDBpath, evalue, pattern):
    """
    :param CPUs:
    :param pDBpath:
    :param evalue:
    :param pattern:
    :return:
    """

    abs_pattern = os.path.join(os.getcwd(), pattern)
    queries = glob.glob(abs_pattern)
    # queries = glob.glob(pattern): i added 2 lines above instead of this one
    label = pattern[0].upper()
    if queries:
        ms = 'plast search in progress ...'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        pSubprocesses = [cmd(CPUs, chunk, pDBpath, evalue) for
                         chunk in queries]
        _ = Parallel(n_jobs=-2)(delayed(my_run)(cmd_line)
                                for cmd_line in pSubprocesses)
        time.sleep(5)
        outpattern = pattern + '.plout'
        res = glob.glob(outpattern)
        res = [(e, label) for e in res]
        return res


def parsingplastsearches(list_of_results, pid, cov, lr):
    """

    :param list_of_results:
    :return:
    """

    time.sleep(5)
    print('list_of_results', list_of_results)
    if len(list_of_results[0]) != 0:
        label = list_of_results[0][1]
        ms = 'parsing plast outputs for %s ....' % label
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        allplast_outs = []
        pouts = [parseplastoutput(out[0], pid, cov, lr) for out in
                 list_of_results if (os.path.exists(out[0]) and
                                     os.stat(out[0]).st_size != 0)]
        if pouts:
            plast_outputs = pd.concat(pouts)
            allplast_outs.append((label, plast_outputs))
        else:
            allplast_outs.append((label, pd.DataFrame()))
        return allplast_outs
    else:
        print('Plast results are empty. exiting program until bug is fixed')
        sys.exit(-1)


def matchmaker(ncbi, regex, readid, taxid):
    """

    :param ncbi:
    :param regex:
    :param readid:
    :param taxid:
    :return:
    """


    lineage = ncbi.get_lineage(taxid)
    Se_lineage = pd.Series(ncbi.get_taxid_translator
                           (lineage), name=taxid)
    nl = Se_lineage.iloc[1]
    match = regex.search(nl)
    exc = ['archaea', 'bacteria', 'eukaryota']
    if match:
        if not match.group().lower() in exc:
            match_word = 'Vir'
        else:
            match_word = match.group()
        return readid, match_word
    else:
        return readid, np.nan


def binningbytaxonomy(report_df):
    """
    slice reports by taxonomic domain
    """

    ncbi = ete3.NCBITaxa()
    regex = re.compile(r'^\bBacteria\b$|^\bArchaea\b$|'
                       r'^\bEukaryota\b$|vir',
                       flags=re.IGNORECASE)
    maintaxidlist = []
    for entry in report_df.itertuples():
        index, readid, seqid, taxid = entry
        taxid = str(taxid)
        if taxid != '0':
            try:
                m = matchmaker(ncbi, regex, readid, taxid)
                if type(m) is tuple:
                    maintaxidlist.append(m)
            except:
                pass
    return maintaxidlist


def taxpickup(binned_lists, main_df):
    """

    :param binned_lists:
    :param main_df:
    :return:
    """

    maintaxidlist = [tupla for sublist in
                     binned_lists for tupla in sublist]

    # creating a df for centrifuge classified reads
    if maintaxidlist:
        centrifuged_taxId = pd.DataFrame(maintaxidlist,
                                         columns=['readID', 'Group'])
        main_df = main_df.merge(centrifuged_taxId, on='readID', how='left')
    else:
        main_df = main_df.assign(Group=np.nan)
    return main_df


def dataframe_collection(df, npartitions):
    """

    :param df:
    :param npartitions:
    :return:
    """

    df_size = len(df.index)

    if npartitions > df_size:
        ms = 'Number of partitions is larger '
        ms += 'than dataframe size. The number of'
        ms += 'partitions will be set to: %s' % df_size
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        npartitions = df_size

    indices = np.array_split(np.arange(len(df)), npartitions)
    df_split = [df.iloc[idx] for idx in indices]

    return df_split


def directplast(new_dir_path, file_name, suffix, name, max_plast, cpus):
    """

    :param new_dir_path:
    :param file_name:
    :param suffix:
    :param name:
    :param max_plast:
    :param cpus:
    :return:
    """

    try:
        cmd = 'seqkit fq2fa %s | seqkit split -p %s -j %s -' % (file_name,
                                                                max_plast,
                                                                cpus)
    except:
        cmd = 'seqkit split -1 %s -p %s -j %s ' % (file_name, max_plast, cpus)

    ms = 'directplast cmd is:\n%s' % cmd
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)

    # change name of output and moving it to a different directory
    new_name = '%s%s' % (name, suffix)
    myfiles = relocating_files(new_name, new_dir_path)
    return myfiles


def relocating_files(suffix, new_dir_path):
    """

    :param suffix:
    :param new_dir_path:
    :return:
    """

    splits = os.path.join(os.getcwd(), 'stdin.split')
    files = glob.glob('%s/*' % splits)
    count = 0
    for fname in files:
        count += 1
        new_name = '%s_%s.query' % (suffix, count)
        new_path = os.path.join(new_dir_path, new_name)
        shutil.move(fname, new_path)
    os.system('rm -r %s' % 'stdin.split')
    myfiles = glob.glob('%s/*' % new_dir_path)
    return myfiles


def isfasta(infile):

    cmd = 'head -1 %s' % infile
    ms = 'cmd isfasta is: %s' % cmd
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    line1 = _.stdout.decode('utf-8')
    if line1.startswith('>'):
        return True
    else:
        return False


def split_reads(new_dir_path, reads_to_select, read_file,
                suffix, name, cpus, max_plast):
    """

    :param new_dir_path:
    :param reads_to_select:
    :param read_file:
    :param suffix:
    :param name:
    :param cpus:
    :param max_plast:
    :return:
    """

    path_list = os.path.join(os.getcwd(), '%s_list.tmp' % suffix)
    handle = open(path_list, 'w')
    for read in reads_to_select:
        handle.write('%s\n' % read)
    handle.close()
    _ = isitready(reads_to_select, path_list)
    if isfasta(read_file[0]):
        cmd = 'seqkit grep -f %s %s |' % (path_list, read_file[0])
        cmd += 'seqkit split2 -p %s -j %s' % (max_plast, cpus)
    else:
        if isinstance(read_file, list) and len(read_file) == 1:
            read_file = read_file[0]
        cmd = 'seqkit grep -f %s %s | seqkit ' % (path_list, read_file)
        cmd += 'fq2fa -|seqkit split2 -p %s -j %s' % (max_plast, cpus)

    ms = 'seqkit cmd:\n %s' % cmd
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)

    # change name of output and moving it to a different directory
    # new_name = '%s%s' % (name, suffix)
    myfiles = relocating_files(suffix, new_dir_path)
    return myfiles


def slicing(new_dir_path, main_df, file_name, df_class_report,
            suffix, name, cpus, max_plast):
    """
    It creates a complete data frame containing pre-classified
    columns by readId and NaN objects post centrifuge
    :param new_dir_path:
    :param main_df:
    :param file_name:
    :param df_class_report:
    :param suffix:
    :param name:
    :param cpus:
    :param max_plast:
    :return:
    """

    start_time = time.time()
    if cpus <= 10:
        n_chunks = 10
    else:
        n_chunks = round(cpus * 0.5)  # use few to avoid memory overhead by
        # taxonomy
    collection = dataframe_collection(df_class_report, n_chunks)
    dfs = Parallel(n_jobs=-2, verbose=1)(delayed(
        binningbytaxonomy)(chunk) for chunk in collection)
    #
    elapsed = round(time.time() - start_time, 3)
    minutes = round((elapsed / 60), 3)
    hours = round((minutes / 60), 3)
    ms = 'Elapsed time at taxonomy binning:\n%s seconds or %s ' % (elapsed,
                                                                   minutes)
    ms += 'minutes or %s hours' % hours
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    #
    full_df = taxpickup(dfs, main_df)
    _ = full_df.to_csv('%sfull_df.tsv' % name, sep='\t', header=True)
    slice = full_df[full_df['Group'].isnull()]
    slice = set(slice['readID'].tolist())

    # double - check that not-nulls go to plast
    if slice:
        file_names = split_reads(new_dir_path, slice, file_name, suffix,
                                 name, cpus, max_plast)
        print(file_names, sep=' ', end='\n', file=sys.stdout, flush=True)
        return full_df, file_names
    else:
        m = 'Nothing to slice!'
        print(m, sep=' ', end='\n', file=sys.stdout, flush=True)
        return pd.DataFrame(), m


def Taxonomy_translation(acc2DBpath, pre_acc2tax_df, suffix, typ):
    """
     Executable of acc2tax MUST be in the environmental path
    :param acc2DBpath:
    :param pre_acc2tax_df:
    :param suffix:
    :param typ:
    :return:
    """

    # writing acc2taxIN as input for acc2tax software
    base_fname = string_out(pre_acc2tax_df, suffix, typ)
    ms = 'Extracting taxonomy with acc2tax...'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    line = 'acc2tax -a -d %s -i %sIN%s -n -o %sOUT%s.raw'
    line %= acc2DBpath, base_fname, suffix, base_fname, suffix
    run(line, stdout=PIPE,
        stderr=PIPE, shell=True)
    time.sleep(5)

    taxdump = pd.read_table('%sOUT%s.raw'
                            % (base_fname, suffix), header=None,
                            sep='\t', names=['Accession', 'Taxonomy'],
                            dtype=object)
    # adding readIDs with taxonomy
    if not taxdump.empty:
        taxdump = pre_acc2tax_df.merge(taxdump,
                                       on='Accession', how='left')

        # getting only the domain
        slice = taxdump['Taxonomy'].str.split(',', 2,
                                              expand=True).loc[:, 1]

        # renaming and creating the taxonomic df
        virpattern = r'\b.{0, 20}vir.{0,100}\b'
        viregex = re.compile(virpattern, flags=re.IGNORECASE)
        slice = slice.replace(viregex, value='Vir', regex=True)
        temp = taxdump.merge(pd.DataFrame(slice),
                             left_index=True, right_index=True).rename(
            columns={1: 'Group'})
        temp = temp.reindex(columns=['readID', 'Accession', 'Group'])
        temp = temp.drop_duplicates()
        return temp
    else:
        ms = 'taxdump is empty'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        return taxdump


def customtaxtranslation(seqidtaxid_map, parseby, plast_output):
    """

    :param seqidtaxid_map:
    :param parseby:
    :param plast_output:
    :return:
    """

    # subject ID and first column of seqidtaxid_map must match
    ms = 'Using custom taxonomy ...'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    colnames = ['subject ID', 'taxID', 'Group']
    taxid_map = pd.read_table(seqidtaxid_map, sep='\t',
                              header=None, names=colnames,
                              dtype={0: str, 1: str, 2: str})
    taxids = plast_output.merge(taxid_map, on='subject ID',
                                how='inner')
    # for now make sure that is 'Group'
    taxids = taxids.loc[:, ['readID', parseby]]
    taxids = taxids.groupby('readID', as_index=False).first()

    ms = '\ttaxids in customtaxtranslation: %s\n' %taxids.head()
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

    return taxids


def string_out(df, suffix, typ='RefSeq'):
    """

    :param df:
    :param suffix:
    :param typ:
    :return:
    """

    outname = 'acc2taxIN'
    if (os.path.isfile('%s%s' % (outname, suffix))
            or os.path.isfile('acc2taxOUT%s.raw' % suffix)):
        os.system('rm *acc2tax*')
    if not typ == 'RefSeq':
        string = set(df['Accession'].tolist())
    else:
        try:
            acc = df['subject ID'].str.split('|', 4,
                                             expand=True).loc[:, 3].str.split(
                '.', 1, expand=True).loc[:, 0]
        except:
            acc = df['subject ID']
        string = set(acc.tolist())

    L = len(string)
    string = '\n'.join(string) + '\n'
    with open('%s%s' % (outname, suffix), 'w') as O:
        O.write(string)
    while True:
        time.sleep(3)
        c = run('wc -l %s%s' % (outname, suffix), stderr=PIPE, stdout=PIPE,
                shell=True, universal_newlines=True)
        lc = int(c.stdout.split()[0])
        if lc == L:
            break
    return outname[:-2]


def other_dfs(df):

    lista = ['Archaea', 'Bacteria',
             'Eukaryota', 'Vir']
    xdf = df[~df['Group'].isin(lista)]
    return xdf['readID'].tolist()


def bin(gdf, group, label):
    """

    :param gdf:
    :param group:
    :param label:
    :return:
    """

    try:
        x = gdf.get_group(group)
        _ = x.to_csv('%s_grouped_by_%s.tsv' % (label, group),
                     sep='\t', header=True)
        val = set(gdf.get_group(group).loc[
                  :, 'readID'].tolist())
    except:
        ms = 'empty bin %s %s' % (label, group)
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        val = set()
    return val


def tokeep(grouped_df, group,
           excludegroup, reject_set, label):
    """
    Subsetting by exclusion of 'rejecting set'
    :param grouped_df:
    :param group:
    :param excludegroup:
    :param reject_set:
    :param label:
    :return:
    """

    try:
        eg = [excludegroup]
        df = grouped_df.get_group(group)
        _ = df.to_csv('%s_grouped_by_%s.tsv' % (label, group),
                      sep='\t', header=True)
        ndf = df[~df['Group'].isin(eg)]
        _ = df.to_csv('%s_grouped_by_%s_excl_%s.tsv' % (
            label, group, excludegroup), sep='\t', header=True)
        ndf = set(ndf['readID'].tolist())
        ndf = ndf.difference(reject_set)
        return ndf
    except:
        return set()


def intersects(set1, set2, set3):
    """
    Checking if intersections among bins exist
    :param set1:
    :param set2:
    :param set3:
    :return:
    """

    isects = [set1.intersection(set2),
              set1.intersection(set3),
              set2.intersection(set3)]
    _ = zip(['Bact-Arch', 'Bact-NonCel', 'Arch-Noncel'], isects)
    nr_inters = set([n for s in isects for n in s])
    set1, set2 = [s.difference(nr_inters)
                  for s in [set1, set2]]
    set3 = set3.union(nr_inters)
    return nr_inters, set1, set2, set3


def grouper(label, pref, df, reads):
    """
    Grouping by major categories and check to deduplicate dfs
    :param label:
    :param pref:
    :param df:
    :param reads:
    :return:
    """

    merged, pre_unknown = df
    if not merged.empty:
        _ = merged.to_csv('%s_merged_df.tsv' % label, sep='\t', header=True)
        merged['Group'] = merged['Group'].str.strip()
        grouped_df = merged.groupby('Group', as_index=False)
        Euk = bin(grouped_df, 'Eukaryota', label)
        Bact = tokeep(grouped_df, 'Bacteria', 'Eukaryota', Euk, label)
        Arch = tokeep(grouped_df, 'Archaea', 'Eukaryota', Euk, label)
        Vir = tokeep(grouped_df, 'Vir', 'Eukaryota', Euk, label)

        # making sure no duplicates are present and
        # forcing the merge of eukaryota and unkown.
        intersect, Bact, Arch, Misc = intersects(Bact, Arch, Vir)
        Euk = Euk.difference(intersect)
        Preset = set().union(*[Bact, Arch, Vir, Euk])
        Unk = pre_unknown.union(other_dfs(merged))
        Unk = Unk.difference(Preset)
        print('len unk', len(Unk), sep=' ', end='\n', file=sys.stdout, flush=True)
        EUnk = set(Euk).union(Unk)
        sel_reads = [('Bact', label, pref, Bact, reads),
                     ('Arch', label, pref, Arch, reads),
                     ('EUnk', label, pref, EUnk, reads),
                     ('Misc', label, pref, Misc, reads),
                     ('Unk', label, pref, Unk, reads),
                     ('Euk', label, pref, Euk, reads)]
        return sel_reads
    else:
        Unk = pre_unknown
        sel_reads = [('Unk', label, pref, Unk, reads)]
        return sel_reads


def nullandmerged(df1, df2):
    """

    :param df1:
    :param df2:
    :return:
    """

    if not 'Group' in df1.columns:
        df1 = df1.assign(Group=np.nan)
    df1.Group = df1.Group.astype(str)
    df1 = df1.set_index('readID')
    df2 = df2.set_index('readID')
    combined = df2.combine_first(df1)
    df1 = df1.reset_index()
    combined = combined.reset_index()
    merged = combined.loc[:, ['readID', 'Group']]
    merged = merged.groupby('readID', as_index=False).first()
    merged = merged.loc[merged['Group'].notnull()]
    allreads = set(df1.loc[:, 'readID'].tolist())
    m = set(merged.loc[:, 'readID'].tolist())
    pre_unknown = allreads.difference(m)
    return merged, pre_unknown


def empty(df):
    """

    :param df:
    :return:
    """

    if not 'Group' in df.columns:
        df = df.assign(Group=np.nan)
    pre_unknown = set(df['readID'].tolist())
    merged = pd.DataFrame()
    ms = 'In emtpy: pre_unknown\n %s, merged\n%s' % (pre_unknown, merged)
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    return merged, pre_unknown


def input_check_and_setup(user_args):
    """

    :param user_args: dictionary containing arguments
    :return: type: list of the checked arguments
    """

    # Check that plast is in path
    if not program_exists('plast'):
        m = 'Eukfinder requirement not met: plast is not installed!\nExiting program...'
        print(m, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)


    # --- Global required values ---
    map_id_path = user_args['plast_id_map']
    plast_path = user_args['plast_database']
    #acc2tax_path = user_args['acc2tax_database']
    pid = user_args['pid']
    cov = user_args['cov']
    e_value = user_args['e_value']
    tax_up = True if user_args['taxonomy_update'] == 'True' else False
    n_cpu = user_args['number_of_threads']
    max_plast = user_args['number_of_chunks']
    base_outname = user_args['out_name']

    # args specific to illumina workflow
    if user_args['func'].__name__ == 'long_seqs':
        reads = [user_args['long_seqs']]
        classification = ['Pass']
    else:
        reads = [user_args['r1'], user_args['r2'], user_args['un']]
        if type(user_args['pclass']) is not None and type(
                user_args['uclass']) is not None:
            classification = [user_args['pclass'],
                              user_args['uclass']]
        else:
            classification = ['None', 'None']

    # check file existence and classification
    redef_reads = [it if (os.path.exists(it)
                          and os.stat(it).st_size != 0)
                   else 'None' for it in reads]

    reads_paths = [os.path.abspath(it) if os.path.exists(it)
                   else 'None' for it in redef_reads]

    redef_class = [os.path.abspath(it) if (os.path.exists(it)
                                           and os.stat(
                it).st_size != 0) else 'None' for it in
                   classification]

    if not user_args['func'].__name__ == 'long_seqs':
        # quick file format check up
        boolean_seqs = [validate_input(f, 'fastq')
                        for f in reads_paths]
        if False in boolean_seqs:
            sys.exit(0)
        # quick file format check up for centrifuge
        boolean_centrif = [validate_input(f, 'centrifuge')
                           for f in redef_class]
        if False in boolean_centrif:
            sys.exit(0)
        declared = reads_paths + redef_class
    else:
        boolean_centrif = [validate_input(f, 'fastx') for f in reads_paths]
        if False in boolean_centrif:
            sys.exit(0)
        print(f'reads_paths is {reads_paths}', sep=' ', end='\n',
              file=sys.stdout, flush=True)
        declared = reads_paths

    if any(i == 'None' for i in declared):  # double check
        ms = '\nIt seems that at least one declared file DOES NOT exist '
        ms += 'in the directory and has been labeled as "None". \n'
        ms += 'Please check your command line.\nDeclared files are:\n'
        ms += '*** Reads and Classification files are : ***\n'
        ms += '%s\n' % '\n'.join(declared)
        ms += '---  Exiting program   ---'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)

    if not os.path.exists(plast_path):
        ms = '\n\nPath to plast database does not exist. Terminating program'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)

    if not os.path.exists(map_id_path):
        format_check = plast_path
        ms = '\n\nPath to map for plast database does not exist.\n'
        ms += 'Terminating program'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)

    #acc2tax_fn = ['nodes.dmp', 'names.dmp', 'acc2tax_prot_all.txt',
    #              'acc2tax_nucl_all.txt']
    #for fn in acc2tax_fn:
    #    apath = os.path.join(, fn)
    #    if not os.path.exists(apath):
    #        ms = '\n\nSomething is wrong with the acc2tax database.\n'
    #        ms += 'Please make sure you specify the path to a '
    #        ms += 'directory that contains the following files: '
    #        ms += 'nodes.dmp, names.dmp, acc2tax_prot_all.txt, '
    #        ms += 'acc2tax_nucl_all.txt\nTerminating program.'
    #        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    #        sys.exit(0)

    if e_value is None:
        e_value = 1.0
    if pid is None:
        pid = 10.0
    if cov is None:
        cov = 0.0

    if len(glob.glob('%s.*.cf' % user_args['cdb'])) != 4:
        cdb = '\nCentrifuge database cannot be found!\n'
        cdb += 'Exiting program...\n'
        print(cdb, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)
    cdb_path = user_args['cdb']
    minhit_len = user_args['mhlen']

    if not user_args['func'].__name__ == 'long_seqs':
        #aplast_path = user_args['ancillary_plast_database']
        #amap_id_path = user_args['ancillary_plast_id_map']
        # max_len = user_args['mlen']
        max_mem = user_args['max_m']
        # kmers
        if user_args['kmers'] is None or user_args['kmers'] == 'None':
            k = '21,33,55'
        else:
            k = suitable_kmers(user_args['kmers'])

        # cdb_path = user_args['cdb']
        # if aplast_path is None or not os.path.exists(aplast_path):
        #     aplast_path = plast_path
        #     amap_id_path = map_id_path
        #     ms = '\nNo ancillary plast database found. Current plast_database '
        #     ms += f'{plast_path} will be used for re-classification\n'
        #     print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

        # if os.path.exists(aplast_path):
        #     if not os.path.exists(amap_id_path):
        #         ms = '\n\nNo map for plast_database was found. Exiting program'
        #         print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        #         sys.exit(0)

        # if len(glob.glob('%s.*.cf' % user_args['cdb'])) != 4:
        #     cdb = '\nCentrifuge database does not exist in path.\n'
        #     cdb += 'Exiting program.\n'
        #     print(cdb, sep=' ', end='\n', file=sys.stdout, flush=True)
        #     sys.exit(0)

        arguments = [tax_up, reads_paths, redef_reads, redef_class,
                     plast_path, map_id_path, n_cpu,
                     max_plast, e_value, pid, user_args['func'].__name__,
                     base_outname, max_mem, cdb_path,
                     cov, minhit_len, k]  #aplast_path, amap_id_path,
    else:
        arguments = [tax_up, reads_paths, redef_reads, plast_path,
                     map_id_path, n_cpu, max_plast,
                     e_value, pid, base_outname, cov,  cdb_path, minhit_len]

    return arguments


def head2(infile):

    o = run(f'head -2 {infile}', stdout=PIPE, stderr=PIPE, shell=True)
    return o.stdout.decode('utf-8').strip('\n').split('\n')


def validator(lines_list, allowed_bases, s, ms):

    if len(lines_list) == 2:
        header, seq = lines_list[0], lines_list[1]
        if header.startswith(s):
            try:
                # check only the 10 first bases (trimmomatic keeps at least 40)
                # adapters file should contain sequences of size >= 19n
                subseq = seq[:19].upper()
                not_allowed = [base for base in subseq if base not in
                               allowed_bases]
                if not_allowed:
                    v = ' '.join(not_allowed)
                    ms = 'At least one base in a sequence is not allowed.\n'
                    ms += 'Offender base(s) is(are): %s' % v
                    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
                    return False
                return True
            except IndexError:
                ms = 'At leas one sequence is smaller than 19 bases'
                print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
                return False
        else:
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            return False
    else:
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        return False


def fastaq_validator(lines_list, infile, seqtype):

    # Recogize file type of very large files relaying on the first sequence.
    # it does not consider sequence duplications or errors in sequences due
    # to file size
    allowed_bases = ['A', 'C', 'G', 'T', 'N', 'X']
    ms = f'Invalid file type. Expected {seqtype} for {infile}'
    if seqtype == 'fastx':
        validation = validator(lines_list, allowed_bases, '@', ms)
        if validation is False:
            validation = validator(lines_list, allowed_bases, '>', ms)
        return validation
    else:
        s = '>' if seqtype == 'fasta' else '@'
        validation = validator(lines_list, allowed_bases, s, ms)
        return validation


def report_validator(lines_list, infile, seqtype):

    expected = ['readID', 'seqID', 'taxID', 'score', '2ndBestScore',
                'hitLength', 'queryLength', 'numMatches']
    if len(lines_list) == 2:
        header, value = lines_list[0], lines_list[1]
        headers = header.split('\t')
        intersect = set(headers).intersection(set(expected))
        if len(intersect) != len(expected):
            ms = ('File does not appear to be a centrifuge report. Expected'
                  f'headers are: {headers}\n')
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            return False
        else:
            return True


def validate_input(infile, seqtype):

    print(f'Quick input format check for {infile}',
          end='\n', file=sys.stdout, flush=True)
    ms = 'Exiting program'
    two_first_lines = head2(infile)
    if infile == 'None':
        return False
    if seqtype == 'fasta' or seqtype == 'fastq':
        fastaq = fastaq_validator(two_first_lines, infile, seqtype)
        if fastaq is False:
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            return False
        return True
    if seqtype == 'fastx':
        fastq = fastaq_validator(two_first_lines, infile, 'fastq')
        if fastq is False:
            # check if file is fasta format
            fasta = fastaq_validator(two_first_lines, infile, 'fasta')
            if fasta is False:
                m = 'File is not fasta format. This workflow requires either '
                m += 'fastq or fasta format files.'
                m += ms
                print(m, sep=' ', end='\n', file=sys.stdout, flush=True)
                return False
            else:
                print('File is fasta format', sep=' ', end='\n',
                      file=sys.stdout, flush=True)
                return True
        else:
            print('File is fastq format', sep=' ', end='\n', file=sys.stdout,
              flush=True)
            return True
    if seqtype == 'centrifuge':
        centrif = report_validator(two_first_lines, infile, seqtype)
        if centrif is False:
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            return False
        return True


def suitable_kmers(k):

    newk = k.replace(' ', '')
    if ',' not in newk:
        ms = 'A single kmer appears to have been specified: %s. ' % k
        ms += 'Please specify at least 3 different \n'
        ms += 'odd kmers under 128. e.g. 21,33,55\n'
        ms += 'Terminating program.'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)
    else:
        v = [int(i) for i in newk.split(',')]
        ktype = {i: (i % 2) == 0 for i in v}
        vs = [i for i in ktype.keys() if ktype[i] is False]
        u = [i for i in ktype.keys() if i < 128]
        if len(vs) != len(v) or len(u) != len(v):
            ms = 'Specified kmers are %s. %s of these ' % (len(v), len(vs))
            ms += 'are odd and %s are under 128' % len(u)
            ms += 'kmer requirements for metaspades have not been met.\n'
            ms += 'Terminating program.'
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            sys.exit(0)
        return newk


def bytearray_sets(infile):
    """

    :param infile:
    :return:
    """

    start_time = time.time()
    mycmd = 'seqkit seq -n %s' % infile
    _ = run(mycmd, stderr=PIPE, stdout=PIPE, shell=True)
    headers = bytearray(_.stdout)
    headers = headers.rstrip().decode('utf-8').split('\n')
    elapsed = round(time.time() - start_time, 3)
    minutes = round((elapsed / 60), 3)
    hours = round((minutes / 60), 3)
    ms = 'Elapsed time at bytearrays:\n%s seconds or %s ' % (elapsed, minutes)
    ms += 'minutes or %s hours' % hours
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    return headers


def mkdir(dir_name):
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    else:
        pass


def cpu_chunk(cpus, maxplast):
    """

    :param cpus:
    :param maxplast:
    :return:
    """

    cpus = int(cpus)
    max_plast = int(maxplast)
    if cpus < max_plast:
        line = 'Not enough CPUs have been allocated... '
        line += 'Terminating program with status 1'
        print(line, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(1)
    else:
        max_cpu = int(cpus / max_plast)
        return max_cpu


def getting_prefix(file_path):
    """

    :param file_paths:
    :return:
    """

    fname = file_path.split(os.sep)[-1]
    fname = fname.lower()
    prefix = [None if fname is None else
              re.split(r'.fastq$|.fq$', fname)[0]][0]
    return prefix


def rename_reads(readfile):
    """

    :param readfile:
    :return:
    """

    read, orient = readfile
    original_path = os.path.abspath(read)
    ftype = re.split(r'.fastq$|.fq$', read.lower())
    if len(ftype) == 1:
        tail = 'fasta'
    else:
        tail = 'fastq'
    newrfile_path = os.path.join(os.getcwd(), 'tmp.%s.%s' % (orient, tail))
    cmd = r'seqkit replace -p "\s.+" %s -o %s' % (original_path, newrfile_path)
    ms = 'rename reads cmd:\n%s' % cmd
    _ = run(cmd, stderr=PIPE, stdout=PIPE, shell=True)
    return newrfile_path, ms


def pair_read_handler(cpus, max_files, files_paths, reports_list,
                      readfile_list, stamp, base_outname):
    """

    :param cpus:
    :param max_files:
    :param files_paths:
    :param reports_list:
    :param readfile_list:
    :param stamp:
    :return:
    """

    ms = 'Pair read handler..'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    reports = [reading_reports(inf) for inf in reports_list]
    reads = [bytearray_sets(readfile) for readfile in readfile_list]
    dirname = 'tmps_%s_%s' % (base_outname, stamp)
    mkdir(dirname)
    dirname_path = os.path.abspath(dirname)
    ms = 'working paths are:\n', files_paths
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    path_r1, path_r2, path_ur = files_paths
    pair_files = [files_paths[0], files_paths[1]]
    unpair_file = [files_paths[2]]

    # main data frame by read name
    df_r1, df_r2, df_ur = [pd.DataFrame(read_list,
                                        columns=['readID'])
                           if bool(read_list) else pd.DataFrame()
                           for read_list in reads]  # read data frames

    dfr_r1r2, dfr_ur = reports  # reports data frames
    ##
    os.chdir(dirname)
    existing_dfs = []
    #  pair reads setup
    # base = pair_prefixes[0]
    maindf_r1r2 = df_r1.merge(df_r2, on='readID', how='inner')
    if not df_r1.empty and not df_r2.empty:
        # merge main r1 and r2 df to make sure only pair reads
        # are left in one single df
        if not dfr_r1r2.empty:
            _ = dfr_r1r2.to_csv('p_report_r1r2.tsv', sep='\t',
                                header=True)

            # creating chunks for only one set of the pair reads
            # selecting only r1
            pfull_df, f_names = slicing(dirname_path, maindf_r1r2,
                                        path_r1, dfr_r1r2, 'pPlastSearch',
                                        'p', cpus, max_files)
            existing_dfs.append(('Pair', base_outname,
                                 pfull_df, pair_files, f_names))
        else:

            # direct plast searches on r1 with no apriori centrifuge info
            # the dataframe returned will be maindf_r1r2 no taxonomy apriori
            f_names = directplast(dirname_path, path_r1, 'pPlastSearch',
                                  'p', max_files, cpus)
            existing_dfs.append(('Pair', base_outname, maindf_r1r2,
                                 pair_files, f_names))
    #  unpair reads must be provided as a concatenated file

    # base = unpair_prefix
    if not df_ur.empty:
        if not dfr_ur.empty:
            _ = dfr_ur.to_csv('uReportNR.tsv', sep='\t', header=True)
            ufull_df, f_names = slicing(dirname_path, df_ur,
                                        path_ur, dfr_ur, 'uPlastSearch',
                                        'u', cpus, max_files)
            existing_dfs.append(('Up', base_outname,
                                 ufull_df, unpair_file, f_names))
        else:
            # direct plast searches with no apriori centrifuge info
            f_names = directplast(dirname_path, path_ur, 'uPlastSearch',
                                  'Up', max_files, cpus)
            existing_dfs.append(('Up', base_outname, df_ur, unpair_file,
                                 f_names))

    return existing_dfs


def single_read_handler(cpus, max_files, files_paths, reports_list,
                        readfile_list, stamp, base_outname):
    """

    :param cpus:
    :param max_files:
    :param files_paths:
    :param reports_list:
    :param readfile_list:
    :param stamp:
    :return:
    """

    reports = [reading_reports(inf) for inf in reports_list]
    reads = [bytearray_sets(readfile) for readfile in readfile_list]
    dirname = 'tmps_%s_%s' % (base_outname, stamp)
    mkdir(dirname)
    dirname_path = os.path.abspath(dirname)
    path1 = files_paths
    df_un = [pd.DataFrame(read_list, columns=['readID'])
             if bool(read_list) else pd.DataFrame()
             for read_list in reads][0]  # reads data frames
    dfr_un = reports[0]  # reports data frames
    os.chdir(dirname)
    existing_dfs = []
    if not df_un.empty:
        jobs = unpair_dfs(dirname_path, files_paths, df_un, path1, dfr_un,
                          cpus, max_files, base_outname)
        existing_dfs.append(jobs)
    return existing_dfs


def unpair_dfs(dirname_path, files_paths, df_un, path1, dfr_un,
               cpus, max_files, base_outname, name='u'):
    """

    :param dirname_path:
    :param files_paths:
    :param df_un:
    :param path1:
    :param dfr_un:
    :param cpus:
    :param max_files:
    :param name:
    :return:
    """

    #
    if not dfr_un.empty:
        _ = dfr_un.to_csv('uReportNR.tsv', sep='\t', header=True)
        ufull_df, f_names = slicing(dirname_path, df_un,
                                    path1, dfr_un, 'uPlastSearch',
                                    name, cpus, max_files)
        record = ('Up', base_outname, ufull_df, files_paths, f_names)
    else:
        # direct plast searches with no apriori centrifuge info
        f_names = directplast(dirname_path, files_paths, 'uPlastSearch',
                              '', max_files, cpus)
        record = ('Up', base_outname, df_un, files_paths, f_names)

    return record


def plast_search_n_results(CPUs, pDBpath, working_dfs, evalue, pid, cov, lr):
    """

    :param CPUs:
    :param pDBpath:
    :param working_dfs:
    :param evalue:
    :return:
    """


    # multi-threaded plast searches
    start_time = time.time()
    tmps_outputs = []
    for entry in working_dfs:
        pattern = '%s*.query' % entry[0][0].lower()
        tmp_outputs = plastsearches(CPUs, pDBpath, evalue, pattern)
        # if os.stat(tmp_outputs).st_size != 0:
        tmps_outputs.append(tmp_outputs)

    # parsing plast outputs
    plast_outputs = []
    for entry in tmps_outputs:
        out = parsingplastsearches(entry, pid, cov, lr)
        plast_outputs.extend(out)

    # checking for taxonomy to dataframes
    elapsed = round(time.time() - start_time, 3)
    minutes = round((elapsed / 60), 3)
    ms = 'Full plast searches took:\n'
    ms += '%s seconds or %s minutes' % (elapsed, minutes)
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    return plast_outputs


def create_df_taxonomy(plast_outputs_list, working_dfs,
                       pid, seqidtaxid_map=None): #taxDBpath,
    """

    :param plast_outputs_list:
    :param taxDBpath:
    :param working_dfs:
    :param percentid:
    :param seqidtaxid_map:
    :return:
    """

    cladfs = []
    for entry in plast_outputs_list:
        etk, plast_output = entry
        if not plast_output.empty:
            outfile_name = f'Parsed_plast_output_{etk.lower()}.plout.tsv'
            if not os.path.isfile(outfile_name):
                _ = plast_output.to_csv(outfile_name, sep='\t', header=True)
            if seqidtaxid_map is not None:
                dfacc2tax = customtaxtranslation(seqidtaxid_map, 'Group',  plast_output)

                if dfacc2tax.empty:
                    ms = 'Something seems wrong with the acc2tax information '
                    ms += 'provided'
                    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            #else:
            #    dfacc2tax = Taxonomy_translation(taxDBpath, plast_output, etk,
            #                                     'RefSeq')

            # --- After plast searches -------------
            # Merging all dataframes by readId using taxonomic info
            # this will fail if there is nothing to merge as it'll generate an
            # empty df

            print("\n working_dfs IN create_df_taxonomy before tmpcladfs,\n", working_dfs)

            tmpcladfs = [(label, filename, nullandmerged(df, dfacc2tax),
                          read_paths) for label, filename, df, read_paths,  other in working_dfs if
                         (not dfacc2tax.empty and
                          label.startswith(etk))]
            cladfs.extend(tmpcladfs)
        else:
            ## re-check outputs... is there a bug here?
            ms = 'plastoutput empty test for unpair'
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            tmpcladfs = [(label, filename, empty(df), read_paths)
                         for label, filename, df, read_paths, other in
                         working_dfs if label.startswith(etk)]
            ms = 'tmpcladfs\n %s' % tmpcladfs
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
            cladfs.extend(tmpcladfs)
    return cladfs


def writinggroups(cladfs, dir_path, boutname):
    """

    :param cladfs:
    :return:
    """

    pre_args = []
    for label, prefix, df, read_paths in cladfs:
        pre_args.extend(grouper(label, prefix, df, read_paths))
    ms = '****\tWriting files by group\t****\n'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    for arg in pre_args:
        _ = writer(arg, dir_path, boutname)
    return 'Done'


def isitready(sel_reads, fpath):
    """

    :param sel_reads:
    :param fpath:
    :return:
    """

    while True:
        lcount = run('wc -l %s' % fpath,
                     stderr=PIPE, stdout=PIPE, shell=True)
        lines = lcount.stdout.decode('utf-8')
        lines = int(lines.split()[0])
        if lines == len(sel_reads):
            break
        time.sleep(1)
    return 'Finished'


def writer(args, dir_path, boutname):
    """

    :param args:
    :return:
    """

    group, label, fname, read_set, path = args

    if len(read_set) > 0:
        path_list = os.path.join(os.getcwd(), '%s_%s_list.tmp' % (group, label))
        start_time = time.time()
        handle = open(path_list, 'w')
        for read in read_set:
            handle.write('%s\n' % read)
        handle.close()
        _ = isitready(read_set, path_list)
        elapsed = round(time.time() - start_time, 3)
        minutes = round((elapsed / 60), 3)
        ms = 'Elapsed time for writing accessions:'
        ms += '%s seconds; %s minutes' % (elapsed, minutes)
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

        start_time = time.time()
        # patch_path0 = os.path.join(path[0], 'Eukfinder_Results')
        # patch_path1 = os.path.join(path[1], 'Eukfinder_Results')

        # use dir_path
        if label == 'Pair':
            outs = [('%s.%s.R1.fq' % (fname, group), path[0]),
                    ('%s.%s.R2.fq' % (fname, group), path[1])]
            for entry in outs:
                outname, path = entry
                out_path = os.path.join(dir_path, 'Classified_reads')
                abs_path_outname = os.path.join(out_path, outname)
                cmd = 'seqkit grep -f %s -o %s %s ' % (path_list,
                                                       abs_path_outname,
                                                       path)
                ms = 'cmd is %s\n' % cmd
                _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
                if _.stdout.decode('utf-8').strip('\n') != '':
                    ms = 'stdout: %s\n' % _.stdout.decode('utf-8')
                elapsed = round(time.time() - start_time, 3)
                minutes = round((elapsed / 60), 3)
                ms += 'Elapsed time for writing Paired outs:'
                ms += '%s seconds; %s minutes' % (elapsed, minutes)
                print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

        # unpaired workflow
        else:
            if fname.startswith('scf_'):
                fname = fname.split('.fasta')
                outname = '%s.%s.fasta' % (fname[0], group) ###1
            else:
                ms = 'output is not list or scaffold. The fname is %s' % fname
                print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
                if isfasta(path):
                    outname = 'Eukfinder_results/%s.%s.fasta' % (fname, group)
                else:
                    outname = 'Classified_reads/%s.%s.un.fq' % (fname, group)

            outname_path = os.path.join(dir_path, outname)

            if 'scf_' in outname_path:
                pathched_path = os.path.join(dir_path, 'Eukfinder_results')
                outname_path = os.path.join(pathched_path, outname)

            if isinstance(path, list):
                cmd = 'seqkit grep -f %s -o %s %s ' % (path_list, outname_path,
                                                       path[0])
            else:
                cmd = 'seqkit grep -f %s -o %s %s ' % (path_list, outname_path,
                                                       path)

            ms = 'scf_ seqkit grep cmd:\n%s\n' % cmd
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
            ms = 'stdout: %s\n' % _.stdout.decode('utf-8')
            elapsed = round(time.time() - start_time, 3)
            minutes = round((elapsed / 60), 3)
            ms += 'Elapsed time for writing Unpaired outs:'
            ms += '%s seconds; %s minutes' % (elapsed, minutes)
            # print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        # new_path = os.path.join(os.getcwd(), 'tmps_%s_%s' % (boutname, stamp))
        # shutil.move(path_list, new_path)
    else:
        m = f'Read_set is empty. Group: {group}, label: {label}, path: {path}, fname: {fname}'
        print(m, sep=' ', end='\n', file=sys.stdout, flush=True)
    return 'Done'


# --- post classification ---

def standard_out(program, outinfo):

    # getting standard output
    o1 = outinfo.stdout.decode('utf-8').strip('\n')
    o2 = outinfo.stderr.decode('utf-8').strip('\n')
    sms = f'{program} stdout is:\n{o1}\n'
    sms += f'{program} stderr is\n{o2}\n'
    print(sms, sep=' ', end='\n', file=sys.stdout, flush=True)
    return 'Done'


def assembly(read_tuple, basename, threads, maxmem, k):

    ms = 'Starting assembly phase. This will take a while...'
    print(ms, sep='\t', end='\n', file=sys.stdout, flush=True)
    R1, R2, uR1R2 = read_tuple
    if os.path.isfile(uR1R2):
        e = '--pe1-s %s' % uR1R2
    else:
        e = ''
    outdir = '%s_metaspades_out' % basename
    cmd = 'metaspades.py -t %s -m %s ' % (threads, maxmem)
    cmd += '--pe1-1 %s --pe1-2 %s %s ' % (R1, R2, e)
    cmd += '-o %s -k %s' % (outdir, k)
    ms = 'Metaspades cmd_line:\n %s' % cmd
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

    scafout = run(cmd, stderr=PIPE, stdout=PIPE, shell=True)

    # getting standard output
    _ = standard_out('Metaspades', scafout)
    # checking assembly
    p1 = os.path.join(os.getcwd(), outdir)
    assembly_path = os.path.join(p1, 'scaffolds.fasta')
    # make new dir
    ccc = 'Centrifuge_contig_classification'
    mkdir(ccc)
    new_dir = os.path.join(os.getcwd(), ccc)
    if os.path.isfile(assembly_path):
        # making a copy of the assembly in the new dir
        shutil.copy(assembly_path, new_dir)
    else:
        # cautionary waiting time
        time.sleep(60)
        # second check
        if os.path.isfile(assembly_path):
            # making a copy of the assembly in new dir
            shutil.copy(assembly_path, new_dir)
        else:
            ms = "The assembly process failed. "
            ms += "There might not be enough reads left that overlap"
            ms += "to create and elongate contigs.\n"
            ms += "Terminating program.\n"
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            sys.exit(0)
    # rename assembly using basename
    scaffolds = change_c_names(new_dir, 'scaffolds.fasta', basename)
    ms = 'Newly named assembly is:\n %s' % scaffolds
    print(ms, sep='\t', end='\n', file=sys.stdout, flush=True)
    return scaffolds, 'scf_%s.fasta' % basename


def change_c_names(new_dir, assembl, basename):

    assembly = os.path.join(new_dir, assembl)
    new = {}
    track = {}
    with open(assembly) as I:
        records = I.read().split('>')[1:]
        count = 0
        for entry in records:
            count += 1
            value = entry.split('\n')
            name = 'scaffold_%s' % count
            if name not in new:
                new[name] = '\n'.join(value[1:])
                track[name] = value[0]
            else:
                ms = 'repeated accession %s' % value[0]
                print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

    outname1 = os.path.join(new_dir, 'scf_%s.fasta' % basename)
    outname2 = os.path.join(new_dir, 'scf_%s_cross_check.txt' % basename)
    with open(outname1, 'w') as O1, \
            open(outname2, 'w') as O2:
        for key in new.keys():
            fasta = '>%s\n%s' % (key, new[key])
            check = '%s\t%s\n' % (key, track[key])
            O1.write(fasta)
            O2.write(check)
    os.remove(assembly)
    return outname1


def post_assembly(threads, out_name, fasta_path, dbpath, mhlen):

    os.chdir('Centrifuge_contig_classification')
    cline = 'centrifuge -f --threads %s -k 1 ' % threads
    cline += '--min-hitlen %s -x %s ' % (mhlen, dbpath)
    cline += '-U %s -S %s_scf.centrifuge_UP ' % (fasta_path, out_name)
    cline += '--report-file %s_scf.centrifuge_UP.tsv ' % out_name
    print('Post assembly classification cmd:', cline)

    centout = run(cline, stdout=PIPE, stderr=PIPE, shell=True)
    report_name = os.path.abspath('%s_scf.centrifuge_UP' % out_name)
    valid, ms = validate_output(centout)
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    if not valid:
        sys.exit(1)
    os.chdir('..')
    return os.path.abspath(report_name)


def name_check(infile):

    cmd = 'head -1 %s' % infile
    head = run(cmd, stderr=PIPE, stdout=PIPE, shell=True)
    myheader = head.stdout.decode('utf-8').strip('\n')
    boolean = True if ' ' in myheader else False
    return boolean


# --- Execute ---

def mini_loop(user_args_lst, ftype):

    for file in user_args_lst:
        if os.path.exists(file) and os.stat(file).st_size > 0:
            valid = validate_input(os.path.abspath(file), ftype)
            if not valid:
                sys.exit(0)
        else:
            ms = f"Declared file {file} does not exist in path or is empty\n"
            ms += "Terminating program."
            print(ms, sep='\t', end='\n', file=sys.stdout, flush=True)
            sys.exit(0)
    return 'Done'


def perform_prep(user_args):
    """
    Full argument set
    :param user_args: arguments for trimmomatic, bowtie2, and centrifuge
    :return: all arguments to trim, map and centrifuge
    """

    start_time = time.time()
    stamp = time.strftime('%Y%m%d', time.gmtime(start_time))
    tmpdirname = 'tmp_readprep_%s%s' % (user_args['out_name'], stamp)
    log = os.path.abspath('Read_prep_%s.log' % stamp)
    sys.stdout = open(log, 'w')

    # Trimming raw reads
    ms = 'Run has started with arguments:\n'
    for key in user_args:
        ms += '%s: %s, ' % (key, user_args[key])
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

    cen = glob.glob('%s*' % user_args['cdb'])
    if len(cen) != 4:
        ms = '\nThere is something wrong with the centrifuge database '
        ms += 'declared. \nMake sure that the '
        ms += 'path declared specifies the shared suffix of '
        ms += 'the database files (see user manual)\n. Terminating program'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)
    plastdb_home = glob.glob('%s*' % user_args['cdb'])

    # fast check for file format
    _ = mini_loop([user_args['r1'], user_args['r2']], 'fastq')
    _ = mini_loop([user_args['hg'], user_args['illumina_clip']],
                  'fasta')

    # declare file existence
    # abspaths
    adapters = os.path.abspath(user_args['illumina_clip'])
    oR1 = os.path.abspath(user_args['r1'])
    oR2 = os.path.abspath(user_args['r2'])
    ori_host_genome = os.path.abspath(user_args['hg'])
    #
    mss = 'Eukfinder v%s is using python %s\n' % (__version__,
                                                 platform.python_version())
    mss += 'Preparing reads for analysis...'

    print(mss, sep=' ', end='\n', file=sys.stdout, flush=True)
    mkdir(tmpdirname)
    tmpdirname_path = os.path.abspath(tmpdirname)
    os.chdir(tmpdirname_path)
    shutil.copy(ori_host_genome, os.getcwd())

    readfile_list = [(oR1, 'R1'), (oR2, 'R2')]
    xR1, xR2 = Parallel(n_jobs=-2)(delayed(rename_reads)(readfile)
                            for readfile in readfile_list)
    nR1, m1 = xR1
    nR2, m2 = xR2
    ms += m1 + m2
    ms += '\nReads have been renamed\n'

    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    R1, R2, uR1R2 = trimming(user_args['out_name'], nR1, nR2, adapters,
                             user_args['wsize'], user_args['qscore'],
                             user_args['hcrop'], user_args['mlen'],
                             user_args['threads'], user_args['leading_trim'],
                             user_args['trail_trim'],user_args['qenc'])
    #
    # Mapping host out
    ms = 'Getting rid of host reads...'
    print(ms, sep='\t', end='\n', file=sys.stdout, flush=True)

    # validate trimmomatic output
    _ = mini_loop([R1, R2], 'fastq')

    b2_index = bowtie2build(user_args['hg'])
    outr1, outr2, outr1r2, cmd_line = bowtie2(user_args['out_name'],
                                              user_args['threads'], R1, R2,
                                              uR1R2, b2_index, 'fastq-host')

    ms = 'bowtie2 cmd line is: %s\noutr1: %s\noutr2: %s\noutr1r2: %s\n'
    ms %= cmd_line, outr1, outr2, outr1r2

    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    _ = run(cmd_line, stdout=PIPE, stderr=PIPE, shell=True)

    # validate bowtie2 output
    _ = mini_loop([outr1, outr2, outr1r2], 'fastq')

    # centrifuge host-less reads
    ccmd_pline, p_report = centrifuge(user_args['out_name'], 
                                      (outr1, outr2),
                                      user_args['threads'],
                                      user_args['mhlen'],
                                      user_args['cdb'], 1, pair=True,
                                      fastq=True)

    pcent_out = run(ccmd_pline, stdout=PIPE, stderr=PIPE, shell=True)
    valid, ms = validate_output(pcent_out)
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    if not valid:
        sys.exit(1)
    ccmd_upline, up_report = centrifuge(user_args['out_name'], outr1r2,
                                        user_args['threads'],
                                        user_args['mhlen'],
                                        user_args['cdb'], 1, pair=False,
                                        fastq=True)
    ucent_out = run(ccmd_upline, stdout=PIPE, stderr=PIPE, shell=True)
    valid, ms = validate_output(ucent_out)
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    if not valid:
        sys.exit(1)

    # ---- Deleting temporary files ----
    os.chdir('..')
    my_preps = (outr1, outr2, outr1r2, p_report, up_report)
    new_abspaths = []
    ms = 'Something went wrong.\n'
    skip = []
    for f in my_preps:
        try:
            myf = os.path.split(f)[1]
            print(f'Moving {f} to {os.getcwd()}', end='\n', file=sys.stdout,
                  flush=True)
            nflocation = os.path.join(os.getcwd(), myf)
            if os.path.isfile(nflocation):
                m = f'WARNING:\nFile {myf} already exists in {os.getcwd()} '
                m += 'and will not be overwritten.\nNewly processed file will '
                m += 'remain in the tmp directory'
                print(m, sep=' ', end='\n', file=sys.stdout, flush=True)
                skip.append(f)
            else:
                shutil.move(f, os.getcwd())
                new_abspaths.append(nflocation)
        except:
            ms += f'ERROR: attempt to move {f} to {os.getcwd()} failed\n'
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            sys.exit(1)
    if not skip:
        try:
            os.system('rm -r %s' % tmpdirname_path)
        except:
            ms += '%s does not seem to exist\n'
            ms += 'Terminating program.'
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            sys.exit(1)
        ms = "Files are ready to use with 'short_seqs' mode:\n%s" % '\n'.join(
            new_abspaths)
        print(ms, flush=True)
    return my_preps

def perform_prep_env(user_args):
    """
    Full argument set
    :param user_args: arguments for trimmomatic and centrifuge
    :return: all arguments to trim and centrifuge
    """

    start_time = time.time()
    stamp = time.strftime('%Y%m%d', time.gmtime(start_time))
    tmpdirname = 'tmp_readprep_%s%s' % (user_args['out_name'], stamp)
    log = os.path.abspath('Read_prep_%s.log' % stamp)
    sys.stdout = open(log, 'w')

    # Trimming raw reads
    ms = 'Run has started with arguments:\n'
    for key in user_args:
        ms += '%s: %s, ' % (key, user_args[key])
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

    cen = glob.glob('%s*' % user_args['cdb'])
    if len(cen) != 4:
        ms = '\nThere is something wrong with the centrifuge database '
        ms += 'declared. \nMake sure that the '
        ms += 'path declared specifies the shared suffix of '
        ms += 'the database files (see user manual)\n. Terminating program'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)
    plastdb_home = glob.glob('%s*' % user_args['cdb'])

    # fast check for file format
    _ = mini_loop([user_args['r1'], user_args['r2']], 'fastq')
    _ = mini_loop([user_args['illumina_clip']], 'fasta')

    # declare file existence
    # abspaths
    adapters = os.path.abspath(user_args['illumina_clip'])
    oR1 = os.path.abspath(user_args['r1'])
    oR2 = os.path.abspath(user_args['r2'])

    #
    mss = 'Eukfinder v%s is using python %s\n' % (__version__,
                                                 platform.python_version())
    mss += 'Preparing reads for analysis...'

    print(mss, sep=' ', end='\n', file=sys.stdout, flush=True)
    mkdir(tmpdirname)
    tmpdirname_path = os.path.abspath(tmpdirname)
    os.chdir(tmpdirname_path)

    readfile_list = [(oR1, 'R1'), (oR2, 'R2')]
    xR1, xR2 = Parallel(n_jobs=-2)(delayed(rename_reads)(readfile)
                            for readfile in readfile_list)
    nR1, m1 = xR1
    nR2, m2 = xR2
    ms += m1 + m2
    ms += '\nReads have been renamed\n'

    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    R1, R2, uR1R2 = trimming(user_args['out_name'], nR1, nR2, adapters,
                             user_args['wsize'], user_args['qscore'],
                             user_args['hcrop'], user_args['mlen'],
                             user_args['threads'], user_args['leading_trim'],
                             user_args['trail_trim'],user_args['qenc'])

    # validate trimmomatic output
    _ = mini_loop([R1, R2, uR1R2], 'fastq')

    # centrifuge classification of reads
    ccmd_pline, p_report = centrifuge(user_args['out_name'], 
                                      (R1, R2),
                                      user_args['threads'],
                                      user_args['mhlen'],
                                      user_args['cdb'], 1, pair=True,
                                      fastq=True)

    pcent_out = run(ccmd_pline, stdout=PIPE, stderr=PIPE, shell=True)
    valid, ms = validate_output(pcent_out)
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    if not valid:
        sys.exit(1)
    ccmd_upline, up_report = centrifuge(user_args['out_name'], uR1R2,
                                        user_args['threads'],
                                        user_args['mhlen'],
                                        user_args['cdb'], 1, pair=False,
                                        fastq=True)
    ucent_out = run(ccmd_upline, stdout=PIPE, stderr=PIPE, shell=True)
    valid, ms = validate_output(ucent_out)
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    if not valid:
        sys.exit(1)

    # ---- Deleting temporary files ----
    os.chdir('..')
    my_preps = (R1, R2, uR1R2, p_report, up_report)
    new_abspaths = []
    ms = 'Something went wrong.\n'
    skip = []
    for f in my_preps:
        try:
            myf = os.path.split(f)[1]
            print(f'Moving {f} to {os.getcwd()}', end='\n', file=sys.stdout,
                  flush=True)
            nflocation = os.path.join(os.getcwd(), myf)
            if os.path.isfile(nflocation):
                m = f'WARNING:\nFile {myf} already exists in {os.getcwd()} '
                m += 'and will not be overwritten.\nNewly processed file will '
                m += 'remain in the tmp directory'
                print(m, sep=' ', end='\n', file=sys.stdout, flush=True)
                skip.append(f)
            else:
                shutil.move(f, os.getcwd())
                new_abspaths.append(nflocation)
        except:
            ms += f'ERROR: attempt to move {f} to {os.getcwd()} failed\n'
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            sys.exit(1)
    if not skip:
        try:
            os.system('rm -r %s' % tmpdirname_path)
        except:
            ms += '%s does not seem to exist\n'
            ms += 'Terminating program.'
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            sys.exit(1)
        ms = "Files are ready to use with 'short_seqs' mode:\n%s" % '\n'.join(
            new_abspaths)
        print(ms, flush=True)
    return my_preps

def validate_output(res_output):

    stderr = res_output.stderr.decode('utf-8').rstrip('\n').split('\n')
    value = [True for e in stderr if e.lower().startswith('error')]
    if value:
        return False, stderr

    # check the output is not empty
    tsvfile = stderr[0].split(' ')[-1]
    cfile = os.path.abspath(tsvfile)
    if os.path.exists(cfile) and os.stat(cfile).st_size > 0:
        fcontent = head2(cfile)
        mylen = len(fcontent)
        if mylen == 1:
            cla_file = cfile[:-4]
            m = f"File {cfile} is empty. It is highly likely that all reads "
            m += "were deemed as 'unclassified' by centrifuge. Eukfinder will "
            m += f"continue unless the report {cla_file} is empty.\n"
            # make sure the classification file is no empty despite empty tsv
            if os.path.exists(cla_file) and os.stat(cla_file).st_size > 0:
                xcontent = head2(cla_file)
                if len(xcontent) != 2:
                    m += f"Report {cla_file} might be empty.\nExiting program"
                    return False, m
            return True, m
        elif mylen > 1:
            m = "No errors reported"
            return True, m
        else:
            m = "Unknown error"
            m += '\n'.join(stderr)
            return True, m
    else:
        m = f"Something is wrong with {cfile}. Terminating program."
        return False, m


def program_exists(name):
    """
    Check whether `name` is on PATH and marked as executable.
    """
    return shutil.which(name) is not None


def perform_short_seqs(user_args):
    """

    :param user_args:
    :return:
    """

    # print(user_args)
    start_time = time.time()
    stamp = time.strftime('%Y%m%d', time.gmtime(start_time))
    log = 'Short_seqs_%s.log' % stamp
    sys.stdout = open(log, 'w')
    ms = 'Eukfinder v%s is using python %s' % (__version__,
                                               platform.python_version())
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    tax_up, reads_paths, redef_reads, redef_class, \
    plast_path, map_id_path, ncpu, max_plast, \
    e_value, pid, mode, base_outname, max_mem, cdb_path,  \
    cov, mhlen, kmers = input_check_and_setup(user_args)  #aplast_path, amap_id_path,

    dirname = 'Intermediate_data'
    mkdir(dirname)
    dirpath = os.path.abspath(dirname)
    os.chdir(dirpath)
    n_cpu = cpu_chunk(ncpu, max_plast)

    ms = 'tmps_%s_%s' % (base_outname, stamp)
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    ms = 'user arguments are:\n %s' % user_args
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

    # --- dumping taxonomy if requested ---
    ncbi = ete3.NCBITaxa()
    if tax_up:
        ncbi.update_taxonomy_database()

    # --- parsing centrifuge output ---
    ms = 'Reading classification reports ...\n'
    ms += 'mode is %s' % mode
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    if mode == 'short_seqs':
        if name_check(reads_paths[0]):
            rr = '\n'.join(redef_reads)
            rp = '\n'.join(reads_paths)
            ms = 'redef_reads:\n%s\nreads paths:\n%s' % (rr, rp)
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
            nrr = [rename_reads((reads_paths[0], 'R1')),
                   rename_reads((reads_paths[1], 'R2')),
                   rename_reads((reads_paths[2], 'U'))]
            ms = '\n'.join([e[1] for e in nrr])
            new_redef_reads = [e[0] for e in nrr]
            ms += 'new redef reads are:\n%s\n' % '\n'.join(new_redef_reads)
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        else:
            new_redef_reads = reads_paths

        existing_dfs = pair_read_handler(n_cpu, max_plast,
                                         new_redef_reads, redef_class,
                                         new_redef_reads, stamp, base_outname)
    else:
        if name_check(redef_reads[0]):
            nrr = [rename_reads((redef_reads[0], 'U'))]
            m = '\n'.join([e[1] for e in nrr])
            new_redef_reads = [e[0] for e in nrr]
            ms += m
            ms += 'new redef reads are:\n%s\n' % '\n'.join(new_redef_reads)
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        else:
            new_redef_reads = reads_paths

        existing_dfs = single_read_handler(n_cpu, max_plast,
                                           new_redef_reads, redef_class,
                                           new_redef_reads, stamp,
                                           base_outname)
    # Management of plast searches
    if existing_dfs:
        print('processing existing df', existing_dfs)
        parsed_plouts = plast_search_n_results(n_cpu, plast_path,
                                               existing_dfs, e_value, pid,
                                               cov, False)

        final_dfs = create_df_taxonomy(parsed_plouts,
                                       existing_dfs, pid,
                                       map_id_path)
        # dirpath
        _ = writinggroups(final_dfs, dirpath, base_outname)

    else:
        ms = 'No files containing reads were declared. Terminating program'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)

    #   Meta-assembly and classification
    ms = 'Re-classification step'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

    cat_reads = cats(base_outname, dirpath)
    scaffolds, new_outname = assembly(cat_reads, base_outname, ncpu, max_mem,
                                      kmers)

    # Contig classification and report validation
    sreport = post_assembly(ncpu, base_outname, scaffolds,
                            cdb_path, mhlen)


    new_dfs = single_read_handler(n_cpu, max_plast,
                                  [scaffolds], [sreport],
                                  [scaffolds], stamp, new_outname)

    if new_dfs:
        parsed_plouts = plast_search_n_results(n_cpu, plast_path,
                                               new_dfs, e_value, pid, cov,
                                               True)

        fdf = create_df_taxonomy(parsed_plouts,
                                 new_dfs, pid,
                                 map_id_path)

        mycwd = os.getcwd()
        os.chdir('..')
        dir_path = os.getcwd()
        os.chdir(mycwd)
        _ = writinggroups(fdf, dir_path, new_outname)
        os.chdir('..')
        print('cwd', os.getcwd(), sep=' ', end='\n', file=sys.stdout,
              flush=True)

        # patching relocation of second centrifuge results
        centrifuge_outs = glob.glob(os.path.abspath('*centrifuge*'))
        patch_dir = os.path.abspath('tmps_%s_%s' % (new_outname, stamp))
        for out in centrifuge_outs:
            shutil.move(out, patch_dir)

        # relocate result directories
        os.chdir('..')
        #dirs_2_move = ['Classified_reads', 'Classified_contigs']
        dirs_2_move = ['Eukfinder_results']
        bdir = os.getcwd()
        for d in dirs_2_move:
            mf = os.path.join(dir_path, d)
            shutil.move(mf, bdir)

        summary_table()        
        ms = '****  RESULTS ARE READY!  ****\n'
        ms += "The 'Eukfinder_results' directory contains the sequences classified "
        ms += "as bacteria, archaea, eukaryota, unknown and Eukaryota-Unknown.\n"
        ms += "Directory Intermediate_data contains "
        ms += "temporary intermediate files to classify the sequences mentioned.\n"  
        
        ms += "\n*****   WARNING !!! *****\n"
        ms += "Each of these classified files may contain MORE than ONE "
        ms += "eukaryotic or prokaryotic taxon or a mixture of them. "
        ms += "Hence, supervised binning or \n"
        ms += "manual inspection on the file of your interest MUST be done.\n"
        ms += "We recommend to apply MyCC software (Lin & Liao, 2016) "
        ms += "(doi:10.1038/srep24175) or \nAnvio (Eren et al 2015)"
        ms += "(10.7717/peerj.1319)\n.  Please see Eukfinder user manual "
        ms += "for binning\n*****\tEND OF WARNING\t****\n"
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)


def perform_long_seqs(user_args):
    #
    start_time = time.time()
    stamp = time.strftime('%Y%m%d', time.gmtime(start_time))
    log = 'Long_seqs_%s.log' % stamp
    sys.stdout = open(log, 'w')
    ms = 'Eukfinder v%s is using python %s' % (__version__,
                                               platform.python_version())
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    tax_up, reads_path, redef_reads, plast_path, \
    map_id_path, n_cpu, max_plast, \
    e_value, pid, base_outname, cov,  cdb_path, mhlen = input_check_and_setup(
        user_args)
    n_cpu = cpu_chunk(n_cpu, max_plast)
    new_reads, m = rename_reads((redef_reads[0], 'LR'))
    print(m, sep=' ', end='\n', file=sys.stdout, flush=True)
    reads = bytearray_sets(new_reads)
    dirname = 'tmps_%s_%s' % (base_outname, stamp)
    mkdir(dirname)
    dirname_path = os.path.abspath(dirname)
    path1 = reads_path[0]
    os.chdir(dirname)

    # centrifugue classification
    if isfasta(new_reads):
        fastq=False
    else:
        fastq=True
    ccmd_upline, up_report = centrifuge(user_args['out_name'], new_reads,
                                        user_args['number_of_threads'],
                                        user_args['mhlen'],
                                        user_args['cdb'], 1, pair=False,
                                        fastq=fastq)

    cstdout = run(ccmd_upline, stdout=PIPE, stderr=PIPE, shell=True)
    # validate output
    valid, ms = validate_output(cstdout)
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    if not valid:
        sys.exit(1)
    df_un = pd.DataFrame(reads, columns=['readID'])
    dfr_un = reading_reports(up_report)

    existing_dfs = []
    # empty df to trigger plast search
    jobs = unpair_dfs(dirname_path, new_reads, df_un, path1, dfr_un,
                      n_cpu, max_plast, base_outname)
    existing_dfs.append(jobs)

    if existing_dfs:
        parsed_plouts = plast_search_n_results(n_cpu, plast_path,
                                               existing_dfs, e_value, pid,
                                               cov, True)
        fdf = create_df_taxonomy(parsed_plouts,
                                 existing_dfs, pid, map_id_path)

        mycwd = os.getcwd()
        os.chdir('..')
        dir_path = os.getcwd()
        os.chdir(mycwd)
        _ = writinggroups(fdf, dir_path, base_outname)
        os.chdir('..')

    # relocate temporary directories
    Interdata_dirname = 'Intermediate_data'
    mkdir(Interdata_dirname)
    Interdata_dir_path = os.path.abspath(Interdata_dirname)
    try:
        shutil.move(dirname_path, Interdata_dir_path)
    except:
        ms = 'Something went wrong when trying to remove:\n'
        ms += '%s\nTerminating program' % dirname_path
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)

    # Folder paths
    old_folder = 'Classified_reads'
    new_folder = 'Eukfinder_results'

    # Check if the folder exists
    if os.path.exists(old_folder) and os.path.isdir(old_folder):
        # Rename the folder
        shutil.move(old_folder, new_folder)
        print(f"Folder '{old_folder}' renamed to '{new_folder}'.")
    else:
        print(f"No folder named '{old_folder}' found.")

    # removing temporary files
    files_to_delete = ["tmp.LR.fasta", "tmp.LR.fastq"]
    for file in files_to_delete:
        if os.path.exists(file):  # Check if the file exists
            os.remove(file)       # Delete the file

    summary_table()    
    # removing temporary files
    """
    try:
        shutil.rmtree(dirname_path)
        os.remove(new_reads)
    except:
        ms = 'Something went wrong when trying to remove:\n'
        ms += '%s\nTerminating program' % dirname_path
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)
    """
    ms = '****  RESULTS ARE READY!  ****\n'
    ms += "The 'Eukfinder_results' directory contains the sequences classified "
    ms += "as bacteria, archaea, eukaryota, unknown and Eukaryota-Unknown.\n"
    ms += "Directory Intermediate_data contains "
    ms += "temporary intermediate files to classify the sequences mentioned.\n" 
    ms += "\n*****   WARNING !!! *****\n"
    ms += "Each of these classified files may contain MORE than ONE "
    ms += "eukaryotic or prokaryotic taxon or a mixture of them. "
    ms += "Hence, supervised binning or \n"
    ms += "manual inspection on the file of your interest MUST be done.\n"
    ms += "We recommend to apply MyCC software (Lin & Liao, 2016) "
    ms += "(doi:10.1038/srep24175) or \nAnvio (Eren et al 2015)"
    ms += "(10.7717/peerj.1319)\n.  Please see Eukfinder user manual "
    ms += "for binning\n*****\tEND OF WARNING\t****\n"
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)
    return 'Done'

# NOTE: It is perhaps worthwhile to verify checksum in the future
def perform_download_db(user_args):
    path = user_args['path']
    name = user_args['name']

    if not os.path.isdir(f"{path}/{name}"):
        os.mkdir(f"{path}/{name}")
        print(f"Created {path}/{name}\n")

    try:
        while True:
            user_input = input(f"Would you like to download all databases ({_all_db[1]})? (yes/no)\n>")

            if user_input == "yes":
                print("\nDownloading...")
                print("Please do not close tab, this may take a while...")
                urllib.request.urlretrieve(_all_db[2], f"{path}/{name}/{_all_db[3]}")

                print("\nDecompressing...")
                file = tarfile.open(f"{path}/{name}/{_all_db[3]}")
                file.extractall(f"{path}/{name}")
                file.close()

                print("\nDecompressing individual databases...")
                for content in _database.values():
                    file = tarfile.open(f"{path}/{name}/{content[3]}")
                    file.extractall(f"{path}/{name}")
                    file.close()
                    os.remove(f"{path}/{name}/{content[3]}")

                print("\nUpdating default database paths in json file...")
                new_json_data = {
                    "centrifuge_db": f"{path}/{name}/{_cdb}",
                    "plast_db": f"{path}/{name}/{_pdb}",
                    "plast_map": f"{path}/{name}/{_pmap}"
                }
                update_json(new_json_data)

                os.remove(f"{path}/{name}/{_all_db[3]}")
                sys.exit(f"\nDatabases downloaded and decompressed in {path}/{name}, exiting...")
            elif user_input == "no":
                while True:
                    print("\nPlease select database(s) which you would like to install, separated by spaces (e.g., 1 2).\n")
                    # TODO: this shouldn't be hardcoded
                    print(f"1. {_database['1'][0]} - {_database['1'][1]}")
                    print(f"2. {_database['2'][0]} - {_database['2'][1]}")
                    print(f"3. {_database['3'][0]} - {_database['3'][1]}")
                    print(f"4. {_database['4'][0]} - {_database['4'][1]}")
                    print(f"5. {_database['5'][0]} - {_database['5'][1]}")
                    print(f"6. {_database['6'][0]} - {_database['6'][1]}")
                    user_input = input("\nOr type exit, if you would like to skip for now:\n>")

                    if user_input == "exit":
                        os.rmdir(f"{path}/{name}")
                        print(f"\nDeleted {path}/{name}\n")
                        sys.exit("No downloads, exiting...")

                    selected = user_input.split(" ") if user_input.strip() else []

                    if not selected:
                        print("\nInvalid option. Enter indices separated by spaces.\n")
                        continue

                    print("\nDownloading...\n")

                    for index in selected:
                        if index in _database.keys():
                            db_path = f"{path}/{name}/{index}_{_database[index][0]}"
                            os.mkdir(db_path)

                            print(f"Downloading {_database[index][0]}...")
                            urllib.request.urlretrieve(_database[index][2],f"{db_path}/{_database[index][3]}")

                            print(f"Decompressing {_database[index][0]}...")
                            file = tarfile.open(f"{db_path}/{_database[index][3]}")
                            file.extractall(db_path)
                            file.close()

                            os.remove(f"{db_path}/{_database[index][3]}")

                            print(f"{_database[index][0]} downloaded and decompressed.\n")
                        else:
                            print(f"WARNING: Unrecognized index {index}, skipping...\n")

                    break

                sys.exit(f"Database(s) downloaded in {path}/{name}, exiting...")

            print("\nInvalid option. Enter yes or no.\n")
    except Exception as e:
        exit(f"Error: {e}")

def short_seqs(args):

    bname = args.o
    reads = args.r1, args.r2, args.ur1, args.ur2
    classification = args.pclass, args.uclass
    paths = args.m, args.a, args.p, args.cdb, args.m2, args.p2
    params = args.n, args.t, args.e, args.z, args.mhlen
    return reads, classification, paths, params, bname


def read_prep(args):

    bname = args.o
    reads = args.r1, args.r2
    host = args.hg
    threads = args.n
    adapters = args.i
    params = args.hcrop, args.l, args.t, args.wsize, \
        args.qscore, args.mlen, args.mhlen, args.qenc
    return bname, reads, threads, adapters, params, host

def read_prep_env(args):

    bname = args.o
    reads = args.r1, args.r2
    threads = args.n
    adapters = args.i
    params = args.hcrop, args.l, args.t, args.wsize, \
        args.qscore, args.mlen, args.mhlen, args.qenc
    return bname, reads, threads, adapters, params

def long_seqs(args):

    bname = args.o
    reads = args.u
    paths = args.m, args.a, args.p
    params = args.n, args.t, args.e, args.z, args.mhlen
    return reads, paths, params, bname

def download_db(args):
    path = args.path
    return path

def summary_table():
    # Change directory to 'Eukfinder_results'
    try:
        os.chdir("Eukfinder_results")
        #print("Changed directory to 'Eukfinder_results'")
    except FileNotFoundError:
        print("Error: 'Eukfinder_results' directory not found. Exiting...")
        sys.exit(1)

    # Generate tables from FASTA/FASTQ files
    for f in glob.glob("*.f*"):
        basename = os.path.splitext(f)[0]  # Improved filename handling
        output = f"{basename}.table"

        cmd = f"seqkit fx2tab --length --name --header-line {f} -o {output}"
        _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)

    # Dictionary to store results
    summary = {
        "Group": [],
        "#Seq": [],
        "Total size(bp)": []
    }

    # Process table files and collect data
    for file in glob.glob("*.table"):
        if file.endswith(".un.table"):
            group = file.split(".")[-3]
        else:
            group = file.split(".")[-2]

        # Read data, skip header
        try:
            data = pd.read_csv(file, sep='\t', header=0)
        except pd.errors.EmptyDataError:
            print(f"Warning: {file} is empty or malformed. Skipping...")
            continue

        # Append results
        summary["Group"].append(group)
        summary["#Seq"].append(len(data))
        summary["Total size(bp)"].append(data['length'].sum())

        # Delete temporary `.table` files
        os.remove(file)
        #print(f"Deleted temporary file: {file}")

    # Create and save summary table
    summary_df = pd.DataFrame(summary)
    output_file = "summary_table.txt"
    summary_df.to_csv(output_file, sep='\t', index=False)

    print(f"Summary table has been created: Eukfinder_results/{output_file}")

def update_json(new_json_data):

    with open(_json_path, "w") as json_file:
        json.dump(new_json_data, json_file, indent=4)

def read_json():

    with open(_json_path, "r") as json_file:
        json_data = json.load(json_file)

    return json_data

def parse_arguments(json_data):

    myargs = {
        '-n': ['--number-of-threads', str, '20', 'Number of threads', False],
        '-z': ['--number-of-chunks', str, '2', 'Number of chunks to split a '
                                          'file', False],
        '-t': ['--taxonomy-update', str, 'False', 'Set to True the first '
               'time the program is used. Otherwise set to False', False],
        '-p': ['--plast-database', str, json_data["plast_db"], 'path to plast database', False],
        '-m': ['--plast-id-map', str, json_data["plast_map"], 'path to taxonomy map for '
                                      'plast database', False],
        '--cdb': ['--centrifuge-database', str, json_data["centrifuge_db"], 'path to centrifuge '
                                                'database', False],
        '-e': ['--e-value', float, 0.01, 'threshold for plast searches', False],
        '--pid': ['--percent_id', float, 60, 'percentage identity for '
                                         'plast searches', False],
        '--cov': ['--coverage', float, 10, 'percentage coverage for '
                                       'plast searches', False],
        '--max_m': ['--max_memory', str, "300", 'Maximum memory allocated to '
                                         'carry out an assembly', False],
        '-k': ['--kmers', str, "21, 33, 55", 'kmers to use during assembly. '
               'These must be odd and less than 128. default is 21,33,55',
               False],
        '--mhlen': ['--min-hit-length', int, 25, 'Maximum memory allocated to '
                                             'carry out an assembly', False],
        '--pclass': ['--p-reads-class', str, None, 'Classification for '
                                             'pair end reads', False],
        '--uclass': ['--u-reads-class', str, None, 'Classification for '
                                             'un-pair end reads', False]
    }

    parser = argparse.ArgumentParser(prog='eukfinder')
    subparsers = parser.add_subparsers()

    #  ---  second level parser for pair mode ---  #
    parser_short_seqs = subparsers.add_parser('short_seqs')
    group1 = parser_short_seqs.add_argument_group('Required arguments',
                                                  'Description')
    group1.add_argument('--r1', '--reads-r1', type=str,
                        help='left reads', required=True)
    group1.add_argument('--r2', '--reads-r2', type=str,
                        help='right reads', required=True)
    group1.add_argument('--un', '--un-pair-reads', type=str,
                        help='orphan reads', required=True)
    group1.add_argument('-o', '--out_name', type=str,
                        help='output file basename', required=True)
    for key in myargs:
        try:
            group1.add_argument(key, myargs[key][0],
                                type=myargs[key][1],
                                default=myargs[key][2],
                                help=myargs[key][3],
                                required=myargs[key][4])
        except:
            parser_short_seqs.add_argument(key, myargs[key][0],
                                           type=myargs[key][1],
                                           default=myargs[key][2],
                                           help=myargs[key][3])

    #  ---  second level parser for unpair mode ---  #
    #  ---  second level parser for read_prep ---  #
    parser_read_prep = subparsers.add_parser('read_prep')
    group2 = parser_read_prep.add_argument_group('Required arguments',
                                                 'Description')
    group2.add_argument('--r1', '--reads-r1', type=str,
                        help='left reads', required=True)
    group2.add_argument('--r2', '--reads-r2', type=str,
                        help='right reads', required=True)
    group2.add_argument('-n', '--threads', type=int,
                        help='number of threads', required=True)
    group2.add_argument('-i', '--illumina-clip', type=str,
                        help='adaptor file', required=True)
    group2.add_argument('--hcrop', '--head-crop', type=int,
                        help='head trim', required=True)
    group2.add_argument('-l', '--leading-trim', type=int,
                        help='leading trim', required=True)
    group2.add_argument('-t', '--trail-trim', type=int,
                        help='trail trim', required=True)
    group2.add_argument('--wsize', '--window-size', type=int,
                        help='sliding window size', required=True)
    group2.add_argument('--qscore', '--quality-score', type=int,
                        help='quality score for trimming', required=True)
    group2.add_argument('--mlen', '--min-length', type=int,
                        help='minimum length', required=True)
    group2.add_argument('--mhlen', '--min-hit-length', type=int,
                        help='minimum hit length', required=True)
    group2.add_argument('--hg', '--host-genome', type=str,
                        help='host genome to get map out', required=True)
    group2.add_argument('-o', '--out_name', type=str,
                        help='output file basename', required=True)
    group2.add_argument('--cdb', '--centrifuge-database', type=str,
                        default= json_data["centrifuge_db"],
                        help='path to centrifuge database', required=False)
    group2.add_argument('--qenc', '--quality-encoding', type=str,
                        help='quality encoding for trimmomatic', default='auto', required=False)


    # --- second level parser for long read mode ---
    parser_long_seqs = subparsers.add_parser('long_seqs')
    group3 = parser_long_seqs.add_argument_group('Required arguments',
                                                 'Description')
    group3.add_argument('-l', '--long-seqs', type=str,
                        help='long sequences file', required=True)
    group3.add_argument('-o', '--out_name', type=str,
                        help='output file basename', required=True)
    group3.add_argument('--mhlen', '--min-hit-length', type=int,
                        help='minimum hit length', required=True)
    group3.add_argument('--cdb', '--centrifuge-database', type=str,
                        default= json_data["centrifuge_db"],
                        help='path to centrifuge database', required=False)

    myargs_lr = {
        '-n': ['--number-of-threads', str, '20', 'Number of threads', False],
        '-z': ['--number-of-chunks', str, '2', 'Number of chunks to split a'
                                          ' file', False],
        '-t': ['--taxonomy-update', str, 'False', 'Set to True the first '
                                         'time the program is used. '
                                         'Otherwise set to False', False],
        '-p': ['--plast-database', str, json_data["plast_db"], 'path to plast database', False],
        '-m': ['--plast-id-map', str, json_data["plast_map"], 'path to taxonomy map for '
                                      'plast database', False],
        '-e': ['--e-value', float, 0.01, 'threshold for plast searches', False],
        '--pid': ['--percent_id', float, 60, 'percentage identity for '
                                         'plast searches', False],
        '--cov': ['--coverage', float, 10, 'percentage coverage for '
                                       'plast searches', False],
    }

    #  ---  second level parser for read_prep_env ---  #
    parser_read_prep_env = subparsers.add_parser('read_prep_env')
    group4 = parser_read_prep_env.add_argument_group('Required arguments',
                                                 'Description')
    group4.add_argument('--r1', '--reads-r1', type=str,
                        help='left reads', required=True)
    group4.add_argument('--r2', '--reads-r2', type=str,
                        help='right reads', required=True)
    group4.add_argument('-n', '--threads', type=int,
                        help='number of threads', required=True)
    group4.add_argument('-i', '--illumina-clip', type=str,
                        help='adaptor file', required=True)
    group4.add_argument('--hcrop', '--head-crop', type=int,
                        help='head trim', required=True)
    group4.add_argument('-l', '--leading-trim', type=int,
                        help='leading trim', required=True)
    group4.add_argument('-t', '--trail-trim', type=int,
                        help='trail trim', required=True)
    group4.add_argument('--wsize', '--window-size', type=int,
                        help='sliding window size', required=True)
    group4.add_argument('--qscore', '--quality-score', type=int,
                        help='quality score for trimming', required=True)
    group4.add_argument('--mlen', '--min-length', type=int,
                        help='minimum length', required=True)
    group4.add_argument('--mhlen', '--min-hit-length', type=int,
                        help='minimum hit length', required=True)
    group4.add_argument('-o', '--out_name', type=str,
                        help='output file basename', required=True)
    group4.add_argument('--cdb', '--centrifuge-database', type=str,
                        help='path to centrifuge database', required=True)
    group4.add_argument('--qenc', '--quality-encoding', type=str,
                        help='quality encoding for trimmomatic', default='auto', required=False)

   
    parser_download_db = subparsers.add_parser("download_db")
    parser_download_db.add_argument("-n", "--name", type=str, default="eukfinder_databases",
                                    help="directory name for storing the databases")
    parser_download_db.add_argument("-p", "--path", type=str, default=f"{os.path.expanduser('~')}/.eukfinder",
                                    help="filesystem path for storing the databases")

    for key in myargs_lr:
        try:
            group3.add_argument(key, myargs_lr[key][0],
                                type=myargs_lr[key][1],
                                default=myargs[key][2],
                                help=myargs_lr[key][3],
                                required=myargs_lr[key][4])
        except:
            parser_long_seqs.add_argument(key, myargs_lr[key][0],
                                          type=myargs_lr[key][1],
                                          default=myargs[key][2],
                                          help=myargs_lr[key][3])

    parser_short_seqs.set_defaults(func=short_seqs)
    parser_long_seqs.set_defaults(func=long_seqs)
    parser_read_prep.set_defaults(func=read_prep)
    parser_read_prep_env.set_defaults(func=read_prep_env)
    parser_download_db.set_defaults(func=download_db)

    return parser.parse_args()

def main():
    # json creation
    if not os.path.exists(_json_path):
        os.makedirs(os.path.dirname(_json_path), exist_ok=True)
        update_json(
            {
                "centrifuge_db": "",
                "plast_db": "",
                "plast_map": ""
            }
        )

    json_data = read_json()
    args = parse_arguments(json_data)

    if len(sys.argv) == 1:
        print('Try Eukfinder.py -h for more information', sep=' ',
              end='\n', file=sys.stdout, flush=True)
        sys.exit(1)

    dic_args = vars(args)
    if dic_args['func'].__name__ == 'read_prep':
        perform_prep(dic_args)
    elif dic_args['func'].__name__ == 'short_seqs':
        perform_short_seqs(dic_args)
    elif dic_args['func'].__name__ == 'long_seqs':
        perform_long_seqs(dic_args)
    elif dic_args['func'].__name__ == 'read_prep_env':
        perform_prep_env(dic_args)
    elif dic_args["func"].__name__ == "download_db":
        perform_download_db(dic_args)


if __name__ == '__main__':
    main()
