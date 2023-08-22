#!/usr/bin/env python
import os
import re
import sys
import ete3
import glob
import time
import shutil
import platform
import numpy as np
import pandas as pd
from subprocess import PIPE, run
from joblib import Parallel, delayed

# ls -lthr /home/dsalas/Shared/Eukfinder
# -rwxrwxrwx 1 dsalas roger  76K Feb 14  2021 eukfinder_Euk-Unk_v1.2.3.py

'''
[diff_v1.2.1_v1.2.3.txt]
14a15,22
> ' ' '
> 919,920c928,929
> <                   'acc2tax_nucl_all.txt', 'gi_taxid_nucl.dmp',
> <                   'gi_taxid_prot.dmp']
> ---
> >                    'acc2tax_nucl_all.txt']#, 'gi_taxid_nucl.dmp',
> >                    #'gi_taxid_prot.dmp']
> ' ' '
799,800c807
<                      ('Misc', label, pref, Misc, reads),
<                      ('Euk', label, pref, Euk, reads)]
---
>                      ('Misc', label, pref, Misc, reads)]
920,921c927,928
<                   'acc2tax_nucl_all.txt']#, 'gi_taxid_nucl.dmp',
<                   #'gi_taxid_prot.dmp']
---
>                    'acc2tax_nucl_all.txt']#, 'gi_taxid_nucl.dmp',
>                    #'gi_taxid_prot.dmp']
1533c1540
<     cen = glob.glob('%s*' % user_args['cdb'])
---
>     cen = glob.glob('%s*.cf' % user_args['cdb'])

'''
# --- preparation ---

def trimming(bn, reads1, reads2, adapath, wsize, qscore, headcrop,
             mlenght, threads, leading_trim, trail_trim):
    r1_out = '%sR1PT.fq %sR1unPT.fq ' % (bn, bn)
    r2_out = '%sR2PT.fq %sR2unPT.fq ' % (bn, bn)
    cmd = 'trimmomatic PE -threads %s -trimlog %s.trim.log ' % (threads, bn)
    cmd += '%s %s %s %s ILLUMINACLIP:%s:2:30:10 HEADCROP:%s LEADING:%s ' % (
        reads1, reads2, r1_out, r2_out, adapath, headcrop, leading_trim)
    cmd += 'TRAILING:%s SLIDINGWINDOW:%s:%s MINLEN:%s' % (trail_trim,
                                                          wsize, qscore,
                                                          mlenght)

    ms = 'trimmomatic cmd_line:\n%s\n' % cmd
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)

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
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    if len(glob.glob('*bt2')) != 0:
        return bn
    else:
        ms = 'bowtie indexes could not be built.\n'
        ms += 'Exiting program\n'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)


def bowtie2(bn, threads, r1, r2, catr1r2, hindex, typ):
    if typ == 'fastq-host':
        cmd = 'bowtie2 --local --phred33 -q --threads %s -x %s' % (threads,
                                                                   hindex)
        cmd += ' -1 %s -2 %s -U %s -S %s.sam ' % (r1, r2, catr1r2, bn)
        cmd += '--un-conc %s_p.fastq --un %s_un.fastq --no-unal' % (bn, bn)
        ms = 'bowtie2 fastq-host cmd_line:\n%s' % cmd
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
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
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        bt2fqR1 = os.path.abspath('%s_bowtie2.1.fastq' % bn)
        bt2fqR2 = os.path.abspath('%s_bowtie2.2.fastq' % bn)
        bt2fqun = os.path.abspath('%s_bowtie2.un.fastq' % bn)
        return bt2fqR1, bt2fqR2, bt2fqun, cmd

    else:
        cmd = 'bowtie2 --local --threads %s -x %s ' % (threads, hindex)
        cmd += '-U %s -S %s.sam --al %s_bowtie2.fasta' % (catr1r2, bn, bn)
        cmd += '--no-unal'
        ms = 'bowtie2 fasta cmd_line:\n%s' % cmd
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
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
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        return gline, report

    else:
        bn_r1r2 = bn_tuple
        report = os.path.join(os.getcwd(),'%s_centrifuge_UP' % bn)

        gline += '-U %s -S %s --report-file %s.tsv ' % (bn_r1r2,
                                                        report, report)
        ms = 'centrifuge unpair cmd_line:\n%s' % gline
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        return gline, report


def cats(b_outname):
    os.chdir('..')
    cwd = os.getcwd()
    # up_dir = cwd.split('/Temp')[0]
    outr1 = os.path.join(cwd, '%s.EUnkR1.fq' % b_outname)
    outr2 = os.path.join(cwd, '%s.EUnkR2.fq' % b_outname)
    outr1r2 = os.path.join(cwd, '%s.EUnk.fq' % b_outname)
    return outr1, outr2, outr1r2


# ---  classification ---

def reading_reports(infile):
    '''

    :param infile:
    :return:
    '''
    ms = 'Processing %s ...' % infile
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
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
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        return table
    else:
        m = 'report file has not been declared '
        m += 'or the declared file is not in the directory. The file or '
        m += 'a symbolic link must exist in the working directory.'
        print(m, sep=' ', end='\n', file=sys.stdout, flush=False)
        return pd.DataFrame()


def minimal(tupla):
    '''

    :param tupla:
    :return:
    '''
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
        print(ms, sep=' ', end='\n', flush=False)
        print('empty chunk is :\n', newchunk)
        return newchunk


def parseplastoutput(single_plout, ident, cov, lr):
    '''

    :param single_plout:
    :return:
    '''
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
    '''

    :param CPUs:
    :param chunk:
    :param pDBpath:
    :param evalid:
    :return:
    '''
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
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    return line


def my_run(cmdline):
    '''

    :param cmdline:
    :return:
    '''
    run(cmdline, stderr=PIPE,
        stdout=PIPE, shell=True, bufsize=1,
        close_fds=True, universal_newlines=True)


def plastsearches(CPUs, pDBpath, evalue, pattern):
    '''

    :param CPUs:
    :param pDBpath:
    :param evalue:
    :param pattern:
    :return:
    '''
    abs_pattern = os.path.join(os.getcwd(), pattern)
    queries = glob.glob(abs_pattern)
    # queries = glob.glob(pattern): i added 2 lines above instead of this one
    label = pattern[0].upper()
    if queries:
        ms = 'plast search in progress ...'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
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
    '''

    :param list_of_results:
    :return:
    '''
    time.sleep(5)
    label = list_of_results[0][1]
    ms = 'parsing plast outputs for %s ....' % label
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    allplast_outs = []
    pouts = [parseplastoutput(out[0], pid, cov, lr) for out in
             list_of_results if os.stat(out[0]).st_size != 0]
    if pouts:
        plast_outputs = pd.concat(pouts)
        allplast_outs.append((label, plast_outputs))
    else:
        allplast_outs.append((label, pd.DataFrame()))
    return allplast_outs


def matchmaker(ncbi, regex, readid, taxid):
    '''

    :param ncbi:
    :param regex:
    :param readid:
    :param taxid:
    :return:
    '''
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
        return readid, np.NaN


def binningbytaxonomy(report_df):
    '''
    slice reports by taxonomic domain
    '''
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
    '''

    :param binned_lists:
    :param main_df:
    :return:
    '''
    maintaxidlist = [tupla for sublist in
                     binned_lists for tupla in sublist]

    # creating a df for centrifuge classified reads
    if maintaxidlist:
        centrifuged_taxId = pd.DataFrame(maintaxidlist,
                                         columns=['readID', 'Group'])
        main_df = main_df.merge(centrifuged_taxId, on='readID', how='left')
    else:
        main_df = main_df.assign(Group=np.NaN)
    return main_df


def dataframe_collection(df, npartitions):
    '''

    :param df:
    :param npartitions:
    :return:
    '''
    df_size = len(df.index)
    if npartitions > df_size:
        ms = 'Number of partitions is larger '
        ms += 'than dataframe size. The number of'
        ms += 'partitions will be set to: %s' % df_size
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        npartitions = df_size
    df_split = np.array_split(df, npartitions)
    return df_split


def directplast(new_dir_path, file_name, suffix, name, max_plast, cpus):
    '''

    :param new_dir_path:
    :param file_name:
    :param suffix:
    :param name:
    :param max_plast:
    :param cpus:
    :return:
    '''
    try:
        cmd = 'seqkit fq2fa %s | seqkit split -p %s -j %s -' % (file_name,
                                                                max_plast,
                                                                cpus)
    except:
        cmd = 'seqkit split -1 %s -p %s -j %s ' % (file_name, max_plast, cpus)

    ms = 'directplast cmd is:\n%s' % cmd
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)

    # change name of output and moving it to a different directory
    new_name = '%s%s' % (name, suffix)
    myfiles = relocating_files(new_name, new_dir_path)
    return myfiles


def relocating_files(suffix, new_dir_path):
    '''

    :param suffix:
    :param new_dir_path:
    :return:
    '''
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
    ms = 'cmd isfasta is %s' % cmd
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    line1 = _.stdout.decode('utf-8')
    if line1.startswith('>'):
        return True
    else:
        return False


def split_reads(new_dir_path, reads_to_select, read_file,
                suffix, name, cpus, max_plast):
    '''

    :param new_dir_path:
    :param reads_to_select:
    :param read_file:
    :param suffix:
    :param name:
    :param cpus:
    :param max_plast:
    :return:
    '''
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
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)

    # change name of output and moving it to a different directory
    # new_name = '%s%s' % (name, suffix)
    myfiles = relocating_files(suffix, new_dir_path)
    return myfiles


def slicing(new_dir_path, main_df, file_name, df_class_report,
            suffix, name, cpus, max_plast):
    '''
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
    '''
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
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    #
    full_df = taxpickup(dfs, main_df)
    _ = full_df.to_csv('%sfull_df.tsv' % name, sep='\t', header=True)
    slice = full_df[full_df['Group'].isnull()]
    slice = set(slice['readID'].tolist())

    # double - check that not-nulls go to plast
    if slice:
        file_names = split_reads(new_dir_path, slice, file_name, suffix,
                                 name, cpus, max_plast)
        print(file_names, sep=' ', end='\n', file=sys.stdout, flush=False)
        return full_df, file_names
    else:
        m = 'Nothing to slice!'
        print(m, sep=' ', end='\n', file=sys.stdout, flush=False)
        return pd.DataFrame(), m


def Taxonomy_translation(acc2DBpath, pre_acc2tax_df, suffix, typ):
    '''
     Executable of acc2tax MUST be in the environmental path
    :param acc2DBpath:
    :param pre_acc2tax_df:
    :param suffix:
    :param typ:
    :return:
    '''
    # writing acc2taxIN as input for acc2tax software
    base_fname = string_out(pre_acc2tax_df, suffix, typ)
    ms = 'Extracting taxonomy with acc2tax...'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
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
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        return taxdump


def customtaxtranslation(seqidtaxid_map, parseby, plast_output):
    '''

    :param seqidtaxid_map:
    :param parseby:
    :param plast_output:
    :return:
    '''
    # subject ID and first column of seqidtaxid_map must match
    ms = 'Using custom taxonomy ...'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    colnames = ['subject ID', 'taxID', 'Group']
    taxid_map = pd.read_table(seqidtaxid_map, sep='\t',
                              header=None, names=colnames,
                              dtype={0: str, 1: str, 2: str})
    taxids = plast_output.merge(taxid_map, on='subject ID',
                                how='inner')
    # for now make sure that is 'Group'
    taxids = taxids.loc[:, ['readID', parseby]]
    taxids = taxids.groupby('readID', as_index=False).first()
    return taxids


def string_out(df, suffix, typ='RefSeq'):
    '''

    :param df:
    :param suffix:
    :param typ:
    :return:
    '''
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
    '''

    :param gdf:
    :param group:
    :param label:
    :return:
    '''
    try:
        x = gdf.get_group(group)
        _ = x.to_csv('%s_grouped_by_%s.tsv' % (label, group),
                     sep='\t', header=True)
        val = set(gdf.get_group(group).loc[
                  :, 'readID'].tolist())
    except:
        ms = 'empty bin %s %s' % (label, group)
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        val = set()
    return val


def tokeep(grouped_df, group,
           excludegroup, reject_set, label):
    '''
    Subsetting by exclusion of 'rejecting set'
    :param grouped_df:
    :param group:
    :param excludegroup:
    :param reject_set:
    :param label:
    :return:
    '''
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
    '''
    Checking if intersections among bins exist
    :param set1:
    :param set2:
    :param set3:
    :return:
    '''
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
    '''
    Grouping by major categories and check to deduplicate dfs
    :param label:
    :param pref:
    :param df:
    :param reads:
    :return:
    '''
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
        # enforcing eukaryota and unkown to merge
        intersect, Bact, Arch, Misc = intersects(Bact, Arch, Vir)
        Euk = Euk.difference(intersect)
        Preset = set().union(*[Bact, Arch, Vir, Euk])
        Unk = pre_unknown.union(other_dfs(merged))
        Unk = Unk.difference(Preset)
        EUnk = set(Euk).union(Unk)
        sel_reads = [('Bact', label, pref, Bact, reads),
                     ('Arch', label, pref, Arch, reads),
                     ('EUnk', label, pref, EUnk, reads),
                     ('Misc', label, pref, Misc, reads),
                     ('Euk', label, pref, Euk, reads)]
        return sel_reads
    else:
        Unk = pre_unknown
        sel_reads = [('Unk', label, pref, Unk, reads)]
        return sel_reads


def nullandmerged(df1, df2):
    '''

    :param df1:
    :param df2:
    :return:
    '''
    if not 'Group' in df1.columns:
        df1 = df1.assign(Group=np.NaN)
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
    '''

    :param df:
    :return:
    '''
    if not 'Group' in df.columns:
        df = df.assign(Group=np.NaN)
    pre_unknown = set(df['readID'].tolist())
    merged = pd.DataFrame()
    ms = 'In emtpy: pre_unknown\n %s, merged\n%s' % (pre_unknown, merged)
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    return merged, pre_unknown


def input_check_and_setup(user_args):
    '''

    :param user_args: dictionary containing arguments
    :return: type: list of the checked arguments
    '''

    # --- Global required values ---
    map_id_path = user_args['plast_id_map']
    plast_path = user_args['plast_database']
    acc2tax_path = user_args['acc2tax_database']
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
        declared = reads_paths + redef_class
    else:
        declared = reads_paths

    if any(i is 'None' for i in declared):
        ms = '\nIt seems that at least one declared file DOES NOT exist '
        ms += 'in the directory and has been labeled as "None". \n'
        ms += 'Please check your command line.\nDeclared files are:\n'
        ms += '*** Reads and Classification files are : ***\n'
        ms += '%s\n' % '\n'.join(declared)
        ms += '---  Exiting program   ---'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        sys.exit(0)

    if not os.path.exists(plast_path):
        ms = '\n\nNo plast database found. Exiting program'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        sys.exit(0)

    if not os.path.exists(map_id_path):
        format_check = plast_path
        ms = '\n\nNo map for plast_database was found. Exiting program'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        sys.exit(0)

    acc2tax_fn = ['nodes.dmp', 'names.dmp', 'acc2tax_prot_all.txt',
                  'acc2tax_nucl_all.txt']#, 'gi_taxid_nucl.dmp',
                  #'gi_taxid_prot.dmp']
    for fn in acc2tax_fn:
        apath = os.path.join(acc2tax_path, fn)
        if not os.path.exists(apath):
            ms = '\n\nSomething is wrong with the acc2tax database.\n'
            ms += 'Please make sure that you are specifiying the path to a '
            ms += 'directory that contains the following files: '
            ms += 'nodes.dmp, names.dmp, prot_all.txt, nucl_all.txt, '
            ms += 'gi_taxid_nucl.dmp, gi_taxid_prot.dmp\nExiting program..'
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
            sys.exit(0)

    if e_value is None:
        e_value = 1.0
    if pid is None:
        pid = 10.0
    if cov is None:
        cov = 0.0

    if len(glob.glob('%s.*.cf' % user_args['cdb'])) != 4:
        cdb = '\nCentrifuge database does not exist in path.\n'
        cdb += 'Exiting program.\n'
        print(cdb, sep=' ', end='\n', file=sys.stdout, flush=False)
        sys.exit(0)
    cdb_path = user_args['cdb']
    minhit_len = user_args['mhlen']

    if not user_args['func'].__name__ == 'long_seqs':
        aplast_path = user_args['ancillary_plast_database']
        amap_id_path = user_args['ancillary_plast_id_map']
        # max_len = user_args['mlen']
        max_mem = user_args['max_m']
        # cdb_path = user_args['cdb']
        if aplast_path is None or not os.path.exists(aplast_path):
            aplast_path = plast_path
            amap_id_path = map_id_path
            ms = '\nNo ancillary plast database found. Current plast_database '
            ms += 'will be used for re-classification\n'
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)

        if os.path.exists(aplast_path):
            if not os.path.exists(amap_id_path):
                ms = '\n\nNo map for plast_database was found. Exiting program'
                print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
                sys.exit(0)

        # if len(glob.glob('%s.*.cf' % user_args['cdb'])) != 4:
        #     cdb = '\nCentrifuge database does not exist in path.\n'
        #     cdb += 'Exiting program.\n'
        #     print(cdb, sep=' ', end='\n', file=sys.stdout, flush=False)
        #     sys.exit(0)

        arguments = [tax_up, reads_paths, redef_reads, redef_class,
                     plast_path, map_id_path, acc2tax_path, n_cpu,
                     max_plast, e_value, pid, user_args['func'].__name__,
                     base_outname, max_mem, cdb_path,
                     aplast_path, amap_id_path, cov, minhit_len]
    else:
        arguments = [tax_up, reads_paths, redef_reads, plast_path,
                     map_id_path, acc2tax_path, n_cpu, max_plast,
                     e_value, pid, base_outname, cov,  cdb_path, minhit_len]

    return arguments


def bytearray_sets(infile):
    '''

    :param infile:
    :return:
    '''
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
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    return headers


def mkdir(dir_name):
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    else:
        pass


def cpu_chunk(cpus, maxplast):
    '''

    :param cpus:
    :param maxplast:
    :return:
    '''
    cpus = int(cpus)
    max_plast = int(maxplast)
    if cpus < max_plast:
        line = 'Not enough CPUs have been allocated... '
        line += 'exiting program with status 1'
        print(line, sep=' ', end='\n', file=sys.stdout, flush=False)
        sys.exit(1)
    else:
        max_cpu = int(cpus / max_plast)
        return max_cpu


def getting_prefix(file_path):
    '''

    :param file_paths:
    :return:
    '''
    fname = file_path.split(os.sep)[-1]
    fname = fname.lower()
    prefix = [None if fname is None else
              re.split(r'.fastq$|.fq$', fname)[0]][0]
    return prefix


def rename_reads(readfile):
    '''

    :param readfile:
    :return:
    '''
    read, orient = readfile
    original_path = os.path.abspath(read)
    ftype = re.split(r'.fastq$|.fq$', read.lower())
    if len(ftype) == 1:
        tail = 'fasta'
    else:
        tail = 'fastq'
    newrfile_path = os.path.join(os.getcwd(), 'tmp.%s.%s' % (orient, tail))
    cmd = 'seqkit replace -p "\s.+" %s -o %s' % (original_path, newrfile_path)
    ms = 'rename reads cmd %s' % cmd
    print(ms, sep='\t', end='\n', flush=False)
    _ = run(cmd, stderr=PIPE, stdout=PIPE, shell=True)
    return newrfile_path


def pair_read_handler(cpus, max_files, files_paths, reports_list,
                      readfile_list, stamp, base_outname):
    '''

    :param cpus:
    :param max_files:
    :param files_paths:
    :param reports_list:
    :param readfile_list:
    :param stamp:
    :return:
    '''
    ms = 'Pair read handler..'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    reports = [reading_reports(inf) for inf in reports_list]
    reads = [bytearray_sets(readfile) for readfile in readfile_list]
    dirname = 'tmps_%s_%s' % (base_outname, stamp)
    mkdir(dirname)
    dirname_path = os.path.abspath(dirname)
    ms = 'working paths are:\n', files_paths
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
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
    '''

    :param cpus:
    :param max_files:
    :param files_paths:
    :param reports_list:
    :param readfile_list:
    :param stamp:
    :return:
    '''
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
    '''

    :param dirname_path:
    :param files_paths:
    :param df_un:
    :param path1:
    :param dfr_un:
    :param cpus:
    :param max_files:
    :param name:
    :return:
    '''
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
    '''

    :param CPUs:
    :param pDBpath:
    :param working_dfs:
    :param evalue:
    :return:
    '''

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
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    return plast_outputs


def create_df_taxonomy(plast_outputs_list, taxDBpath, working_dfs,
                       pid, seqidtaxid_map=None):
    '''

    :param plast_outputs_list:
    :param taxDBpath:
    :param working_dfs:
    :param percentid:
    :param seqidtaxid_map:
    :return:
    '''
    cladfs = []
    for entry in plast_outputs_list:
        etk, plast_output = entry
        if not plast_output.empty:
            if not os.path.isfile('summary_%s.plout.tsv' % etk.lower()):
                _ = plast_output.to_csv('summary_%s.plout.tsv' % etk.lower(),
                                        sep='\t', header=True)
            if seqidtaxid_map is not None:
                dfacc2tax = customtaxtranslation(seqidtaxid_map, 'Group',
                                                 plast_output)
                if dfacc2tax.empty:
                    ms = 'Something seems wrong with the acc2tax information '
                    ms += 'provided'
                    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
            else:
                dfacc2tax = Taxonomy_translation(taxDBpath, plast_output, etk,
                                                 'RefSeq')

            # --- After plast searches -------------
            # Merging all dataframes by readId using taxonomic info
            # this will fail if there is nothing to merge as it'll generate an
            # empty df
            tmpcladfs = [(label, filename, nullandmerged(df, dfacc2tax),
                          read_paths) for label, filename, df, read_paths,
                                          other in working_dfs if
                         (not dfacc2tax.empty and
                          label.startswith(etk))]
            cladfs.extend(tmpcladfs)
        else:
            ## re-check outputs... is there a bug here?
            ms = 'plastoutput empty test for unpair'
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
            tmpcladfs = [(label, filename, empty(df), read_paths)
                         for label, filename, df, read_paths, other in
                         working_dfs if label.startswith(etk)]
            ms = 'tmpcladfs\n %s' % tmpcladfs
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
            cladfs.extend(tmpcladfs)
    return cladfs


def writinggroups(cladfs, stamp, boutname):
    '''

    :param cladfs:
    :return:
    '''
    pre_args = []
    for label, prefix, df, read_paths in cladfs:
        pre_args.extend(grouper(label, prefix, df, read_paths))
    ms = '****\tWriting files by group\t****\n'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)

    for arg in pre_args:
        _ = writer(arg, stamp, boutname)
    return 'Done'


def isitready(sel_reads, fpath):
    '''

    :param sel_reads:
    :param fpath:
    :return:
    '''
    while True:
        lcount = run('wc -l %s' % fpath,
                     stderr=PIPE, stdout=PIPE, shell=True)
        lines = lcount.stdout.decode('utf-8')
        lines = int(lines.split()[0])
        if lines == len(sel_reads):
            break
        time.sleep(1)
    return 'Finished'


def writer(args, stamp, boutname):
    '''

    :param args:
    :return:
    '''
    group, label, fname, read_set, path = args
    if len(read_set) > 0:
        if isinstance(path, list) and len(path) > 0:
            rel_path = os.path.split(path[0])[0]
        else:
            rel_path = os.path.split(path)[0]
        if not 'TempEukfinder' in rel_path:
            rel_path = os.path.join(rel_path, 'TempEukfinder')
        end_path = os.path.join(os.getcwd(), '%s_%s_list.tmp' % (group, label))
        path_list = os.path.join(rel_path, end_path)
        start_time = time.time()
        handle = open(path_list, 'w')
        for read in read_set:
            handle.write('%s\n' % read)
        handle.close()
        _ = isitready(read_set, path_list)
        elapsed = round(time.time() - start_time, 3)
        minutes = round((elapsed / 60), 3)
        ms = '\nElapsed time for writing headers:\n'
        ms += '%s seconds; %s minutes' % (elapsed, minutes)
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)

        start_time = time.time()
        # patch_path0 = os.path.join(path[0], 'TempEukfinder')
        # patch_path1 = os.path.join(path[1], 'TempEukfinder')
        if label == 'Pair':
            outs = [('%s.%sR1.fq' % (fname, group), path[0]),
                    ('%s.%sR2.fq' % (fname, group), path[1])]
            for entry in outs:
                outname, path = entry
                outname = os.path.join(rel_path, outname)
                cmd = 'seqkit grep -f %s -o %s %s ' % (path_list, outname,
                                                       path)
                ms = 'seqkit grep cmd:\n%s' % cmd
                print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
                _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
                ms = 'stdout: %s\n' % _.stdout.decode('utf-8')
                elapsed = round(time.time() - start_time, 3)
                minutes = round((elapsed / 60), 3)
                ms += '\nElapsed time for writing paired outs:\n'
                ms += '%s seconds; %s minutes' % (elapsed, minutes)
                print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        else:
            # if isinstance(fname, list):
            #     fname = fname[0]
            #     if isfasta(path[0]):
            #         tail = '.fasta'
            #     else:
            #         tail = 'fq'
            #     outname = '%s.%s.%s' % (fname, group, tail)
            #     print('When in Rome:', fname, 'outname', outname)
            if fname.startswith('scf_'):
                fname = fname.split('.fasta')
                outname = '%s.%s.fasta' % (fname[0], group) ###1
            else:
                ms = 'output is not list or scaffold. The fname is %s' % fname
                print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
                if isfasta(path):
                    tail = 'fasta'
                else:
                    tail = 'fq'
                outname = '%s.%s.%s' % (fname, group, tail)

            outname_path = os.path.join(rel_path, outname)

            if 'scf_' in outname_path:
                pathched_path = os.path.join(rel_path, 'Classified_contigs')
                outname_path = os.path.join(pathched_path, outname)

            if isinstance(path, list):
                cmd = 'seqkit grep -f %s -o %s %s ' % (path_list, outname_path,
                                                       path[0])
            else:
                cmd = 'seqkit grep -f %s -o %s %s ' % (path_list, outname_path,
                                                       path)

            ms = 'seqkit grep cmd:\n%s\n' % cmd
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
            _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
            ms = 'stdout: %s\n' % _.stdout.decode('utf-8')
            elapsed = round(time.time() - start_time, 3)
            minutes = round((elapsed / 60), 3)
            ms += 'Elapsed time for writing unpaired outs:\n'
            ms += '%s seconds; %s minutes' % (elapsed, minutes)
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        # new_path = os.path.join(os.getcwd(), 'tmps_%s_%s' % (boutname, stamp))
        # shutil.move(path_list, new_path)
    return 'Done'


# --- post classification ---


def assembly(read_tuple, basename, threads, maxmem):
    ms = 'Starting assembly phase. This will take a while...'
    print(ms, sep='\t', end='\n', file=sys.stdout, flush=False)
    R1, R2, uR1R2 = read_tuple
    outdir = '%s.out' % basename
    cmd = 'metaspades.py -t %s -m %s ' % (threads, maxmem)
    cmd += '--pe1-1 %s --pe1-2 %s --pe1-s %s ' % (R1, R2, uR1R2)
    cmd += '-o %s' % outdir
    ms = 'Metaspades cmd_line:\n %s' % cmd
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    _ = run(cmd, stderr=PIPE, stdout=PIPE, shell=True)
    p1 = os.path.join(os.getcwd(), outdir)
    assembly_path = os.path.join(p1, 'scaffolds.fasta')

    if os.path.isfile(assembly_path):
        # making a copy of the assembly in the cwd
        shutil.copy(assembly_path, os.getcwd())
    else:
        # cautionary waiting time
        time.sleep(30)
        # second check
        if os.path.isfile(assembly_path):
            # making a copy of the assembly in cwd
            shutil.copy(assembly_path, os.getcwd())
        else:
            ms = "The assembly process failed. "
            ms += "There might not be enough reads left that overlap"
            ms += "to create and elongate contigs."
            ms += "Exiting program.\n"
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
            sys.exit(0)
    # rename assembly using basename
    scaffolds = change_c_names(os.path.abspath('scaffolds.fasta'),
                               basename)
    ms = 'Newly named assembly is:\n %s' % os.path.abspath(scaffolds)
    print(ms, sep='\t', end='\n', file=sys.stdout, flush=False)
    return os.path.abspath(scaffolds), scaffolds


def change_c_names(assembly, basename):
    new = {}
    track = {}
    with open(assembly) as I:
        records = I.read().split('>')[1:]
        count = 0
        for entry in records:
            count += 1
            value = entry.split('\n')
            name = 'scaffold_%s' % count
            if not name in new:
                new[name] = '\n'.join(value[1:])
                track[name] = value[0]
            else:
                ms = 'repeated accession %s' % value[0]
                print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)

    outname1 = 'scf_%s.fasta' % basename
    outname2 = 'scf_%s_cross_check.txt' % basename
    with open(outname1, 'w') as O1, \
            open(outname2, 'w') as O2:
        for key in new.keys():
            fasta = '>%s\n%s' % (key, new[key])
            check = '%s\t%s\n' % (key, track[name])
            O1.write(fasta)
            O2.write(check)
    original_scaff = os.path.abspath(assembly)
    os.remove(original_scaff)
    return outname1


def post_assembly(threads, out_name, fasta_path, dbpath, mhlen):
    # fasta_path = os.path.abspath(fasta)
    cline = 'centrifuge -f --threads %s -k 1 ' % threads
    cline += '--min-hitlen %s -x %s ' % (mhlen, dbpath)
    cline += '-U %s -S %s_scf.centrifuge_UP ' % (fasta_path, out_name)
    cline += '--report-file %s_scf.centrifuge_UP.tsv ' % out_name
    print('post assembly classification cmd:', cline)
    _ = run(cline, stdout=PIPE, stderr=PIPE, shell=True)
    report_name = '%s_scf.centrifuge_UP' % out_name
    return os.path.abspath(report_name)


def name_check(infile):
    cmd = 'head -1 %s' % infile
    head = run(cmd, stderr=PIPE, stdout=PIPE, shell=True)
    myheader = head.stdout.decode('utf-8').strip('\n')
    boolean = True if ' ' in myheader else False
    return boolean


# --- Execute ---


def Perform_prep(user_args):
    '''
    Full argument set
    :param user_args: trimmomatic arguments
    :return: all arguments to trim, map and centrifuge
    '''
    # Trimming raw reads
    ms = 'Run has started with arguments:\n'
    for key in user_args:
        ms += '%s: %s, ' % (key, user_args[key])
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)

    cen = glob.glob('%s*' % user_args['cdb'])
    if len(cen) != 4:
        ms = '\nThere is something wrong with the centrifuge database '
        ms += 'declared. Exiting program'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        sys.exit(0)

    if (not os.path.exists(user_args['r1']) or
            not os.path.exists(user_args['r2'])
            or not os.path.exists(user_args['hg']) or not
            os.path.exists(user_args['illumina_clip'])):

        ms = 'At least one of the mandatory input files do not exist in path\n'
        ms += 'Declared files are:\nReads1 %s\nReads2: %s\nHost genome: %s\n'
        ms %= user_args['r1'], user_args['r2'], user_args['hg'],
        ms += 'Adapters file: %s\n' % user_args['illumina_clip']
        ms += '\nExiting program.'
        print(ms, sep='\t', end='\n', file=sys.stdout, flush=True)
        sys.exit(0)

    # abspaths
    adapters = os.path.abspath(user_args['illumina_clip'])
    oR1 = os.path.abspath(user_args['r1'])
    oR2 = os.path.abspath(user_args['r2'])
    #
    ms = 'Check logs for additional information on run progress or errors'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    ori_host_genome = os.path.abspath(user_args['hg'])
    start_time = time.time()
    stamp = time.strftime('%Y%m%d%H%M%S', time.gmtime(start_time))
    tmpdirname = 'tmp_readprep_%s%s' % (user_args['out_name'], stamp)
    log = os.path.abspath('Read_prep_%s.log' % stamp)
    sys.stdout = open(log, 'w')
    mkdir(tmpdirname)
    tmpdirname_path = os.path.abspath(tmpdirname)
    os.chdir(tmpdirname_path)
    shutil.copy(ori_host_genome, os.getcwd())
    ms = 'Eukfinder is using python %s\n' % platform.python_version()
    ms += 'Preparing reads for analysis.'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)

    readfile_list = [(oR1, 'R1'), (oR2, 'R2')]

    nR1, nR2 = Parallel(n_jobs=-2)(delayed(rename_reads)(readfile)
                            for readfile in readfile_list)

    ms = 'Reads have been renamed'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    R1, R2, uR1R2 = trimming(user_args['out_name'], nR1, nR2, adapters,
                             user_args['wsize'], user_args['qscore'],
                             user_args['hcrop'], user_args['mlen'],
                             user_args['threads'], user_args['leading_trim'],
                             user_args['trail_trim'])
    #
    # Mapping host out
    ms = 'Getting rid of host reads...'
    print(ms, sep='\t', end='\n', file=sys.stdout, flush=False)
    b2_index = bowtie2build(user_args['hg'])
    outr1, outr2, outr1r2, cmd_line = bowtie2(user_args['out_name'],
                                              user_args['threads'], R1, R2,
                                              uR1R2, b2_index, 'fastq-host')

    ms = 'bowtie2 cmd line is:\n %s\noutr1: %s, outr2: %s and outr1r2:%s\n'
    ms %= cmd_line, outr1, outr2, outr1r2

    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    _ = run(cmd_line, stdout=PIPE, stderr=PIPE, shell=True)
    # centrifuge host-less reads
    ccmd_pline, p_report = centrifuge(user_args['out_name'], (outr1, outr2),
                                      user_args['threads'],
                                      user_args['mhlen'],
                                      user_args['cdb'], 1, pair=True,
                                      fastq=True)

    _ = run(ccmd_pline, stdout=PIPE, stderr=PIPE, shell=True)

    ccmd_upline, up_report = centrifuge(user_args['out_name'], outr1r2,
                                        user_args['threads'],
                                        user_args['mhlen'],
                                        user_args['cdb'], 1, pair=False,
                                        fastq=True)
    _ = run(ccmd_upline, stdout=PIPE, stderr=PIPE, shell=True)

    # ---- Deleting temporary files ----

    os.chdir('..')
    my_preps = (outr1, outr2, outr1r2, p_report, up_report)

    new_abspaths = []
    ms = 'Something went wrong during the read preparation.\n'
    for f in my_preps:
        try:
            shutil.move(f, os.getcwd())
            new_abspaths.append(os.path.abspath(f))
        except:
            ms += 'Exiting program.'
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
            sys.exit(1)
    try:
        os.system('rm -r %s' % tmpdirname_path)
    except:
        ms += '%s does not seem to exist\n'
        ms += 'Exiting program.'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        sys.exit(1)
    ms = "Files are ready to apply 'short_seqs' mode:\n%s" % '\n'.join(
        new_abspaths)
    print(ms, flush=True)
    return my_preps


def Perform_eukfinder(user_args):
    '''

    :param user_args:
    :return:
    '''
    # print(user_args)
    tax_up, reads_paths, redef_reads, redef_class, \
    plast_path, map_id_path, acc2tax_path, ncpu, max_plast, \
    e_value, pid, mode, base_outname, max_mem, cdb_path, aplast_path, \
    amap_id_path, cov, mhlen = input_check_and_setup(user_args)
    dirname = 'TempEukfinder'
    mkdir(dirname)
    dirpath = os.path.abspath(dirname)
    os.chdir(dirpath)
    start_time = time.time()
    stamp = time.strftime('%Y%m%d%H%M%S', time.gmtime(start_time))
    n_cpu = cpu_chunk(ncpu, max_plast)
    log = 'Class_%s.log' % stamp
    sys.stdout = open(log, 'w')
    ms = 'Eukfinder is using python %s' % platform.python_version()
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    ms = 'tmps_%s_%s' % (base_outname, stamp)
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    ms = 'user arguments are:\n %s' % user_args
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)

    # --- dumping taxonomy if requested ---
    ncbi = ete3.NCBITaxa()
    if tax_up:
        ncbi.update_taxonomy_database()

    # --- parsing centrifuge output ---
    ms = 'Reading classification reports ...\n'
    ms += 'mode is %s' % mode
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    if mode == 'short_seqs':
        if name_check(reads_paths[0]):
            rr = '\n'.join(redef_reads)
            rp = '\n'.join(reads_paths)
            ms = 'redef_reads\n%s reads paths\n%s' % (rr, rp)
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
            new_redef_reads = [rename_reads((reads_paths[0], 'R1')),
                               rename_reads((reads_paths[1], 'R2')),
                               rename_reads((reads_paths[2], 'U'))]
            ms = 'new redef reads are:\n%s\n' % '\n'.join(new_redef_reads)
            print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)

        else:
            new_redef_reads = reads_paths

        existing_dfs = pair_read_handler(n_cpu, max_plast,
                                         new_redef_reads, redef_class,
                                         new_redef_reads, stamp, base_outname)
    else:
        if name_check(redef_reads[0]):
            new_redef_reads = [rename_reads((redef_reads[0], 'U'))]
        else:
            new_redef_reads = reads_paths

        existing_dfs = single_read_handler(n_cpu, max_plast,
                                           new_redef_reads, redef_class,
                                           new_redef_reads, stamp,
                                           base_outname)
    # Management of plast searches
    if existing_dfs:
        parsed_plouts = plast_search_n_results(n_cpu, plast_path,
                                               existing_dfs, e_value, pid,
                                               cov, False)

        final_dfs = create_df_taxonomy(parsed_plouts, acc2tax_path,
                                       existing_dfs, pid,
                                       map_id_path)
        _ = writinggroups(final_dfs, stamp, base_outname)

    else:
        ms = 'No files containing reads were declared. Exiting program'
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        sys.exit(0)

    #   Meta-assembly and classification
    ms = 'Re-classification step'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    cat_reads = cats(base_outname)
    scaffolds, new_outname = assembly(cat_reads, base_outname, ncpu, max_mem)

    # Clean up -> Storing reads
    new_dir = os.path.join(os.getcwd(), 'Classified_reads')
    up_dir = os.getcwd().split('/Temp')[0]
    fqlist = glob.glob('%s/*fq' % up_dir) + glob.glob('*.fq')

    mkdir(new_dir)
    for fq in fqlist:
        shutil.move(fq, new_dir)

    # Contig classification
    sreport = post_assembly(ncpu, base_outname, scaffolds, cdb_path, mhlen)

    new_dfs = single_read_handler(n_cpu, max_plast,
                                  [scaffolds], [sreport],
                                  [scaffolds], stamp, new_outname)

    if new_dfs:
        parsed_plouts = plast_search_n_results(n_cpu, plast_path,
                                               new_dfs, e_value, pid, cov,
                                               True)

        fdf = create_df_taxonomy(parsed_plouts, acc2tax_path,
                                 new_dfs, pid,
                                 map_id_path)
        mycwd = os.getcwd()
        os.chdir('..')
        reclass_dir = 'Classified_contigs'
        mkdir(reclass_dir)
        reclass_dir_path = os.path.join(os.getcwd(), reclass_dir)
        os.chdir(mycwd)
        _ = writinggroups(fdf, stamp, new_outname)
        os.chdir('..')
        flist = glob.glob(os.path.abspath('scf_%s*' % base_outname))
        for f in flist:
            shutil.move(f, reclass_dir_path)

        # patching relacation of second centrifuge results
        centrifuge_outs = glob.glob(os.path.abspath('*centrifuge*'))
        patch_dir = os.path.abspath('tmps_%s_%s' % (new_outname, stamp))
        for out in centrifuge_outs:
            shutil.move(out, patch_dir)

        ms = "\n*****\tWARNING BEGINS\t*****\n"
        ms += "The directory 'Classified_reads' contains the read subsets "
        ms += "used to obtain an assembly.\nSuch assembly was re-classified "
        ms += "in bacteria, archaea and eukaryota/unknown files which were "
        ms += "placed in the directory 'Classified_contigs'. "
        ms += "Each of these files may contain MORE than ONE "
        ms += "eukaryotic or prokaryotic taxon or a mixture of them."
        ms += "Hence, supervised binning or "
        ms += "manual inspection on the file of your interest MUST be done.\n"
        ms += "We recommend to apply MyCC software (Lin & Liao, 2016) "
        ms += "(doi:10.1038/srep24175) or the Anvio platform (Eren et al 2015)"
        ms += "(10.7717/peerj.1319)\n"
        ms += "***\tEND OF WARNING***\n"
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)




def Perform_long_seqs(user_args):
    #
    start_time = time.time()
    stamp = time.strftime('%Y%m%d%H%M%S', time.gmtime(start_time))
    log = 'Long_seqs_%s.log' % stamp
    sys.stdout = open(log, 'w')
    ms = 'Eukfinder is using python %s' % platform.python_version()
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
    tax_up, reads_path, redef_reads, plast_path, \
    map_id_path, acc2tax_path, n_cpu, max_plast, \
    e_value, pid, base_outname, cov,  cdb_path, mhlen = input_check_and_setup(
        user_args)
    n_cpu = cpu_chunk(n_cpu, max_plast)
    new_reads = rename_reads((redef_reads[0], 'LR'))
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

    _ = run(ccmd_upline, stdout=PIPE, stderr=PIPE, shell=True)

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
        fdf = create_df_taxonomy(parsed_plouts, acc2tax_path,
                                 existing_dfs, pid, map_id_path)
        _ = writinggroups(fdf, stamp, base_outname)

    os.chdir('..')

    # removing temporary files
    '''
    try:
        shutil.rmtree(dirname_path)
        os.remove(new_reads)
    except:
        ms = 'Something went wrong when trying to remove:\n'
        ms += '%s\nExiting program' % dirname_path
        print(ms, sep=' ', end='\n', file=sys.stdout, flush=False)
        sys.exit(0)
    '''
    ms = '\n*****\tWARNING BEGINS\t*****\n'
    ms += 'The input file has been split in several files: bacteria, '
    ms += 'archaea, eukaryota and miscelaneous. Each may contain '
    ms += 'more than one eukaryotic or prokaryotic taxon.\nSupervised '
    ms += 'binning or manual inspection on the file of your interest '
    ms += 'is MANDATORY\n*****\tEND OF WARNING\t*****\n'
    print(ms, sep=' ', end='\n', file=sys.stdout, flush=True)

    return 'Done'


if __name__ == '__main__':
    import argparse

    if len(sys.argv) == 1:
        m = 'Try Eukfinder.py -h for more information'
        print(m, sep=' ', end='\n', file=sys.stdout, flush=True)
        sys.exit(1)


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
                 args.qscore, args.mlen, args.mhlen
        return bname, reads, threads, adapters, params, host


    def long_seqs(args):
        bname = args.o
        reads = args.u
        paths = args.m, args.a, args.p
        params = args.n, args.t, args.e, args.z, args.mhlen
        return reads, paths, params, bname


    myargs = {
        '-n': ['--number-of-threads', str, 'Number of threads', True],
        '-z': ['--number-of-chunks', str, 'Number of chunks to split a'
                                          ' file', True],
        '-t': ['--taxonomy-update', str, 'Set to True the first '
                                         'time the program is used. '
                                         'Otherwise set to False', True],
        '-p': ['--plast-database', str, 'path to plast database', True],
        '-m': ['--plast-id-map', str, 'path to taxonomy map for '
                                      'plast database', True],
        '-p2': ['--ancillary-plast-database', str, 'path to plast '
                                                   'database'],
        '-m2': ['--ancillary-plast-id-map', str, 'path to taxonomy map '
                                                 'for plast database'],
        '--force-pdb': ['--force_plast_database', str,
                            'impose the declared plast_database'],
        '-a': ['--acc2tax-database', str, 'path to acc2tax database', True],
        '--cdb': ['--centrifuge-database', str, 'path to centrifuge '
                                                'database', True],
        '-e': ['--e-value', float, 'threshold for plast searches', True],
        '--pid': ['--percent_id', float, 'percentage identity for '
                                         'plast searches', True],
        '--cov': ['--coverage', float, 'percentage coverage for '
                                       'plast searches', True],
        '--max_m': ['--max_memory', str, 'Maximum memomry allocated to '
                                         'carry out an assembly', True],
        '--mhlen': ['--min-hit-length', int, 'Maximum memomry allocated to '
                                         'carry out an assembly', True],
        '--pclass': ['--p-reads-class', str, 'Classification for '
                                             'pair end reads', True],
        '--uclass': ['--u-reads-class', str, 'Classification for '
                                             'un-pair end reads', True]
    }

    #  top-level parser
    parser = argparse.ArgumentParser(prog='Eukfinder')
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
                        help='out name', required=True)
    for key in myargs:
        try:
            group1.add_argument(key, myargs[key][0], type=myargs[key][1],
                                help=myargs[key][2], required=myargs[key][3])
        except:
            parser_short_seqs.add_argument(key, myargs[key][0],
                                           type=myargs[key][1],
                                           help=myargs[key][2])

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
                        help='out name', required=True)
    group2.add_argument('--cdb', '--centrifuge-database', type=str,
                        help='path to centrifuge database', required=True)

    # --- second level parser for long read mode ---
    parser_long_seqs = subparsers.add_parser('long_seqs')
    group3 = parser_long_seqs.add_argument_group('Required arguments',
                                                 'Description')
    group3.add_argument('-l', '--long-seqs', type=str,
                        help='long sequences file', required=True)
    group3.add_argument('-o', '--out_name', type=str,
                        help='out name', required=True)
    group3.add_argument('--mhlen', '--min-hit-length', type=int,
                        help='minimum hit length', required=True)
    group3.add_argument('--cdb', '--centrifuge-database', type=str,
                        help='path to centrifuge database', required=True)

    myargs_lr = {
        '-n': ['--number-of-threads', str, 'Number of threads', True],
        '-z': ['--number-of-chunks', str, 'Number of chunks to split a'
                                          ' file', True],
        '-t': ['--taxonomy-update', str, 'Set to True the first '
                                         'time the program is used. '
                                         'Otherwise set to False', True],
        '-p': ['--plast-database', str, 'path to plast database', True],
        '-m': ['--plast-id-map', str, 'path to taxonomy map for '
                                      'plast database', True],
        '-a': ['--acc2tax-database', str, 'path to acc2tax database', True],
        '-e': ['--e-value', float, 'threshold for plast searches', True],
        '--pid': ['--percent_id', float, 'percentage identity for '
                                         'plast searches', True],
        '--cov': ['--coverage', float, 'percentage coverage for '
                                       'plast searches', True],
    }

    for key in myargs_lr:
        try:
            group3.add_argument(key, myargs_lr[key][0],
                                type=myargs_lr[key][1],
                                help=myargs_lr[key][2],
                                required=myargs_lr[key][3])
        except:
            parser_long_seqs.add_argument(key, myargs_lr[key][0],
                                          type=myargs_lr[key][1],
                                          help=myargs_lr[key][2])

    parser_short_seqs.set_defaults(func=short_seqs)
    parser_read_prep.set_defaults(func=read_prep)
    parser_long_seqs.set_defaults(func=long_seqs)
    args = parser.parse_args()
    dic_args = vars(args)
    # parser.print_help()

    if dic_args['func'].__name__ == 'read_prep':
        Perform_prep(dic_args)
    if dic_args['func'].__name__ == 'short_seqs':
        Perform_eukfinder(dic_args)
    if dic_args['func'].__name__ == 'long_seqs':
        Perform_long_seqs(dic_args)

    # else:
    #     ms = 'Please select among read_prep, pair or unpair tasks'
    #     ms += '\n exiting program'
    #     print(ms, sep='\t', end='\n', file=sys.stdout, flush=False)