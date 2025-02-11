#!/usr/bin/env python
import gc
import os
import sys
import glob
import time
import shutil
import urllib3
import dask.delayed
import dask.dataframe as dd
from bs4 import BeautifulSoup
from urllib3.util import Retry
from dask import config as cfg
from subprocess import PIPE, run
from multiprocessing.dummy import Pool  # use threads for I/O bound tasks
from multiprocessing import freeze_support
from distributed import Client, LocalCluster
from urllib3.exceptions import MaxRetryError, ProtocolError


#   Info  #
__author__ = 'Dayana E. Salas-Leiva'
__email__ = 'ds2000@cam.ac.uk'
__version__ = '1.0.0'
#   End Info   #


def mkdir(dir_name):
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    else:
        print(f'Directory {dir_name} already exists. Entering to directory')



def ncbi_index(f):
    ncbi_path = 'https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/'
    htmlfile = open(f, "r")
    index = htmlfile.read()
    parse = BeautifulSoup(index, 'html.parser')
    exclude = ['Parent Directory', 'README', 'prot.accession2taxid.FULL.gz',
               'HHS Vulnerability Disclosure', 'prot.accession2taxid.gz']
    prot_list = []
    nucl_list = []
    for tag in parse.find_all('a'):
        if not tag.text.endswith('md5') and tag.text not in exclude:
            if 'prot' in tag.text or 'pdb' in tag.text:
                abs_path = os.path.join(ncbi_path, tag.text)
                prot_list.append(abs_path)
            else:
                abs_path = os.path.join(ncbi_path, tag.text)
                nucl_list.append(abs_path)
    htmlfile.close()
    return prot_list, nucl_list


def parse_acc2taxid_readme():
    base_path = 'https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/'
    cmd = f'wget {base_path} --no-check-certificate -O my_index.html'
    _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    index = os.path.join(os.getcwd(), 'my_index.html')
    prot_files, nucl_files = ncbi_index(index)
    return prot_files, nucl_files


def unzip_downloads(fpath):
    print('unzipping:', fpath)
    p = fpath.split(' ')[-1]
    cfile = os.path.split(p)
    myfile = os.path.join(os.getcwd(), cfile[-1])
    _ = run(f'gunzip {myfile}', stdout=PIPE, stderr=PIPE, shell=True)
    _ = run(f'rm {myfile}', stdout=PIPE, stderr=PIPE, shell=True)
    unzipped = myfile.replace('.gz', '')
    return unzipped


def fetch_url(url):
    filename = os.path.split(url)[-1]
    print('retrieving:', url)
    out_file = os.path.join(os.getcwd(), filename)
    # urlretrieve(url, out_file)
    http = urllib3.PoolManager()
    retry = Retry(5, raise_on_status=True,
                  status_forcelist=range(500, 600))
    while True:
        try:
            with (http.request('GET', url, retries=retry,
                               preload_content=False) as r, open(out_file,
                                                                 'wb') as f):
                shutil.copyfileobj(r, f)
        # except MaxRetryError as m_err:
        except ProtocolError:
            print(f'ProtocolError. Retrying')
            continue
        return out_file


def download_taxonomy(urls):
    # download all available files
    cwd = os.getcwd()
    fnames = [os.path.split(f)[-1].replace('.gz', '') for f in urls]
    fnames = [os.path.join(cwd, f) for f in fnames]
    # print('expected files', '\n'.join(fnames))
    file_exist = [f for f in fnames if os.path.isfile(f)
                  is True and os.stat(f).st_size != 0]
    if file_exist:
        print('pre-existent files:', '\n'.join(file_exist))
    if len(fnames) == len(file_exist):
        lst_files = '\n'.join(file_exist)
        print(f'The following files have already been downloaded '
              f'and unzipped:\n{lst_files}')
        unzipped = fnames
    else:
        all_files = Pool(len(urls)).map(fetch_url, urls)
        unzipped = Pool(len(all_files)).map(unzip_downloads, all_files)
    return unzipped


@dask.delayed
def read_file(infile):
    print('parsing:', infile)
    ddf = dd.read_csv(infile, sep='\t', header=0, dtype='object')
    return ddf


@dask.delayed
def prepare_df(ddf):
    if 'gi' in ddf.columns:
        ddf = ddf.drop('gi', axis=1)
    if 'accession' not in ddf.columns:
        ddf = ddf.assign(accession=ddf['accession.version'].apply(
                                    lambda x: x.split('.')[0],
                         meta=('accession.version', 'object')),
                         gi=ddf['accession.version'].apply(
                                    lambda x: '',
                         meta=('gi', 'object')))
    # force empty gi field
    if 'gi' not in ddf.columns:
        ddf = ddf.assign(gi=ddf['accession.version'].apply(
                                    lambda x: '',
                         meta=('gi', 'object')))
    return ddf


def df_concat(lst_dfs):
    print('concatenating')
    cats = dd.concat(lst_dfs, interleave_partitions=True)
    cats = cats.dropna(subset=['accession'])
    print('cats are ready, sorting in progress')
    cats = cats.set_index('accession', sort=True)
    return cats


def dask_workflow(files_list, f_type):
    lst = []
    for f in files_list:
        if os.path.isfile(f) and os.stat(f).st_size != 0:
            mydf = read_file(f)
            lst.append(mydf)
        else:
            print(f'{f} is missing or is empty')
            sys.exit(-1)
    prep_dfs = [prepare_df(ddf).compute() for ddf in lst]
    fdf = df_concat(prep_dfs)
    fdf.to_csv(f'acc2tax_{f_type}_all.txt', sep='\t',
               index=True, single_file=True)
    return 'Done'


def get_nodes_names():
    fpath = 'https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'
    print(f'getting nodes and names from {fpath}')
    cmd = f'wget {fpath}'
    _ = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    ofile = os.path.join(os.getcwd(), 'taxdump.tar.gz')
    _ = run(f'tar -zxvf {ofile}', stdout=PIPE, stderr=PIPE, shell=True)
    return 'Done'


def clean_up():
    keep = ['nodes.dmp', 'names.dmp', 'acc2tax_prot_all.txt',
            'acc2tax_nucl_all.txt']
    to_keep = [os.path.join(os.getcwd(), fname) for fname in keep]
    print('files to keep after clean up:', '\n'.join(to_keep))
    all_files = glob.glob(f'{os.getcwd()}/*')
    for rmf in all_files:
        if rmf not in to_keep:
            _ = run(f'rm {rmf}', stdout=PIPE, stderr=PIPE, shell=True)
    print('Workflow has ended without errors')
    return 'Done'


def host_directory(dirname):
    start_time = time.time()
    stamp = time.strftime('%d_%m_%Y', time.gmtime(start_time))
    dirname = f'{dirname}_{stamp}'
    dirname_path = os.path.abspath(dirname)
    mkdir(dirname_path)
    os.chdir(dirname_path)
    return dirname_path


def prepare(f_type, lst):
    # downloading files
    down_files = download_taxonomy(lst)
    # concatenating proteins
    _ = dask_workflow(down_files, f_type)
    gc.collect()
    return 'Done'


#  ---   execute   ---   #
if __name__ == '__main__':

    if len(sys.argv) == 1 or len(sys.argv) > 2:
        myscript = os.path.basename(__file__)
        m = (f'\nThis program will download and prepare the databases needed for'
             f'Acc2tax while using Eukfinder.\n\n'
             f'Usage:\n'
             f'{myscript} <db_type>\n'
             f'For db_type type:\n'
             f'prot to download the taxonomy for protein sequences, or\n'
             f'nucl to download the taxonomy for nucleotide sequences\n\n'
             f'*** IMPORTANT ***\n\n'
             f'Both db_types MUST be downloaded\n\n**************')

        print(m)
        sys.exit(0)
    if sys.argv[1] == 'nucl' or sys.argv[1] == 'prot':
        d_type = sys.argv[1]
    else:
        print('The db_type specified is incorrect.\n'
              'Please type either: prot or nucl')
        sys.exit(0)
    freeze_support()
    cluster = LocalCluster(n_workers=20)
    client = Client(cluster)
    cfg.set({'distributed.scheduler.worker-ttl': None})
    # create host directory
    _ = host_directory('Acc2Tax')
    # available files
    prot_list, nucl_list = parse_acc2taxid_readme()
    # ----------------------------------------------------------------
    # nucleotide info files
    if sys.argv[1] == 'nucl':
        npath = os.path.join(os.getcwd(), f'acc2tax_nucl_all.txt')
        if not os.path.isfile(npath):
            nucl_prep = prepare('nucl', nucl_list)
    # ----------------------------------------------------------------
    # protein files
    if sys.argv[1] == 'prot':
        ppath = os.path.join(os.getcwd(), f'acc2tax_prot_all.txt')
        if not os.path.isfile(ppath):
            prot_prep = prepare('prot', prot_list)
    # get nodes.dmp and names.dmp files from tax dump
    _ = get_nodes_names()
    _ = clean_up()
    os.chdir('..')
    