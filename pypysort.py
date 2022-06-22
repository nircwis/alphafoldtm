import pandas as pd
import glob
import multiprocessing
import sys
import os
import warnings
import urllib.parse
import urllib.request
import io
import argparse
import datetime
import textwrap
from pandas.core.common import SettingWithCopyWarning

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

def getpar():
    global THRESHOLD, LenghTHRESHOLD, RMSDTHRESHOLD, LaliTHRESHOLD, c
    class C:
        pass
    c=C()
    parser = argparse.ArgumentParser(
        prog='USalign sorter',
        formatter_class=argparse.RawDescriptionHelpFormatter,
       description=textwrap.dedent('''\
       This script is used to process and sort the output of the USalign screen.
       it can be used for both a singal protein results and multipul results.
       Nir Cohen 6.2022
       Enzyme Hunting SubGroup(much better than the peroxisome subgroup).
       MS Lab
       '''),
      epilog=textwrap.dedent('''\
       a normal run for a folder of proteins would be: python3 uspypysort.py -i /home/nir/Desktop/USalign/ -o /home/nir/Desktop/USalign/output.csv -t 0.5 30 2 30 -y -p -s -v -d
       a normal run for a single protein would be: python3 uspypysort.py -i /home/nir/Desktop/USalign/protein.txt -o /home/nir/Desktop/USalign/output.csv -t 0.5 30 2 30 -y -p -c
       good luck!
       '''))

    parser.add_argument('-t', '--threshold', type=float , nargs='*', help='thresholds set for the structural similarity [similarity threshold, length threshold, RMSD threshold, overlap length threshold]')
    parser.add_argument('-o', '--output', type=str, default='', help='output file name and path')
    parser.add_argument('-i', '--input', type=str, default='', help='input folder path')
    parser.add_argument('-y', '--yeast', action='store_true', help='yeast data from SGD')
    parser.add_argument('-s', '--short', action='store_true', help='make a short list of EC numbers?')
    parser.add_argument('-p', '--parallel', action='store_true', help='run in parallel?')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
    parser.add_argument('-b', '--branda', action='store_true', help='get EC numbers from brenda?')
    parser.add_argument('-u', '--update', action='store_true', help='update databases?')
    parser.add_argument('-c', '--cluster', action='store_true', help='run in cluster mode')
    parser.add_argument('-d', '--drop_peroxisome', action='store_true', help='drop all peroxisome proteins?')
    args = parser.parse_args(args=None, namespace=c)

    if c.cluster is True:
        c.branda is False
        c.short = False
        c.verbose = False

    if c.branda is False: c.short=False

    if c.input=='' and c.cluster is False:
        c.input=input('Enter input folder path: ')
        if c.input=='':
            sys.exit('No input folder path given')
    elif c.input=='' and c.cluster is True:
        sys.exit('No input folder path given')

    if c.output=='' and c.cluster is False:
        c.output=input('Enter output file name and path(if leave empty files will be saved to the input path): ')
        if c.output=='':
            if os.path.exists(c.input):
                c.output=os.path.dirname(os.path.abspath(c.input))+os.sep
            else:
                c.output=c.input
    elif c.output=='' and c.cluster is True:
        if os.path.exists(c.input):
            c.output = os.path.dirname(os.path.abspath(c.input)) + os.sep
        else:
            sys.exit('No output file name given')

    if c.threshold is None:
        THRESHOLD = 0.5
        LenghTHRESHOLD = 30
        RMSDTHRESHOLD = 2
        LaliTHRESHOLD = 30
    else:
        THRESHOLD = c.threshold[0]
        LenghTHRESHOLD = c.threshold[1]
        RMSDTHRESHOLD = c.threshold[2]
        LaliTHRESHOLD = c.threshold[3]


def extpat():
    global THRESHOLD, LenghTHRESHOLD, RMSDTHRESHOLD, LaliTHRESHOLD, c
    return THRESHOLD, LenghTHRESHOLD, RMSDTHRESHOLD, LaliTHRESHOLD, c


def update_db():
    today = datetime.datetime.today()
    scriptpath = os.path.abspath(os.path.dirname(__file__)) + os.sep
    modified_date_brnddb = datetime.datetime.fromtimestamp(os.path.getmtime(scriptpath + 'brnddb.h5'))
    modified_date_yeastprot = datetime.datetime.fromtimestamp(os.path.getmtime(scriptpath + 'yeastprot.h5'))
    modified_date_SGDDB = datetime.datetime.fromtimestamp(os.path.getmtime(scriptpath + 'SGDDB.h5'))
    duration_brnddb = today - modified_date_brnddb
    duration_yeastprot = today - modified_date_yeastprot
    duration_SGDDB = today - modified_date_SGDDB

    if os.path.exists(scriptpath + 'brnddb.h5') & os.path.exists(scriptpath + 'yeastprot.h5') & os.path.exists(scriptpath + 'SGDDB.h5') |duration_brnddb.days > 90|duration_yeastprot.days > 90|duration_SGDDB.days > 90| c.update is True: doupdate=False
    else: doupdate=True
    if doupdate is False:
        brnddb = pd.read_hdf(scriptpath + 'brnddb.h5','brnddb')
        yeastprot = pd.read_hdf(scriptpath + 'yeastprot.h5','yeastprot')
        SGDDB = pd.read_hdf(scriptpath + 'SGDDB.h5','SGDDB')

    if doupdate is True:
        brnddb = pd.read_csv('https://raw.githubusercontent.com/nircwis/alphafoldtm/main/brDB.csv', sep=',', names=['EC_PDB2', 'PDBchain2'], skiprows=1)
        yeastprot = pd.read_csv('https://raw.githubusercontent.com/nircwis/alphafoldtm/main/proteome.csv', sep=',', names=['PDBchain1', 'ORF'], skiprows=1)
        SGDDB = pd.read_csv('https://raw.githubusercontent.com/nircwis/alphafoldtm/main/SGDDB.csv', sep=',', names=['ORF', 'name', 'description'], skiprows=1)
        brnddb.to_hdf(scriptpath + 'brnddb.h5','brnddb')
        yeastprot.to_hdf(scriptpath + 'yeastprot.h5', 'yeastprot')
        SGDDB.to_hdf(scriptpath + 'SGDDB.h5', 'SGDDB')
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return brnddb, yeastprot, SGDDB


def getuniprotdatav2(listofuniprot):
    url = 'https://www.uniprot.org/uploadlists/'
    listofuniprot = listofuniprot.to_string(index=False).split('\n')
    listofuniprot = [' '.join(ele.split()) for ele in listofuniprot]
    listofuniprot = '\r\n'.join(listofuniprot)
    listofuniprot = listofuniprot.encode('utf8')
    params = {
        'from': 'ACC',
        'to': 'ACC',
        'format': 'tab',
        'columns': 'id,genes,organism,protein names,genes(PREFERRED),genes(ALTERNATIVE),genes(OLN),genes(ORF)',
        'query': listofuniprot
            }
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read().decode('utf-8').replace('\t', '^').replace(',', '_')
    result = pd.read_csv(io.StringIO(response),delimiter='^',header=0).iloc[: , :-1]
    result.rename(columns={'Entry': 'PDBchain2'}, inplace=True)
    #print(result)
    return result


def sorta(fname):
    THRESHOLD, LenghTHRESHOLD, RMSDTHRESHOLD, LaliTHRESHOLD, c = extpat()
    #read data from files filter and get DB data
    data = pd.read_csv(fname, sep='\t', names=['PDBchain1', 'PDBchain2', 'TM1', 'TM2', 'RMSD', 'ID1', 'ID2', 'IDali', 'L1', 'L2', 'Lali'], skiprows=1, usecols=['PDBchain1', 'PDBchain2', 'TM1', 'TM2', 'RMSD', 'ID1', 'ID2', 'IDali', 'L1', 'L2', 'Lali'], low_memory=False)
    data['max']=data.apply(lambda x: max(x['TM1'], x['TM2']), axis=1)
    data['minlen']=data.apply(lambda x: min(x['L1'], x['L2']), axis=1)
    data = data[data['max'] >=THRESHOLD].reset_index(drop=True)
    data = data[data['minlen'] >=LenghTHRESHOLD].reset_index(drop=True)
    data = data[data['RMSD'] <=RMSDTHRESHOLD].reset_index(drop=True)
    data = data[data['Lali'] >=LaliTHRESHOLD].reset_index(drop=True)
    if len(data) > 0:
        tp= data['PDBchain1'].str.split("-", n =2, expand = True)
        data['PDBchain1'] =tp[1]
        tp2= data['PDBchain2'].str.split("-", n =2, expand = True)
        data['PDBchain2'] =tp2[1]
        data['PDBchain2']=data['PDBchain2'].astype(str)
        data['PDBchain1']=data['PDBchain1'].astype(str)
        unidata = getuniprotdatav2(data['PDBchain2'])
        unidata['PDBchain2']=unidata['PDBchain2'].astype(str)
        data=pd.merge(unidata,data, on='PDBchain2', how='inner')
        if c.branda is True:
            brnddb['PDBchain2'] = brnddb['PDBchain2'].astype(str)
            data = pd.merge(data, brnddb, left_on='PDBchain2', right_on='PDBchain2', how='inner')
        if c.yeast is True:
            yeastprot['PDBchain1'] = yeastprot['PDBchain1'].astype(str)
            data=pd.merge(data, yeastprot, on='PDBchain1', how='inner')
    data=data.drop(['max', 'minlen'], axis=1)
    return data


def makeshortlist(newset):
    trulist = pd.DataFrame({'ORF': [], 'EC_numbers': []})
    newsetpot = newset['ORF'].unique()
    for unic in newsetpot:
        newuni = str(newset[newset['ORF'] == unic]['EC_PDB2'].unique()).replace('[', '').replace(']', '').replace("''",',').replace("'", '')
        tp = pd.DataFrame({'ORF': [unic], 'EC_numbers': [newuni]})
        trulist = pd.concat([trulist, tp])
    trulist = pd.merge(trulist, SGDDB, on='ORF', how='inner')
    return trulist

global  THRESHOLD, SGDDB, brnddb, yeastprot, LenghTHRESHOLD, RMSDTHRESHOLD, LaliTHRESHOLD, c
getpar()
brnddb,yeastprot,SGDDB = update_db()

if __name__ == '__main__':
#files location
    files = glob.glob(c.input+'*.txt')
    if c.verbose is True: print('starting with threshold of',c.threshold,'input file',c.input,'and output file', c.output)
    if c.parallel is True:
        if c.verbose is True: print('starting US align sorting with ' + str(multiprocessing.cpu_count()) + ' cores')
        with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
            data = pool.map(sorta, files)
    else:
        if c.verbose is True: print('starting US align sorting')
        data = [sorta(f) for f in files]
    if c.cluster is True:
        data=[sorta(c.input)]
        new_name=c.input.split('/')[-1].split('.')[0]+' '
    else:
        new_name=''
    if c.verbose is True : print('done with USalign sorting. saiving data')
    newset = pd.concat(data, axis=0)
    newset.to_csv(c.output + new_name + 'reduced_datasets.csv', index=False)
    if c.short is True:
        if c.verbose is True: print('runing with SGD and UNIPROT to make short list')
        trulist=makeshortlist(newset)
        if c.verbose is True: print('saiving short list')
        trulist.to_csv(c.output+'shortlist.csv', index=False, sep=',')
