#!/usr/bin/env python
from Bio import AlignIO
import os
import pandas as pd

def clustal_to_ama(clustal_file,species_dic=None,id_prefix=""):
    ama_file = clustal_file+'.ama'
    clustal_file_renamed = clustal_file+'.species_named.aln' 
    print("Converting", clustal_file, 'to:', ama_file, 'and:', clustal_file_renamed)
    clustal_rec = AlignIO.read(clustal_file,'clustal')
    
    
    with open(ama_file,'w') as amaout:
        ama_header = "a score=0.000000 ENTRY={}\n".format(id_prefix+os.path.basename(clustal_file))
        #print(ama_header)
        amaout.write(ama_header)
        
        for ir, rec in enumerate(clustal_rec):
            name, seq = rec.id,str(rec.seq)
            if species_dic is not None:
                new_name = species_dic[name]
            else:
                new_name = name
            ungapped_seq = seq.replace('-','')
            ama_entry = ("s\t{}\t0\t{}\t+\t{}\t{}\n".format(new_name,len(ungapped_seq),len(ungapped_seq),seq))
            clustal_rec[ir].id = new_name
            amaout.write(ama_entry)
        amaout.write('\n')
    AlignIO.write(clustal_rec,clustal_file_renamed,format='clustal')
    return ama_file, clustal_file_renamed

def get_species_id_dic(clusterall_file):
    df_clusterall = pd.read_csv(clusterall_file,delim_whitespace=True,
                                names=['id','RESULTH','RESULT','CM_SCOREH','CMSCORE','CMSEARCHH',
                                       'CMSEARCH','ORIGIDH','ORIGID','ORIGHEADH','ORIGH']
                               )
    df_clusterall['id'] = df_clusterall['id'].str.replace('#','_')
    df_clusterall.set_index('id',inplace=True)
    return df_clusterall['ORIGH'].to_dict()

import sys
print(sys.argv, len(sys.argv))
if (len(sys.argv)!=3):
    raise RuntimeError('Expect exactly 2 argument, found {}'.format(len(sys.argv)-1))
clustal = sys.argv[1]
cluster_all = sys.argv[2]
ama, clustal_rename = clustal_to_ama(clustal, get_species_id_dic(cluster_all))
print('{} is converted to {}'.format(clustal, ama))