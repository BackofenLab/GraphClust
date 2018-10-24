#!/usr/bin/env python

import os
from subprocess import Popen, PIPE
import pandas as pd

from Bio import AlignIO, SeqIO

def call_RNAz(aln_file, use_cache=True):
    rnzout_file = aln_file+'.rnaz.out'
    if use_cache and os.path.exists(rnzout_file):
        print("Using cache", rnzout_file)
        with open (rnzout_file) as rzin:
            outstr = ''.join(rzin.readlines())
    else:
        cmd = 'RNAz --locarnate "{}"'.format(aln_file)
        print (cmd)
        p = Popen( cmd , stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if err:
            raise RuntimeError("Error in calling RNAz\n{}\n{}\n".format(out, err))
        outstr = out.decode('utf8')
        with open(rnzout_file,'w') as rzout:
            rzout.write(outstr)
    
    #print (out)
    rnaz_dict = dict()
    rnaz_dict['alignment'] = aln_file
    for line in outstr.split('\n'):
        splits = line.split(":")
        if len(splits) > 1:
            assert (len(splits)==2)
            rnaz_dict[splits[0].strip()] = splits[1].strip()
    rnaz_dict['num_seqs'],rnaz_dict['multihit per_species'] = check_multihit_per_species(aln_file)
    return rnaz_dict

def RNAz_out_to_df(aln_file, rnaz_file):
    print ('RNAz_out_to_df:',aln_file, rnaz_file)
    rnaz_dict = dict()
    rnaz_dict['alignment'] = aln_file

    with open(rnaz_file) as rnazin:
        for line in rnazin:
            print('line:',line)
            splits = line.split(":")
            if len(splits) > 1:
                assert (len(splits)==2)
                rnaz_dict[splits[0].strip()] = splits[1].strip()
    rnaz_dict['num_seqs'],rnaz_dict['multihit per_species'] = check_multihit_per_species(aln_file)
    return rnaz_dict


def check_multihit_per_species(align_file):
    align_rec = AlignIO.read(align_file,'clustal')
    align_ids = [r.id for r in align_rec]
    align_ids_fastaid = [i.split('_')[0] for i in align_ids]
    if len(align_ids_fastaid) != len(set(align_ids_fastaid)):
        print ('Warning possbily multiple hits of same species: {}\n\t{}'.format(align_file,align_ids))

    return len(align_ids_fastaid),(len(align_ids_fastaid) - len(set(align_ids_fastaid)))

def get_Evofold_df(evoout_file):
    try:
        df_evofold = pd.read_csv(evoout_file,sep='\t')
    except pd.io.common.EmptyDataError:
        print('Empty evofold out')
        return dict()
        
    
    print(df_evofold)
    if(df_evofold.empty):
        return {col:None for col in df_evofold.columns}
    return df_evofold.iloc[df_evofold['score'].idxmax()].to_dict() # If multiples, return one with max score, 1-row dataframe return
#     return df_evofold[df_evofold['score']==df_evofold['score'].max()] # If multiples, return one with max score, 1-row dataframe return

def parse_rscape_sumout(rout_file):
    print('Parse rscape:', rout_file)
    try:
        df_rscape = pd.read_csv(rout_file,delim_whitespace=True)
    except:
        print("Warning Rscape file is empty")
        return dict()
    assert(len(df_rscape)==1)
    dict_rscape = df_rscape.iloc[0].to_dict()
    dict_rscape['out'] = rout_file
    return {'rscape_{}'.format(k):dict_rscape[k] for k in dict_rscape}
    
    
def get_species_id_dic(clusterall_file):
    df_clusterall = pd.read_csv(clusterall_file,delim_whitespace=True,
                                names=['id','RESULTH','RESULT','CM_SCOREH','CMSCORE','CMSEARCHH',
                                       'CMSEARCH','ORIGIDH','ORIGID','ORIGHEADH','ORIGH']
                               )
    df_clusterall['id'] = df_clusterall['id'].str.replace('#','_')
    df_clusterall.set_index('id',inplace=True)
    return df_clusterall['ORIGH'].to_dict()


def extract_species_seq(align_file, key='SEQ1_'):
    align_rec = AlignIO.read(align_file,'clustal')
    #     align_ids = [r.id for r in align_rec]
    #     align_ids_fastaid = [idd.split('_')[0] for idd in align_ids]
    found_recs = [(r.id,str(r.seq.ungap('-'))) for r in align_rec if key in r.id]
    if (len(found_recs)>1):
        print ("WARNING multiple seq in alignment found for key ", key, align_file)
    if (len(found_recs)==0):
        print ("WARNING NO seq in alignment found for key ", key, align_file)
    return found_recs
    


def get_alignment_df(align_file,
                     rnazout_file,
                       sum_file,
                     evofold_file,
                     cluster_all_file,
                     full_transcript_seq=None,
                     full_transcript_bedstr='',
                       sepecies_assmbly_id='hg38',
                       out_tsv_file=None,
                       
                       
                      ):
    # Location human sequence on full lncRNA transcript
    if len(align_file.split('/'))>4:
        cluster_id = align_file.split('/')[3]
    else:
        cluster_id = 'X'
    if(full_transcript_seq):
        loc_dicts = cluster_locs_on_transcript(align_file,full_transcript_seq)
    else:
        loc_dicts = []
    hit_ranges = [d['location_range'] for d in loc_dicts]
    hit_ids = [d['graphclust_id'] for d in loc_dicts]
    hit_genomic_beds = [relative_location_to_genomic(rang,full_transcript_bedstr,
                                                     cluster_id) for rang in hit_ranges] 
    locations_dict= {
        'cluster_human_ids':','.join(hit_ids),
        'cluster_human_loc':','.join(hit_ranges),
        'cluster_bed':'\n'.join(hit_genomic_beds)
        }

    rnaz_result_dict = RNAz_out_to_df(align_file,rnazout_file)

    # use the evo df to add .all human and spec counts
    spec_orighead_dic = get_species_id_dic(cluster_all_file)
    num_all_seqs = len(spec_orighead_dic)
    num_human_seqs = sum( v == sepecies_assmbly_id for v in spec_orighead_dic.values() )
    rnaz_result_dict['num_seqs_all'] = num_all_seqs
    rnaz_result_dict['num_human_seqs_all'] = num_human_seqs

    print('rnaz_result_dict:',rnaz_result_dict)
    # Call and read RNAz
    df_rnaz = pd.DataFrame.from_dict([rnaz_result_dict])
    
    
  
    print(rnaz_result_dict)

    print('df_rnaz:',df_rnaz)
    #print(df_rnaz['alignment'])
    #df_rnaz['cluster']  = (df_rnaz['alignment'].str.split('/',n=4,expand=True)[3]).astype(int)   

    for col in [#'Columns',
    'Combinations/Pair', 'Consensus MFE',
       'Covariance contribution',  'Energy contribution',
       'G+C content', 'Mean pairwise identity', 'Mean single sequence MFE',
       'Mean z-score', 
       'SVM RNA-class probability', 'SVM decision value', 'Sequences',
       'Shannon entropy', 'Structure conservation index']:
        print('col:',col)
        df_rnaz[col] = df_rnaz[col].astype(float)
    
    
    # Call Evofold on 
    df_evo_single = get_Evofold_df(evofold_file)
    
    
 
    df_evofold = pd.DataFrame.from_dict([df_evo_single])
    df_rnaz_evofold = df_rnaz.join(df_evofold)

    # Parse precomputed R-scape output
    rscape_result_dict = parse_rscape_sumout(sum_file)
    df_rscape = pd.DataFrame.from_dict([rscape_result_dict])
    df_rnaz_evofold_rscape = df_rnaz_evofold.join(df_rscape)
    df_rnaz_evofold_rscape = df_rnaz_evofold_rscape.join(pd.DataFrame.from_dict([locations_dict]))
    
    
    if out_tsv_file is not None:
            df_rnaz_evofold_rscape.to_csv(out_tsv_file,sep='\t',float_format='%.3f')
    
    #return df_rnaz_evofold_rscape


def cluster_locs_on_transcript(align_file,trascript_seq):
    human_seqs = extract_species_seq(align_file)
    #     print(human_seqs)
    #     print(trascript_seq.upper().replace('T','U'))
    import re
    locations = list()
    for seqid, seq in human_seqs:  # There can be multi seqs of same species per cluster
        print(seqid, seq)
        hits_in_transcript = re.findall(seq, trascript_seq.upper().replace('T','U'))
        assert(len(hits_in_transcript)==1)
        hit_location =trascript_seq.upper().replace('T','U').find(seq) + 1
        assert(hit_location != 0) # 1 based location
#         locations.append({'graphclust_id':seqid,'start':hit_location, 'end':hit_location+len(seq)-1})
        locations.append({'graphclust_id':seqid,'location_range':'{}-{}'.format(hit_location,hit_location+len(seq)-1)})
    return locations



def relative_location_to_genomic(rel_range,transcript_bedstr,cluster_id='X',cluster_score=0):
    bed_row = bed_str_to_dict(transcript_bedstr)
    rel_start, rel_end =  rel_range.split('-')
    rel_start, rel_end = int(rel_start), int(rel_end)
    newbed_elements = [bed_row['chrom']]
    if(bed_row['strand'] == '+'):
        genom_start = bed_row['start'] + rel_start
        genom_end = bed_row['start'] + rel_end
        newbed_elements += [str(genom_start),str(genom_end)]
    elif (bed_row['strand'] == '-'):
        genom_start = bed_row['end'] - rel_start
        genom_end = bed_row['end'] - rel_end
        newbed_elements += [str(genom_end),str(genom_start)]
    else:
        raise(RuntimeError('Invalid strand in row {}'.format(bed_row)))

    newbed_elements += ['cluster-'+cluster_id,str(cluster_score),bed_row['strand']]
    return ('\t'.join(newbed_elements))

# LNCRNA_BED = './lncRNAs_hg38.bed'
# df_LNCRNA_BED = pd.read_csv(LNCRNA_BED,sep='\t',comment='#',names=['chrom','start','end','name','score','strand'])
def bed_str_to_dict(bed_str):
    splits = bed_str.split()
    if len(splits)<6:
        raise RuntimeError('Bed string for the trasccript must have at least 6 elements (chomr, start, end, name, score, strand)')
    if not (splits[1].isdigit() and splits[2].isdigit()):
        raise RuntimeError('Bed string start or end are not integer')
    if splits[5] != '+' and splits[5] != '-':
        raise RuntimeError('Bed string strand unknow expected + or - found: "{}"'.format(splits[5]))
    bed_dict = {'chrom':splits[0],'start':int(splits[1]),'end':int(splits[2]),'strand':splits[5]}
    return bed_dict


def read_transcript(transcript_fasta):
    recs = list(SeqIO.parse(transcript_fasta,'fasta'))
    if len(recs) != 1:
        raise RuntimeError('Transcript fasta must contain exactly 1 sequence, found: {}'.format(len(recs)))
    print('Transcript of length {} read'.format(len(recs[0])))
    return str(recs[0].seq)

import sys
help_str = 'Format:\nextract_conservation_metrics.py result.aln result.aln.rnaz.out ' + \
            'result.aln.sum result.aln.ama.evoout cluster.all transcript.fasta "chrom start end name 0 strand" "hg38" out.tsv' 
if (len(sys.argv)!=10):
    raise RuntimeError('Expect exactly 9 argument, found {}\n{}'.format(len(sys.argv)-1, help_str))
clustal, rnaz_out, rscape_out, evofold_out, cluster_all, full_reference_transcript, transcript_bed_str, assembly_hg, out_tsv = sys.argv[1:10]


print('Extarcting conservation metrics..')

if (len(full_reference_transcript)==0):
    full_seq = None
else:
    full_seq = read_transcript(full_reference_transcript)
get_alignment_df(clustal, rnaz_out, rscape_out, evofold_out, cluster_all, 
                full_seq, transcript_bed_str, assembly_hg, out_tsv)
                
print('Conservation metrics were extracte for', clustal)
