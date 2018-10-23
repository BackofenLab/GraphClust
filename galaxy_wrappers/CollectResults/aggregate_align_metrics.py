#!/usr/bin/env python

import pandas as pd
def filter_df_rnaz(df_rnaz, min_SCI,min_prob, min_num_seqs, max_multihit=0):
    return df_rnaz[(df_rnaz['multihit per_species']<=max_multihit)&(df_rnaz['num_seqs']>=min_num_seqs) & 
#                    (df_rnaz['Covariance contribution'].astype(float)<=0) & 
                   ((df_rnaz['Prediction']=='RNA') | 
                   ((df_rnaz['Structure conservation index']>=min_SCI) & (df_rnaz['SVM RNA-class probability']>=min_prob))              | 
        (df_rnaz['Structure conservation index']*100>=df_rnaz['Mean pairwise identity']))]
def clean_columns(df):
    return df[[#'cluster',
               'rscape_TP','score','Prediction','Structure conservation index',
               'SVM RNA-class probability','SVM decision value','multihit per_species',
               #'Mean z-score',
               'Mean pairwise identity','Covariance contribution',
               'num_seqs',
               'rscape_Found',
               'rscape_True',
               'beginPos','endPos', 'basePairCount',
               'cluster_human_loc',#'cluster_human_ids',
               'num_human_seqs_all', 'num_seqs_all',
               'cluster_bed',
               'alignment','fold'
               #'rscape_out',
               
              ] ]
# clean_columns(filter_df_rnaz(dfs_dictdup['HOTAIR']))

def add_itemRGB_to_bed(df,min_rscape_bp, min_rnaz_prob, skip_rest=False, simple_name=True, cluster_prefix='',cluster_suffix=''):
    color1 = '0,158,115'
    color1_darker = '0,118,95'
    color2 =  '0,114,178'
    color3=' 213,94,0'
    color4 = '86,180,233'
    color_blue_dark = '0,71,185'
    color_red_dark = '192,0,0'
    color_blue = color_blue_dark #color2
    color_green = color1_darker
    color_red = color3
    color_yellow = '240,228,66'
    color_pink = '204,121,167'

    new_multibeds = []
    for ir, r in df.iterrows():
        itemRGB = '255,255,255'
        
    
        if r['score']>1:
            if r['rscape_TP']>=min_rscape_bp and r['SVM RNA-class probability']>=min_rnaz_prob:
                itemRGB = color_red
            elif r['rscape_TP']>=min_rscape_bp:
                itemRGB = color_green
            #elif r['rscape_TP']>0:
            #    itemRGB = color_pink
            elif r['SVM RNA-class probability']>=min_rnaz_prob:
                itemRGB = color_blue
            else:
                itemRGB = '10,10,10'
        elif r['SVM RNA-class probability']>=min_rnaz_prob:
            if r['rscape_TP']>1:
                 itemRGB = color_pink
            else:
                 itemRGB = color_yellow
                    
    #         elif r['rscape_TP']>0:
    #             itemRGB = color2
            
            
        multibed = []
        for bedline in r['cluster_bed'].split('\n'):
            if (skip_rest and itemRGB=='255,255,255'):
                continue
            bed6  = bedline.split('\t')
            if(len(bed6)<3):
                #print("Warning: skip entry, missing bed6",bed6)
                continue
            bed9 = bed6 + [bed6[1], bed6[2], itemRGB]
            bed9[4] = str(int(r['SVM RNA-class probability']*100)) # Use Prob as bed score
            if simple_name is True:
                bed9[3] = '-'.join(bed9[3].split('-')[1:]) # remove gene name section for shorter names
                bed9[3] = bed9[3].replace('cluster-','C')
                bed9[3] = bed9[3].replace('Xtend5UTR-','')
            if(len(r['cluster_human_loc'].split(','))>1):
                bed9[3]+='-paralog'
            elif r['num_human_seqs_all'] >= 2:
                bed9[3]+= '-Hparalog'
            
            bed9[3]= cluster_prefix + bed9[3] + cluster_suffix
            multibed.append('\t'.join(bed9))
        new_multibeds.append('\n'.join(multibed))
        
    df['cluster_bed9'] = new_multibeds
        
def df_to_ucsc_track(df,min_rscape_bp, min_rnaz_prob, track_name='coords',skip_header=False,skip_rest=True, cluster_prefix='',cluster_suffix=''):
    add_itemRGB_to_bed(df,min_rnaz_prob,skip_rest,cluster_prefix=cluster_prefix,cluster_suffix=cluster_suffix)
    bed_str = ""
    if not skip_header:
        bed_str +='track name={} description="{}" visibility=3 itemRgb=On'.format(track_name.replace(' ','-'), 'GraphClust-2.0 ' + track_name)+ "\n"

    if len(df) > 0:
        bed_str += 'browser position ' + '\t'.join(df['cluster_bed'].values[0].split()[0:3]) + "\n"
    bed_str += '\n'.join([s for s in df.sort_values('SVM RNA-class probability', ascending=False)['cluster_bed9'].values if len(s)>0]) + "\n" # Skip empty beds due to skip_rest=True
    #     print('\n'.join(df['SVM RNA-class probability'].astype(str).values))
    return bed_str

# clean_columns(
def filter_by_min_seqnum(df, min_seq):
    return(df[df['num_seqs_all']>=min_seq])
def filter_evo_bpcount(df,min_bp):
    return(df[df['fold'].str.count('\(')>=min_bp])
def print_bed_simple(df):
    print('\n'.join([s for s in df.sort_values('SVM RNA-class probability', ascending=False)['cluster_bed'].values if len(s)>0])) # Skip empty beds due to skip_rest=True
    return df

import argparse
import glob
parser = argparse.ArgumentParser(description='Aggregate and filter alignment metrics of individual clusters, like the output of graphclust_align_cluster')
parser.add_argument('--exclude-spurious-structs', action='store_true')
parser.add_argument('--spurious-SCI', default=0.01)

parser.add_argument('--RNAz-prob-threshold', default=0.50)
parser.add_argument('--rscape-bp-threshold', default=2)

parser.add_argument('--min-seq-num', type=int, required=True)
parser.add_argument('--clusters-tsv-pattern', default="./*.tsv")
parser.add_argument('--filtered-tsv-out', type=str, required=False)
parser.add_argument('--bed-out', type=str, required=False)


args = parser.parse_args()

tsv_files = glob.glob(args.clusters_tsv_pattern)
if len(tsv_files)==0:
    raise RuntimeError('No tsv file found')


print("Info: {} tsv files obtained".format(len(tsv_files)))

dfs_tsvs = [pd.read_csv(f,sep='\t') for f in tsv_files]

df_tsvs = pd.concat(dfs_tsvs, sort=False)

print('Info: {} rows of cluster info retrived'.format(len(df_tsvs)))
 
df_filtered = clean_columns(df_tsvs)
df_filtered = filter_by_min_seqnum(df_filtered, args.min_seq_num)
if args.exclude_spurious_structs is True:
    df_filtered = df_filtered[df_filtered['Structure conservation index']>args.spurious_SCI]

print('Info: Clusters were reduced to {}'.format(len(df_filtered)))

if args.filtered_tsv_out:
    df_filtered.to_csv(args.filtered_tsv_out,sep='\t')

if args.bed_out:
    df_filtered_with_human_loc = df_filtered[~ (df_filtered['cluster_human_loc'].isnull())].copy()
    bed_str = df_to_ucsc_track(df_filtered_with_human_loc, args.rscape_bp_threshold,args.RNAz_prob_threshold,'GENE')
    with open(args.bed_out,'w') as outb:
        outb.write(bed_str)
    print('BED was written to file')
