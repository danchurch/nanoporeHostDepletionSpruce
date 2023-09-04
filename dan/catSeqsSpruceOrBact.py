#!/usr/bin/env/python3

import pandas as pd
import os, random
from Bio import SeqIO

os.chdir('/vol/piceaNanopore/dan/chrisHostDepletionTrial/compareAdapt')

def bestGuessBactSpruce(aligned2spruceDF, aligned2bactDF):
  Spruce2bactSet = set(aligned2spruceDF['qseqid']).intersection(set(aligned2bactDF['qseqid']))
  aligned2spruceDF.set_index(['qseqid'], inplace=True)
  aligned2bactDF.set_index(['qseqid'], inplace=True)
  aa = aligned2spruceDF.loc[list(Spruce2bactSet)]
  aa['aligned2'] = 'spruce'
  bb = aligned2bactDF.loc[list(Spruce2bactSet)]
  bb['aligned2'] = 'bact'
  cc = pd.concat([aa,bb])
  ff = pd.DataFrame(columns=['qseqid','aligned2'])
  for seq in list(Spruce2bactSet):
      dd = cc.loc[seq]
      dd.reset_index(inplace=True)
      ee = dd.iloc[dd['bitscore'].idxmax(), ][['qseqid','aligned2']]
      ff = pd.concat([ff,ee.to_frame().T], ignore_index=True)
  return(ff)

head=[ 'qseqid','sseqid','pident','length','mismatch',
  'gapopen','qstart','qend','sstart','send','evalue','bitscore' ]

noAdapt2spruce = pd.read_csv("noAdaptiveSamp_Spruce_blastn.csv", names=head)
adapt2spruce = pd.read_csv("depletionTrial_Spruce_blastn.csv", names=head)
noAdapt2bact = pd.read_csv("noAdaptiveSamp_bacterial_blastn.csv", names=head)
adapt2bact = pd.read_csv("depletionTrial_bacterial_blastn.csv", names=head)

noAdaptSeqCategorized = bestGuessBactSpruce(noAdapt2spruce, noAdapt2bact)
adaptSeqCategorized = bestGuessBactSpruce(adapt2spruce, adapt2bact)


noAdaptSeqCategorized.to_csv('noAdaptSeqCategorized.csv')
adaptSeqCategorized.to_csv('adaptSeqCategorized.csv')
