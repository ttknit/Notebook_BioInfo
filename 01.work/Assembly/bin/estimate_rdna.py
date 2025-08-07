# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import seaborn as sns
from icecream import ic
import numpy as np
import pandas as pd
import sys

from Bio.SeqUtils import GC

from matplotlib import (pyplot as plt, lines)

def parse_depth(depth_input, coverage):
    kmer_copies = []
    total_count = []
    with open(depth_input) as fh:
        for row in fh:
            pos, kmerseq, count = row.strip().split()
            gc = GC(kmerseq)
            kmer_copies.append([int(pos), round(int(count) / coverage), gc],)
            total_count.append(round(int(count) / coverage))
    average = sum(total_count) / len(kmer_copies)
    df = pd.DataFrame(data=kmer_copies,
                 columns = ['pos','copies','gc'],
                 index = range(0, len(kmer_copies)))

    return df, average

from matplotlib import ticker
import matplotlib.pyplot as plt

def plot_depth(df, average, output, max_length=None, max_copy=None):
    plot_title = "estimated rDNA copies"
    plt.title(plot_title)
    fig, ax = plt.subplots(figsize = (15, 4))
    g = sns.lineplot(data=df, x="pos", y="copies")
    g.set(xlabel='Genome Position (bp)', ylabel="Copies", ylim=(0, max_copy), xlim=(0, max_length))
   '''
    sns.lineplot(data=df, x="pos", y="copies", ax=ax)
    ax.set(xlabel='Genome Position (bp)', ylabel="Copies", ylim=(0, max_copy), xlim=(0, max_length))
    ax2 = ax.twinx()
    sns.lineplot(data=df, x="pos", y="gc", ax=ax2, color='red', label='GC Content', alpha=0.7)
    ax2.set_ylabel('GC Content', color='red')
    ax2.set_ylim(0, 1)
    '''
    # add average of depth line
    # gene features
    '''
    features = {(3658,5526): '18S',
                (6597,6753): '5.8S',
                (7921,12971): '28S'}
    for region in features:
        start, end = region
        gene_name = features[region]
        ax.add_line(lines.Line2D([start], [df['copies'][start-1], 0], color='purple', ls='-', lw=0.2))
        ax.add_line(lines.Line2D([end], [df['copies'][end-1], 0], color='purple', ls='-', lw=0.2))
        x1, y1 =  (start + end) /2 , average * 1.2
        ax.text(x1, y1, gene_name, ha="center")
    '''
    #ax.axhline(y=average, color='red', linestyle='--', label=f'Average Depth ({average:.2f})')
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(base=10))  # 主刻度，每隔10个单位一个主刻度
    #ax.yaxis.set_major_locator(ticker.MultipleLocator(base=1)) 

    #plt.figure(figsize=(8, 4))
    # bbox_inches='tight'
    plt.savefig(output, format='pdf')
    plt.close()


if __name__ == '__main__':
    if len(sys.argv) < 5:
        sys.exit(f"python3 {sys.argv[0]} rdnaRef.31mer.freq_in_rawreads.txt output coverage max_copies")
    else:
        coverage = int(sys.argv[3])
        output = sys.argv[2]
        mc = int(sys.argv[4])
        kmer_copies, average = parse_depth(sys.argv[1], coverage)
        plot_depth(kmer_copies, average, output, max_length=kmer_copies['pos'].max(), max_copy=mc)
        #plot_depth(kmer_copies, average, output, max_length=14000, max_copy=300)
