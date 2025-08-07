注意还要加上yah和hic中的使用
写上一些可能有用的工具但是还没有使用


QuarTeT ，TRFill，ScatTR
Snapbam igvsnap(?没有记录全名)
ScatTR 估计不太长的串联重复的拷贝数

着丝粒注释的工具
GAVISUNK 利用特异性kmer比对到基因组的工具

It's totally a draft.

align seq
```bash
bwa-mem2 index PolarBear.fa
bwa-mem2 mem -t 36 $genome $mgi1 $mgi2 |samtools view -@12 -bS -h - > aligned.bam
samtools sort -@48 -o aligned_sorted.bam aligned.bam
samtools index aligned_sorted.bam
```

check copy number by sequencing depth
	use: rDNA_depth_boxplot.R
```bash
# rDNA: chr8a:33751809-38015069 123 copies
# control: chr16:14037806-14466285
# hifi rDNA
samtools depth -a -b <(echo -e "chr8a\t33751809\t38015069") hifi.sort.bam -o rDNA_hifi.depth
echo $(($(awk '{sum += $1; count++} END {print sum/count}' <(cut -f3 rDNA_hifi.depth))*123/62.83)) # 65.168465701098199
# ngs rDNA
samtools depth -a -b <(echo -e "chr8a\t33751809\t38015069") aligned_sorted.bam -o rDNA_ngs.depth
echo $(($(awk '{sum += $1; count++} END {print sum/count}' <(cut -f3 rDNA_ngs.depth))*123/30.3432)) # 400.5311338289963
# ctrl hifi
samtools depth -a -b <(echo -e "chr16\t14037806\t14466285") hifi.sort.bam -o ctrl_hifi.depth
echo $(($(awk '{sum += $1; count++} END {print sum/count}' <(cut -f3 ctrl_hifi.depth))*1/62.83)) # 1.2219783542893523
# ctrl ngs
samtools depth -@48 -a -b <(echo -e "chr16\t14037806\t14466285") ngs.sort.bam -o ctrl_ngs.depth
echo $(($(awk '{sum += $1; count++} END {print sum/count}' <(cut -f3 ctrl_ngs.depth))*1/30.3432)) # 0.92588125181259728
# dp.txt
bedtools makewindows -b <(echo -e "chr16\t14037806\t14466285") -w 1000 > ctrl_win.bed
bedtools makewindows -b <(echo -e "chr8a\t33751809\t38015069") -w 1000 > rDNA_win.bed
for f in *.depth;do
        awk -v OFS='\t' '{print $1, $2-1, $2, $3}' ${f} > ${f}.bed
done
# bedtools map -a ctrl_win.bed -b ctrl_hifi.depth.bed -c 4 -o mean > ctrl_hifi_win_depth.txt
bedtools map -a ctrl_win.bed -b ctrl_hifi.depth.bed -c 4 -o mean | awk '{print $4*1/62.83}' > ./dp_file/ctrl_hifi
bedtools map -a rDNA_win.bed -b rDNA_hifi.depth.bed -c 4 -o mean | awk '{print $4*123/62.83}' > ./dp_file/rDNA_hifi
bedtools map -a ctrl_win.bed -b ctrl_ngs.depth.bed -c 4 -o mean | awk '{print $4*1/30.3432}' > ./dp_file/ctrl_ngs
bedtools map -a rDNA_win.bed -b rDNA_ngs.depth.bed -c 4 -o mean | awk '{print $4*123/30.3432}' > ./dp_file/rDNA_ngs

Rscript rDNA_depth_boxplot.R
```

check copy number by kmer
use all kmers: kmc
	use: estimate_rdna.py
```bash
:<<reads
declare -A seq
seq[ngs]="/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/02.polish/00.dateset/yak/BJX-xin_dna.clean.fastq"
seq[hifi]="/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/00.dataset/03.10k_acc/hifi_merge.fa"

for s in "${!seq[@]}"; do
        cd /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/rDNA_CopyNum/kmc/DB
        p_seq="${seq[$s]}"
        [ -d ${s}   ] || mkdir ${s}
        cd ${s}
        for n in 21 31 51 71; do
                echo '#!/bin/bash
#SBATCH --job-name=cnt_kmc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=cpu64,cpu128
#SBATCH --mem=36g' > ${s}_k${n}.sh
                echo "[ -d ${s}_k${n}_tmp   ] || mkdir ${s}_k${n}_tmp
#kmc -k${n} -t48 -m36 -ci1 -cs100000 ${p_seq} kmcDB_${s}_k${n} ${s}_k${n}_tmp
rm -r ${s}_k${n}_tmp
kmc_tools transform kmcDB_${s}_k${n} histogram kmcDB_${s}_k${n}.histo -cx100000
kmc_dump -ci2 -cx100000  kmcDB_${s}_k${n} kmcDB_${s}_k${n}_freq.txt" >> ${s}_k${n}.sh
                sbatch ${s}_k${n}.sh
        done
done
reads

rdna_path="/share/home/zhanglab/user/qiulingxin/script/source/rDNA_copy_estimation/bin"
rdna_fa="/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/rDNA/rDNA_unit.fa"
ctrl_fa="/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/02.polish/02.assess/BUSCO/Spectrin_repeat.fa"
for n in 21 31 51 71;do
        echo '#!/bin/bash
#SBATCH --job-name=rdna
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=cpu64,cpu128
#SBATCH --mem=36g' > rdna_k${n}.sh
        cat << EOF_SCRIPT >> rdna_k${n}.sh
source /share/home/zhanglab/user/yangchentao/miniconda3/bin/activate base
# rDNA
python3 ${rdna_path}/generateKmer.py ${rdna_fa} ${n} > ./rDNA/rdnaRef_k${n}.txt
python3 ${rdna_path}/extractKmerFreqFromDb.py ./ngs/kmcDB_ngs_k${n}_freq.txt ./rDNA/rdnaRef_k${n}.txt > ./rDNA/rdna_k${n}_ngs.txt
awk '{if (\$3/64>500) print \$1"\t500"; else print \$1"\\t"\$3/126}' ./rDNA/rdna_k${n}_ngs.txt > ./result/rdna_k${n}_ngs_plot.txt
python3 ${rdna_path}/estimate_rDNA_copy.py  ./rDNA/rdna_k${n}_ngs.txt ./result/rdna_k${n}_ngs.pdf 126
python3 ${rdna_path}/extractKmerFreqFromDb.py ./hifi/kmcDB_hifi_k${n}_freq.txt ./rDNA/rdnaRef_k${n}.txt > ./rDNA/rdna_k${n}_hifi.txt
awk '{if (\$3/64>500) print \$1"\t500"; else print \$1"\\t"\$3/147}' ./rDNA/rdna_k${n}_hifi.txt > ./result/rdna_k${n}_hifi_plot.txt
python3 ${rdna_path}/estimate_rDNA_copy.py  ./rDNA/rdna_k${n}_hifi.txt ./result/rdna_k${n}_hifi.pdf 147
EOF_SCRIPT
        sbatch rdna_k${n}.sh
        echo '#!/bin/bash
#SBATCH --job-name=rdna
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=cpu64,cpu128
#SBATCH --mem=36g' > ctrl_k${n}.sh
        cat << EOF_SCRIPT_2 >> ctrl_k${n}.sh
source /share/home/zhanglab/user/yangchentao/miniconda3/bin/activate base
# control
python3 ${rdna_path}/generateKmer.py ${ctrl_fa} ${n} > ./rDNA/ctrlRef_k${n}.txt
python3 ${rdna_path}/extractKmerFreqFromDb.py ./ngs/kmcDB_ngs_k${n}_freq.txt ./rDNA/ctrlRef_k${n}.txt > ./rDNA/ctrl_k${n}_ngs.txt
awk '{if (\$3/64>500) print \$1"\\t"500; else print \$1"\\t"\$3/126}' ./rDNA/ctrl_k${n}_ngs.txt > ./result/ctrl_k${n}_ngs_plot.txt
python3 ${rdna_path}/estimate_rDNA_copy.py ./rDNA/ctrl_k${n}_ngs.txt ./result/ctrl_k${n}_ngs.pdf 126
python3 ${rdna_path}/extractKmerFreqFromDb.py ./hifi/kmcDB_hifi_k${n}_freq.txt ./rDNA/ctrlRef_k${n}.txt > ./rDNA/ctrl_k${n}_hifi.txt
awk '{if (\$3/64>500) print \$1"\\t"500; else print \$1"\\t"\$3/146}' ./rDNA/ctrl_k${n}_hifi.txt > ./result/ctrl_k${n}_hifi_plot.txt
python3 ${rdna_path}/estimate_rDNA_copy.py ./rDNA/ctrl_k${n}_hifi.txt ./result/ctrl_k${n}_hifi.pdf 147
EOF_SCRIPT_2
        sbatch ctrl_k${n}.sh
done
```

```bash
rdna_path="/share/home/zhanglab/user/qiulingxin/script/source/rDNA_copy_estimation/bin"
rdna_fa="/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/rDNA/rDNA_unit.fa"
ctrl_fa="/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/02.polish/02.assess/BUSCO/Spectrin_repeat.fa"
for n in 21 31 51 71;do
        for seq in rdna ctrl;do
                echo '#!/bin/bash
#SBATCH --job-name=rdna
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=cpu64,cpu128
#SBATCH --mem=36g' > ${seq}_k${n}_draw.sh
                cat << EOF_SCRIPT >> ${seq}_k${n}_draw.sh
source /share/home/zhanglab/user/yangchentao/miniconda3/bin/activate base
a=$(sort -k3,3nr ./rDNA/${seq}_k${n}_ngs.txt|cut -f3|uniq -c|sort -k1,1nr|awk 'NR==3{print $2*1.5}')
python3 ${rdna_path}/estimate_rDNA_copy.py  <(awk -v aa=\$a '\$3<aa' ./rDNA/${seq}_k${n}_ngs.txt) ./result_2/${seq}_k${n}_ngs.pdf 126
a=$(sort -k3,3nr ./rDNA/${seq}_k${n}_hifi.txt|cut -f3|uniq -c|sort -k1,1nr|awk 'NR==3{print $2*1.5}')
python3 ${rdna_path}/estimate_rDNA_copy.py  <(awk -v aa=\$a '\$3<aa' ./rDNA/${seq}_k${n}_hifi.txt) ./result_2/${seq}_k${n}_hifi.pdf 147
EOF_SCRIPT
                sbatch ${seq}_k${n}_draw.sh
        done
done

# NEB
kmc -k31 -t48 -m36 -ci1 -cs100000 -fm /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/rDNA/kmer/data/NEB.fasta kmcDB_NEB masked_tmp
kmc_tools transform kmcDB_masked histogram kmcDB_NEB.histo -cx100000
kmc_dump -ci2 -cx100000  kmcDB_NEB kmcDB_NEB_freq.txt
rdna_path="/share/home/zhanglab/user/qiulingxin/script/source/rDNA_copy_estimation/bin"
python3 ${rdna_path}/generateKmer.py /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/rDNA/kmer/data/NEB.fasta 31 > ./k31_draw/NEB/NEB_k31.txt
python3 ${rdna_path}/extractKmerFreqFromDb.py ./ngs/kmcDB_ngs_k31_freq.txt ./k31_draw/NEB/NEB_k31.txt > ./k31_draw/NEB/NEB_k31_ngs.txt
python3 ${rdna_path}/extractKmerFreqFromDb.py ./hifi/kmcDB_hifi_k31_freq.txt ./k31_draw/NEB/NEB_k31.txt > ./k31_draw/NEB/NEB_k31_hifi.txt
mm=3
cd /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/rDNA/kmer/kmc/k31_draw
a=$(sort -k3,3nr ./NEB/NEB_k31_ngs.txt|cut -f3|uniq -c|sort -k1,1nr|awk -v mm=$m 'NR==3{print $2*mm}')
python3 estimate_rdna.py  <(awk -v aa=$a '$3<aa' ./NEB/NEB_k31_ngs.txt) ./NEB_k31_ngs.pdf 126 5
a=$(sort -k3,3nr ./NEB/NEB_k31_hifi.txt|cut -f3|uniq -c|sort -k1,1nr|awk -v mm=$m 'NR==3{print $2*mm}')
python3 estimate_rdna.py  <(awk -v aa=$a '$3<aa' ./NEB/NEB_k31_hifi.txt) ./NEB_k31_hifi.pdf 63 5
```

```bash
rdna_path="/share/home/zhanglab/user/qiulingxin/script/source/rDNA_copy_estimation/bin"
rdna_fa="/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/rDNA/rDNA_unit.fa"
ctrl_fa="/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/02.polish/02.assess/BUSCO/Spectrin_repeat.fa"
n=31
m=3
source /share/home/zhanglab/user/yangchentao/miniconda3/bin/activate base

a=$(sort -k3,3nr ../rDNA/ctrl_k${n}_ngs.txt|cut -f3|uniq -c|sort -k1,1nr|awk -v mm=$m 'NR==3{print $2*mm}')
python3 estimate_rdna.py  <(awk -v aa=$a '$3<aa' ../rDNA/ctrl_k${n}_ngs.txt) ./ctrl_k${n}_ngs.pdf 126 5
a=$(sort -k3,3nr ../rDNA/ctrl_k${n}_hifi.txt|cut -f3|uniq -c|sort -k1,1nr|awk -v mm=$m 'NR==3{print $2*mm}')
python3 estimate_rdna.py  <(awk -v aa=$a '$3<aa' ../rDNA/ctrl_k${n}_hifi.txt) ./ctrl_k${n}_hifi.pdf 63 5

a=$(sort -k3,3nr ../rDNA/rdna_k${n}_ngs.txt|cut -f3|uniq -c|sort -k1,1nr|awk -v mm=$m 'NR==3{print $2*mm}')
python3 estimate_rdna.py  <(awk -v aa=$a '$3<aa' ../rDNA/rdna_k${n}_ngs.txt) ./rdna_k${n}_ngs.pdf 126 400
a=$(sort -k3,3nr ../rDNA/rdna_k${n}_hifi.txt|cut -f3|uniq -c|sort -k1,1nr|awk -v mm=$m 'NR==3{print $2*mm}')
python3 estimate_rdna.py  <(awk -v aa=$a '$3<aa' ../rDNA/rdna_k${n}_hifi.txt) ./rdna_k${n}_hifi.pdf 63 150

convert ctrl_k31_hifi.pdf rdna_k31_hifi.pdf ctrl_k31_ngs.pdf rdna_k31_ngs.pdf merge_3.pdf
```

use all kmers: kat
```bash
declare -A region
declare -A seq
seq[ngs]="/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/02.polish/00.dateset/yak/BJX-xin_dna.clean.fastq"
seq[hifi]="/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/00.dataset/03.10k_acc/hifi_merge.fa"
region[rDNA]="/share/home/zhanglab/user/yangchentao/projects/panda/Carnivora/polar_bear/01.assembly/hifiasm_simplex/rDNA/one_copy.polarBear.rDNA.fa"
region[NEB]="/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/rDNA_CopyNum/data/NEB.fasta"

for s in "${!seq[@]}"; do
        for r in "${!region[@]}"; do
                p_seq="${seq[$s]}"
                p_region="${region[$r]}"
                echo '#!/bin/bash
#SBATCH --job-name=sect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=cpu64,cpu128
#SBATCH --mem=150g
source /share/home/zhanglab/user/qiulingxin/miniconda3/bin/activate kat' > ${r}_${s}_sect.sh
                echo "kat sect -o ${r}_${s}_sect -t 48 -m 21 ${p_region} ${p_seq}
kat plot profile ${r}_${s}_sect-counts.cvg -o ${r}_${s}_sect-counts.png" >> ${r}_${s}_sect.sh
                #sbatch ${r}_${s}_sect.sh
        done
done
```

use rDNA uniq kmers
1 get rDNA_specific_k31
```bash
makeblastdb -in pb.fa -dbtype nucl -out db
blastn -task blastn -db ./db -query rDNA_unit.fa  -out blastn.tab -evalue 1e-05 -max_target_seqs 1000000000 -num_threads 64 -outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qlen qstart qend slen sstart send evalue score bitscore sstrand stitle" -max_hsps 1000000000
rm db*
awk '$4>20000' blastn.tab|sort -k1,1V -k12,12n |cut -f2,12-13 | awk -v OFS='\t' '{ if ($2 <= $3) print $1, $2, $3; else print $1, $3, $2  }'| bedtools merge -d 2 -i - > rDNA_filtered.bed
awk '$4>20000' blastn.tab|sort -k1,1V -k12,12n  > blastn_filtered.tab
bedtools maskfasta -fi pb.fa -bed rDNA_filtered.bed -fo pb_masked.fa

# get masked kmer
[ -d masked_tmp    ] || mkdir masked_tmp
kmc -k31 -t48 -m36 -ci1 -cs100000 -fm /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/rDNA/find_rDNA/pb_masked.fa kmcDB_masked masked_tmp
kmc_tools transform kmcDB_masked histogram kmcDB_masked.histo -cx100000
kmc_dump -ci2 -cx100000  kmcDB_masked kmcDB_masked_freq.txt
rm -r masked_tmp
# get rDNA kmer
[ -d masked_tmp    ] || mkdir masked_tmp
kmc -k31 -t48 -m36 -ci1 -cs100000 -fm /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/rDNA/rDNA_unit.fa kmcDB_rDNA masked_tmp
rm -r masked_tmp
kmc_tools transform kmcDB_rDNA histogram kmcDB_rDNA.histo -cx100000
kmc_dump -ci2 -cx100000  kmcDB_rDNA kmcDB_rDNA_freq.txt
# get all rDNA_specific_kmers
kmc_tools simple kmcDB_rDNA kmcDB_masked kmers_subtract rDNA_specific
kmc_tools transform rDNA_specific histogram rDNA_specific.histo -cx100000
kmc_dump -ci1 -cx100000  rDNA_specific rDNA_specific_freq.txt
# kmc_dump -ci1 -cx100000  kmcDB_masked kmcDB_masked_freq_ci1.txt
# kmc_dump -ci1 -cx100000  kmcDB_rDNA kmcDB_rDNA_freqci1.txt
```
2 plot
```bash
for f in $(ls *_k31_*);do
        grep -f <(cut -f1 ../rDNA_specific_k31/rDNA_specific_freq.txt) ${f} > $(basename ${f} .txt)_filtered.txt
done
python3 estimate_rdna.py <(sort -k1,1n rdna_k31_hifi_filtered.txt|awk 'BEGIN{OFS="\t"}{print NR,$2,$3}') ./rdna_k31_hifi.pdf 63 80
python3 estimate_rdna.py  <(sort -k1,1n rdna_k31_ngs_filtered.txt|awk 'BEGIN{OFS="\t"}{print NR,$2,$3}') ./rdna_k31_ngs.pdf 126 400
convert rdna_k31_hifi.pdf rdna_k31_ngs.pdf rdna_k31.pdf
```

get rDNA consensus
```bash
# msa
cut -f2,12-13  blastn_filtered.tab | awk -v OFS='\t' '{ if ($2<=$3) print $1,$2,$3; else print $1,$3,$2   }' > rDNA.bed
samtools faidx -@48 /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/freeze/PolarBear.fa -r <(awk '{print $1":"$2"-"$3}' rDNA.bed) > rDNA.fa
mafft --thread 64 --auto rDNA.fa > rDNA_msa.fa
# consensus
n=0.5 # 0.1-0.5 all the outputs are the same
plu=$(($n*123))
trimal -in rDNA_msa.fa -out /dev/stdout -gt 0.8 | cons -sequence /dev/stdin -plurality ${plu} -name "rDNA_consensus_${n}"-outseq rDNA_consensus.fasta
# diff FDM_ref/rDNA_check/rDNA_unit.fa == 134 bp
# blastn
makeblastdb -in pb.fa -dbtype nucl -out db
blastn -task blastn -db ./db -query ../../rDNA_consensus.fasta -out blastn.tab -evalue 1e-05 -max_target_seqs 1000000000 -num_threads 64 -outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qlen qstart qend slen sstart send evalue score bitscore sstrand stitle" -max_hsps 1000000000
awk '$4>20000' blastn.tab|sort -k2,2V -k12,12n |cut -f2,12-13 | awk -v OFS='\t' '{ if ($2 <= $3) print $1, $2, $3; else print $1, $3, $2   }'| bedtools merge -d 2 -i - > rDNA_filtered.bed
awk '$4>20000' blastn.tab|sort -k2,2V -k12,12n  > blastn_filtered.tab
```

check two ends of rDNA reads
```bash
# minimap2 -d rDNA_consensus.fasta.mmi rDNA_consensus.fasta
fa=hifi_merge_30k.fa # depth=6.888
minimap2 -t 48 -x map-pb rDNA_consensus.fasta hifi_merge_30k.fa > map.paf
awk '$10/$11>0.9 && $11>20000 && ($8 < 100 || $9 > 34597)' map.paf > map_edge.paf
python ~/script/alignmentStatFromPaf.py map_edge.paf | tr ' ' '\t' | awk 'NF==7 && $5 < 0.7 && $6 > 0.9' > map_edge.paf.stat
grep -f <(cut -f1 map_edge.paf.stat) map_edge.paf  > map_edge_reads.paf
# del have space between rDNA
sort -k1,1V -k3,4n map_edge_reads.paf|cut -f1,3-4|bedtools merge -d 10 -i -|cut -f1 |sort |uniq -c |sort -k1,1nr|sed 's/^[[:space:]]*//'|awk '$1==1{print $2}' > right.id
sort -k1,1V -k3,4n map_edge_reads.paf|cut -f1,3-4|bedtools merge -d 10 -i -|grep -v right.id > reads_rDNA.bed
sort -k1,1V -k3,4n map_edge_reads.paf|awk 'BEGIN{OFS="\t"}{print $1,"1",$2}'|grep -f <(cut -f1 reads_rDNA.bed)|uniq > reads_rDNA.bed.len
# del span rDNA
bedtools subtract -a reads_rDNA.bed.len -b reads_rDNA.bed|awk '$3-$2>100' > reads_extract.bed
cut -f1 reads_extract.bed|uniq -c |sort -k1nr|sed 's/^[[:space:]]*//'|awk '$1==1{print $2}' > right.id
grep -f right.id reads_extract.bed > a && mv a reads_extract.bed
# get reads and map
bedtools getfasta -fi ${fa} -bed reads_extract.bed -fo reads_extract.fa
minimap2 -x ava-pb -t 48 reads_extract.fa reads_extract.fa > overlaps.paf
awk '$1!=$6 && $10/$11>0.9 ' overlaps.paf > overlaps_filt.paf
# python ~/script/alignmentStatFromPaf.py overlaps_filt.paf > overlaps_filt.paf.stat
# filt identity>0.9 reads
python ../alignmentStatFromPaf_avamap.py overlaps_filt.paf > overlaps_filt.paf.stat
awk '$9>0.9 && $7>0.9' overlaps_filt.paf.stat|cut -f1,3|sed 's/\([^[:space:]]*\):[^[:space:]]*/\1/g'>similar_reads.txt
python ../cluster_by_similarity.py similar_reads.txt > cluster_result.txt
cat cluster_result.txt
# use hifi reads would be better i think

# bam check
samtools view simplex_hifi.bam | awk 'BEGIN{OFS="\t"}{print $1, $3, $4, $4+length($10)-1}' > simplex_hifi.tsv
samtools view new_hifi.bam | awk 'BEGIN{OFS="\t"}{print $1, $3, $4, $4+length($10)-1}' > new_hifi.tsv

samtools view -@48 -bh -F 256 -F 4 -F 2048 simplex_hifi.bam > simplex_hifi_pri.bam
samtools view -@48 -bh -F 256 -F 4 -F 2048 new_hifi.bam > new_hifi_pri.bam
python extract_info.py new_hifi_pri.bam new_hifi.tsv
python extract_info.py simplex_hifi_pri.bam simplex_hifi.tsv

for i in 1 2;do
        grep -f cluser${i}.id simplex_hifi.tsv > simplex_hifi_${i}.tsv
        grep -f cluser${i}.id new_hifi.tsv > new_hifi_${i}.tsv
done
```