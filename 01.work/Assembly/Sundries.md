scaffold annotation
```bash
exonerate \
--model est2genome \
--softmasktarget \
--showtargetgff \
--percent 70 \
--geneseed 200 \
--minintron 50 \
--maxintron 2000 \
--fsmmemory 1000 \
"gpT2T_annot_v2.0_CDS.fa" \
"scf.fa" \
> "exonerate_scf_cds.gff"

awk 'NF==19 && $3~/exon/ && NR>1' exonerate_scf_cds.gff > exonerate_scf_cds_exon.txt
sort -k1,1V -k4,4n exonerate_scf_cds_exon.txt|cut -f1,4-5 > exonerate_scf_cds_exon.bed
awk 'BEGIN{OFS="\t"}{print $1,1,$2}' scf.fa.fai > scf.fa.fai.bed
bedtools coverage -a scf.fa.fai.bed -b exonerate_scf_cds_exon.bed>exonerate_scf_cds_exon.bed.coverage
```

haphic
```bash
bwa=/share/home/zhanglab/user/yangchentao/software/Aligner/bwa/bwa
$bwa index simplex_split_rename_filt.fa
#$bwa mem -5SP -t 28 simplex_split_rename_filt.fa BJX-xin_Hic.1.clean.fastq.gz BJX-xin_Hic.2.clean.fastq.gz | samblaster | samtools view - -@ 48 -S -h -b -F 3340 -o HiC.bam
$bwa mem -5SP -t 28 simplex_split_rename_filt.fa BJX-xin_Hic.1.clean.fastq.gz BJX-xin_Hic.2.clean.fastq.gz  > hic_bwa.sam
samblaster -i hic_bwa.sam | samtools view - -@ 48 -S -h -b -F 3340 -o HiC.bam
filter_bam HiC.bam 1 --nm 3 --threads 48 | samtools view - -b -@ 48 -o HiC.filtered.bam
haphic pipeline simplex_split_rename_filt.fa HiC.filtered.bam 38 --threads 4 --processes 8 --correct_nrounds 2 # default:RE="GATC"
```

break seq by gci
```bash
bedtools maskfasta -fi simplex_gci/PolarBear.hifiasm_simplex.fasta -bed <(cat simplex_gci/GCI_hifi.0.depth.bed simplex_gci/GCI_nano.0.depth.bed |sort -k1,1V -k2,2n) -fo ../masked_lowDP_simplex.fa
awk '$5~/W/' ragtag.splitasm.agp|cut -f1,6 > ptg_seq.txt
awk 'BEGIN{OFS="\t"} NR==FNR { lines[$1]=$2; next  } NR!=FNR{ if ($1 in lines) $1=lines[$1]; print $1,$2  }' rename_1.txt ptg_seq.txt > updated_file.txt
paste <(cut -f2 updated_file.txt) <(cut -f1 updated_file.txt -d' '|awk '{counts[$1]++;print $1"_"counts[$1]}')  > rename_2.txt
seqkit replace --ignore-case --kv-file rename_2.txt --pattern "^(\w+)" --replacement "{kv}" masked_lowDP_simplex_split.fa -o simplex_split_rename.fa
seqkit grep simplex_split_rename.fa -f f > simplex_split_rename_filt.fa
# change to lines
awk '/^>/ {if (NR > 1) {print quality_string} print /bin/zsh; quality_string = ; next} {quality_string = quality_string sprintf(%c, ( + 33))} END {print quality_string}' ontASM1_nano.depth > ontASM1_nano.depth.lines
awk '/^>/ {if (NR > 1) {print quality_string} print /bin/zsh; quality_string = ; next} {quality_string = quality_string sprintf(%c, ( + 33))} END {print quality_string}' ontASM1_hifi.depth > ontASM1_hifi.depth.lines
```

split ragtag
```bash
# patch for each chromosome
for i in {1..35} X Y;do
        mkdir -p /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/02.gapfill/04.patch_chr/chr${i}
        mkdir -p /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/02.gapfill/04.patch_chr/01.out_patch
        echo "#!/bin/bash
#SBATCH --job-name=chrRAG
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=cpu64
#SBATCH --mem=30g
ragtag.py patch -o chr${i} -t 48 --aligner minimap2 chr${i}_with_scaffold.fa /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/01.ont/PolarBear.hifiasm_ul_hic.asm.hic.p_ctg.fasta
ln -s /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/02.gapfill/04.patch_chr/chr${i}/ragtag.patch.fasta /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/02.gapfill/04.patch_chr/01.out_patch/chr${i}_ragtag.fa" > job_chr${i}.sh
        sbatch job_chr${i}.sh
done
```

```bash
# split fasta
draft=/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/00.draft_v1.0/PolarBear_scaffolds_final.renameID.fa
for i in {1..35} X Y;do
        echo "chr${i}" > chr${i}_list.txt
        cat draft_v1.0_scaffoldNAME.txt >> chr${i}_list.txt
done
for i in {1..35} X Y;do
        seqkit grep -f chr${i}_list.txt $draft > chr${i}_with_scaffold.fa
done
# patch for each chromosome
bash work.sh
# rename
bash /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/02.gapfill/01.patch/02.patch_chr/02.out_rename/log
```

`centIER.py draft_v1.0.fa -o ./centIER_out`

linkview2
```bash
minimap2 -t 16  -x asm10  *.asm.hic.p_ctg.fasta *_scaffolds_final.renameID.fa >map.paf
python3 alignmentStatFromPaf.py map.paf > map.paf.stat
# rDNA telomere gaps
cat ../rDNA/rDNA.filt.withColor.bed ../rDNA/v1_rDNA.filt.withColor.bed simplex.full.telomere.200_0.5_m300.bed v1.full.telomere.200_0.5_m300.bed draft_v1.0_gaps_withColor.bed > highlight.bed
sh AssignPlot1asmLinkview.sh *_ul_hic.asm.hic.p_ctg.fasta highlight.bed v1.fa map.paf
```

craq
```bash
craq -g PolarBear_scaffolds_final.renameID.splitchr8.fa -sms v1.0_hifi_pri.bam -t 48 -o craq_out


#src=/share/home/zhanglab/user/qiulingxin/miniconda3/bin/../src/
#ref_fa_size=/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/02.gapfill/03.hic_scaf/craq/craq_out/seq.size
#p=/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/02.gapfill/03.hic_scaf/craq/craq_out
#OUTPUT=/share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/02.gapfill/03.hic_scaf/craq/craq_plot/out_circos.pdf
echo -e "[M::worker_pipeline:: Plot CRAQ metrics]"
cd /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/02.gapfill/03.hic_scaf/craq/craq_plot/craq_out
/share/home/zhanglab/user/qiulingxin/miniconda3/bin/../src/runAQI_SMS.sh -g  PolarBear.hifiasm_ul_hic.asm.hic.p_ctg.fasta  -z seq.size   -e LRout/LR_eff.size  -C LRout/LR_putative.SE.SH -D LRout/LR_sort.depth  -r 0.75 -p 0.4 -q 0.6 -R 0.75 -P 0.4 -Q 0.6  -n 10 -s 249765 -w 1000000  -j 1 -b T -y T -x /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/02.gapfill/03.hic_scaf/craq/craq_plot/craq_out/seq.size -v F
```

purge_dup
```bash
pbcstat simplex.pri.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log

pri_asm=./PolarBear.hifiasm_ul_hic.asm.hic.p_ctg.fasta
split_fa ${pri_asm} >  $pri_asm.split
minimap2 -xasm5 -DP -t 48 $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed $pri_asm

# stat
awk 'BEGIN{FS=" ";OFS="\t"} $5>0.8{print $1,$3,$5,$6} ' PolarBear.hifiasm_ul_hic.asm.hic.p_ctg.fasta.split.self_dropSelf.paf.stat > v1_self.stat
```

`snakemake -s workflow/Snakefile_test --configfile=${p}/config/config.yaml --config sample=polarBear fasta=./PolarBear.hifiasm_ul_hic.asm.hic.p_ctg.fasta samtools_mem=3G --cores 48 --use-conda --rerun-incomplete cooler
`
yahs juicer
```bash
# yahs
genome='draft.fasta'
hic1=''
hic2=''
# 比对
bwa mem -t 24 $genome $hic1 $hic2 |samtools view -@12 -bS -h - > aligned.bam
#samtools sort -@36 -o aligned_sorted.bam aligned.bam
#samtools index aligned_sorted.bam

yahs --no-contig-ec --telo-motif "TTAGGG" -o PolarBear $genome aligned_sorted.bam

#cat output_prefix_scaffolds_final.agp
#cat output_prefix_scaffolds_final.fa

# juicer
p=/share/home/zhanglab/user/yangchentao/projects/panda/Carnivora/polar_bear/01.assembly/hifiasm_simplex/yahs/hic_maps
# generate alignment
juicer pre ${p}/../PolarBear.bin ${p}/../PolarBear_scaffolds_final.agp ${p}/../draft.fasta.fai |sort -k2,2d -k6,6d -T ./ --parallel=4 -S24G | awk 'NF' > alignments_sorted.txt.part && mv alignments_sorted.txt.part alignments_sorted.txt
# generate hic
ln -s ${p}/../PolarBear_scaffolds_final.fa .
samtools faidx PolarBear_scaffolds_final.fa
cut -f 1,2 PolarBear_scaffolds_final.fa.fai > scaffolds_final.chrom.sizes
java -jar -Xmx24G juicer_tools_1.22.01.jar pre alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes && mv out.hic.part out.hic
```

clip_map
```bash
file=$1
chr=$2
start=$3
end=$4
length=$(expr $end - $start)
echo "$chr $start $end"
#awk -v a="$chr" -v b="$start" -v c="$end" '$6==a && $8 >= b && $9 <= c' ${file} | sort -k8n > temp_${file}_0.paf
#awk 'BEGIN{OFS="\t"}{print $6,$7,$8,$9,$5,$1,$2,$3,$4,$10, $11, $12, $13, $14, $15, $16, $17}' temp_${file}_clip.paf > temp_${file}.paf
#awk -v len=$length 'BEGIN{OFS="\t"}{print $6,len,$8,$9,$5,$1,$2,$3,$4,$10, $11, $12, $13, $14, $15, $16, $17}' temp_${file}_0.paf > temp_${file}.paf
if [[ $chr == ptg*  ]]; then
        awk -v a="$chr" -v b="$start" -v c="$end" '$6==a && $8 >= b && $9 <= c' ${file} | sort -k8n > temp_${file}_0.paf
        awk -v len=$length 'BEGIN{OFS="\t"}{print $6,len,$8,$9,$5,$1,$2,$3,$4,$10, $11, $12, $13, $14, $15, $16, $17}' temp_${file}_0.paf > temp_${file}.paf
elif [[ $chr == chr*   ]]; then
        awk -v a="$chr" -v b="$start" -v c="$end" '$1==a && $3 >= b && $4 <= c' ${file} | sort -k8n > temp_${file}_0.paf
        awk -v len=$length 'BEGIN{OFS="\t"}{print $1,len,$3,$4,$5,$6,$7,$8,$9,$10, $11, $12, $13, $14, $15, $16, $17}' temp_${file}_0.paf > temp_${file}.paf
fi
python3 alignmentStatFromPaf.py temp_${file}.paf > temp_${file}.paf.stat
cat temp_${file}.paf.stat
cut -f3 temp_${file}.paf.stat|awk '{print "most likely: ",$0}'
#cut -f3 temp_${file}.paf.stat|awk 'BEGIN{FS=" "}{print $1}'| xargs -I {} grep -e {} temp_${file}_0.paf > temp_${file}_clip.paf
cut -f3 temp_${file}.paf.stat|awk 'BEGIN{FS=" "}{print $1}'| xargs -I {} grep -e {} temp_${file}.paf|sort -k3n > temp_${file}_clip.paf
head -n 1 temp_${file}_clip.paf|cut -f3|awk '{print "start:",$0}'
tail -n 1 temp_${file}_clip.paf|cut -f4|awk '{print "end:",$0}'

sort -k10nr temp_${file}_clip.paf|head -n1|awk '{print "longest_mapping:",$0}'

file2=temp_${file}_clip.paf
a=$(awk 'NR==1{print $1}' $file2)
b=$(awk 'NR==1{print $2}' $file2)
c=$(awk 'NR==1{print $6}' $file2)
d=$(awk 'NR==1{print $7}' $file2)
echo "${a}:1:${b}" > temp_${file}_kary.txt
echo "${c}:1:${d}" >> temp_${file}_kary.txt

cut -f1,3-6,8-9 temp_${file}_clip.paf|awk 'BEGIN {FS=OFS="\t"} $4 == "-" {t=$6; $6=$7; $7=t} 1'|cut -f1-3,5-7  > temp_${file}_alignments.txt
PYTHON=python3
$PYTHON /share/home/zhanglab/user/qiulingxin/software/01.draw/LINKVIEW/LINKVIEW.py -t 0 temp_${file}_alignments.txt --chro_axis --svg_height 300 -o ${chr}_${start}-${end} -k temp_${file}_kary.txt

#rm temp_${file}* ${chr}_${start}-${end}.png
```

read hifiasm graph
```bash
minimap2  -t 48 -x asm10 out_rename.fa PolarBear.hifiasm_ul_hic.asm.hic.p_utg.fa > map.paf
python3 /share/home/zhanglab/user/yangchentao/pipeline/customize/asm2ref/alignmentStatFromPaf.py map.paf > map.paf.stat
cut -f1 map.paf.stat >a
cut -f1,3 map.paf.stat|cut -d' ' -f1-4|cut -f2|awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,$4}' >b
paste a b > map.paf.stat.most.txt
awk 'BEGIN{OFS="";print ",v1_name"}{print ",",$2}' map.paf.stat.most.txt > a
paste ../p_utg.dp.csv a -d '' > p_utg.dp.v1name.csv

awk 'BEGIN{OFS="";print ",length"}{print ",",$2}' PolarBear.hifiasm_ul_hic.asm.hic.p_utg.fa.fai > a
paste p_utg.dp.v1name.csv a -d '' >b
awk 'BEGIN{OFS="\t"} {gsub(/,/, OFS); print}' b |sort -k3,3V -k4,4nr| awk '!seen[$3]++'|awk 'BEGIN{print "Name,Depth,v1_name";OFS=","}{print $1,$2,$3}' > p_utg.dp.v1name_uniq.csv
awk 'BEGIN{OFS="\t"} {gsub(/,/, OFS); print}' b |sort -k4,4nr| head -n 300 > head.csv

grep ptg000037l p_utg.dp.v1name.csv > 37l_28l.csv
grep chr8a p_utg.dp.v1name.csv >> 37l_28l.csv
grep chr8a map.paf.stat.most.txt|cut -f1 > chr8a_37l_utg.txt
grep ptg000037l map.paf.stat.most.txt|cut -f1 >> chr8a_37l_utg.txt
seqkit grep -f chr8a_37l_utg.txt -j 48 PolarBear.hifiasm_ul_hic.asm.hic.p_utg.fa > chr8a_37l_utg.fa
minimap2  -t 48 -x asm10 PolarBear.hifiasm_ul_hic.asm.hic.p_ctg.fasta chr8a_37l_utg.fa > map_chr8a.paf
bash get_most.sh map_chr8a.paf
awk 'BEGIN{print "Name,v1_name";OFS=","}{print $1,$2}' map_chr8a.paf.stat.most.txt > chr8a_v1name.csv

minimap2  -t 48 -x asm10 PolarBear.hifiasm_ul_hic.asm.hic.p_ctg.fa PolarBear.hifiasm_ul_hic.asm.hic.p_utg.fa > map_ctg.paf
bash get_most.sh map_ctg.paf
grep  -E  'ptg000003l|ptg000050l' map_ctg.paf.stat.most.txt | awk 'BEGIN{print "Name,v1_name";OFS=","}{print $1,$2}' > ctgName_03_50.csv

# 通过reads确定utg属于哪个ctg. 是消除noise最终选择的结果。确认了不存在一个utg对应多个ctg的情况，比map的结果真实
join -1 5 -2 5 <(awk '$1~/A/' ../../PolarBear.hifiasm_ul_hic.asm.hic.p_utg.noseq.gfa| sort -k5,5 ) <(awk '$1~/A/' ../../PolarBear.hifiasm_ul_hic.asm.hic.p_ctg.noseq.gfa|sort -k5,5)|cut -f3,11 -d' '|sort -k2|uniq |awk 'BEGIN{OFS=",";print "Name,ptg"}{print $1,$2}' > utgptg_link.csv
```