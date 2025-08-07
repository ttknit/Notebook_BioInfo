Some scripts come from [comery (Chentao Yang)](https://github.com/comery)
```bash
# asm.fa: 有gaps的chromosomes
# material.fa: 其他版本中可以和asm.fa比对上的contig + backbone中没有用在组装中的contig
# =============
# ues assembly
# ------
# Ragtag patch
ragtag.py patch asm.fa material.fa -o ragtag_output --aligner minimap2 -t 48
# ------
# GPatch
minimap2 -t 48 -ax asm5 material.fa asm.fa > map.sam
GPatch -q ./map.sam -r material.fa -b out.bam
# ------
# TGSgapcloser: ues assembly
tgsgapcloser --minmap_arg ' -x asm5'\
        --scaff  asm.fa \
        --reads  material.fa \
        --output tgs_out \
        --samtools /path/samtools  \
        --java    /path/java \
        --thread 48 --ne\
        > tgsgapcloser.log 2 > tgsgapcloser.err
# =============================
# ues reads
# TGSgapcloser: ues reads
# 非常耗内存，所以最好只提出来gap附近的reads和比对到组装上质量低的reads
# 3G的reads用了250g内存，供参考
tgsgapcloser --min_idy 0.9\
        --scaff  asm.fa \
        --reads  ont_reads.fa \
        --output tgs_out \
        --samtools /path/samtools  \
        --java    /path/java \
        --thread 48 \
        --ne \
        > tgsgapcloser.log 2 > tgsgapcloser.err
```

grep reads
```bash
perl /path/gaps /path/sv_polished.np_2nd.fasta > /path/gaps.bed
awk 'BEGIN{OFS="\t"}{print $1,$2-10000,$3+10000}' /path/gaps.bed > /path/gaps_broaden.bed
samtools view -@48 -L ./path/gaps_broaden.bed /path/ont.sort.bam |cut -f1 > /path/gap_ont.id
seqkit grep -j 48 -f <(cat /path/simplex_ont_unMap_lowMAPQ.id /path/v1.0_ont_unMap_lowMAPQ.id /path/gap_ont.id|sort |uniq) /path/ONT_10k.fa > /path/ont_reads.fa
seqkit grep -j 48 -f <(echo "chr24\nchr33\nchr8a\nchr8b\nchrY") /path/sv_polished.np_2nd.fasta > query.fa
samtools faidx -@48 query.fa
```

another grep reads script
```bash
# get unMap_lowMAPQ.id
cat << 'EOF' > process_bam.sh
#!/bin/bash
bam_file="$1"
output_file=$(basename "$bam_file" .bam)_unMap_lowMAPQ.id
samtools view "$bam_file" -@ 48 | awk '$2==4 || $5 <= 5{print $1}' > "$output_file"
EOF
chmod +x process_bam.sh
parallel -j 4 ./process_bam.sh ::: bam/*.bam
rm process_bam.sh

for f in /path/ONT_50k.fastq.gz /path/added.ONT_50k.fastq.gz;do
        pigz -dc $f|awk 'NR % 4 == 1 {print $1}' |sed 's/^@//' >> temp/ont_50k.id
done

mv v1.0_ont_unMap_lowMAPQ.id temp/v1.0_ont_all_unMap_lowMAPQ.id
comm -12 <(sort temp/ont_50k.id) <(sort temp/v1.0_ont_all_unMap_lowMAPQ.id) > v1.0_ont_unMap_lowMAPQ.id

# get MissTeloTail
p=/path/
join -1 1 -2 2 <(sort -k1 ${p}/genomeTelo.teloMiss.txt|grep chr) <(sort -k2 ${p}/rename.txt) |awk 'BEGIN{FS=" "}{print $3,$2}' > genomeTelo_newname.teloMiss.txt
grep -f <(cut -f1 -d' ' genomeTelo_newname.teloMiss.txt) ${p}/manual.agp > teloMiss.agp
join -1 1 -2 1 <(sort -k1 teloMiss.agp) <(sort -k1 genomeTelo_newname.teloMiss.txt) |cut -f6,10 -d' '|sort -k1V > origName_teloMiss.txt
rm genomeTelo_newname.teloMiss.txt teloMiss.agp origName_teloMiss.txt
join -1 1 -2 1 <(sort -k1,1 origName_teloMiss.txt) <(sort -k1,1 ${p}/comp.fa.fai)|cut -f1-3 -d' '|awk 'BEGIN{OFS="\t"}$2=="head"{print $1":1-20000"}$2=="tail"{v=$3-20000;print $1":"v"-"$3}' > missTelo.region

# get no_telomere_tail
threads=48
regions_file="/path/missTelo.region_2"
hifi="/path/hifi_merge.fa.gz"
ont_fq="/path/merge.ONT_10k.fastq.gz"
ont_fa="/path/ONT_10k.fa"
v1_ont_bam="/path/v1.0_ont_10k.bam"
simplex_ont_bam="/path/simplex_ont_10k.bam"

while IFS= read -r region; do
        chrom=$(echo "$region" | cut -d':' -f1)
        start=$(echo "$region" | cut -d':' -f2 | cut -d'-' -f1)
        end=$(echo "$region" | cut -d':' -f2 | cut -d'-' -f2)
        case "$chrom" in ptg*) bam_file_ont="$simplex_ont_bam" ;; chr*) bam_file_ont="$v1_ont_bam" ;; esac
        case "$chrom" in ptg*) out_file="simplex_tail.id" ;; chr*) out_file="v1_tail.id" ;; esac
        samtools view -@ ${threads} "$bam_file_ont" "${chrom}:${start}-${end}" | cut -f1 > $out_file
done < "$regions_file"

seqkit grep -f simplex_tail.id -j ${threads} $ont_fq > simplex_tail.fq
seqkit grep -f v1_tail.id -j ${threads} $ont_fq > v1_tail.fq

seqkit grep -f simplex_ont_10k_unMap_lowMAPQ.id -j ${threads} $ont_fq >> simplex_tail.fq
seqkit grep -f v1.0_ont_10k_unMap_lowMAPQ.id -j ${threads} $ont_fq >> v1_tail.fq
bgzip simplex_tail.fq
bgzip v1_tail.fq
```

ragtag rename for agp
```bash
python_code='
rename_dict = {}
with open("target_name.txt", "r") as rename_file:
    for line in rename_file:
        parts = line.strip().split("\t")
        if len(parts) == 3:
            _, old_name, new_name = parts
            rename_dict[old_name] = new_name
with open("ragtag.patch.agp", "r") as data_file, open("ragtag.patch_chrNAME.agp", "w") as output_file:
    for line in data_file:
        parts = line.strip().split("\t")
        if len(parts) >= 6:
            old_value = parts[5]
            if old_value in rename_dict:
                parts[5] = rename_dict[old_value]
        output_file.write("\t".join(parts) + "\n")'


awk '$1 !~/^#/{OFS="\t";print $1,$6}' ragtag.patch.ctg.agp|awk '$2 !~/100/{print}' |sort -k1V|less > target_name_1.txt
awk '{if ($1 != prev_chr) {counter = 1};identifier = $1 "." counter;print $0 "\t" identifier;prev_chr = $1;counter++}' target_name_1.txt > target_name.txt
rm target_name_1.txt
python -c "$python_code"
rm target_name.txt
```

