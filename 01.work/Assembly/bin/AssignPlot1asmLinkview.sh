#!/usr/bin/bash
if [[ $# != 4 ]] ; then
	echo "Usage : $0 ref ref_highlight(cent/rDNA) asm asm2ref.paf"
	exit 1
fi

ref=$1
ref_highlight=$2
asm=$3
paf=$4

src=$(cd $(dirname /share/home/zhanglab/user/yangchentao/pipeline/customize/asm2ref/AssignPlot1asmLinkview.sh);pwd)
wk=`pwd`

asm_file=`basename $asm`
asm_prefix=${asm_file%.*}
PYTHON=/share/home/zhanglab/user/yangchentao/miniconda3/bin/python3

# paf.stat
#$PYTHON $src/alignmentStatFromPaf.py $paf > $paf.stat


#assign unitigs to chrs
$PYTHON $src/assignUnitig2Chrs.py $paf.stat $asm asm

cd splitbyChrs
[ -d output  ] || mkdir output
ls |grep -v 'output' | sort -k1,1V|while read chr
do
	# generate karyotype file
	cd $wk/splitbyChrs/$chr
	if [ -s $chr.asm.besthit ];then
		echo "$chr looks good"
		$PYTHON $src/getLine_by_ColumnInTab.py $paf 1 $chr.asm.besthit 1 > $chr.asm.paf
		cat $chr.asm.paf|sort -k6,6V -k8n > tmp
		direction=$(sort -k11nr tmp|head -n 1|cut -f5)
		$PYTHON /share/home/zhanglab/user/qiulingxin/projects/GP/01.polar/01.assembly/02.gapfill/03.hic_scaf/linkview/v1_ptg/MakeLinkViewKaryotypeFromPaf.py tmp $chr $direction > karyotype.txt && rm tmp
		# linkview
		if [ -s karyotype.txt ] && [ -s $chr.asm.paf ];then
			#~/software/01.draw/LINKVIEW2/node-v14.15.0-linux-x64/bin/LINKVIEW2 --style classic -k karyotype.txt --svg_width 2000 --svg_height 500 --label_font_size 10 --label_angle 30 --chro_axis --gap_length 0.01 --svg2png_dpi 600 --no_dash  -hl $ref_highlight $chr.asm.paf 
			~/software/01.draw/LINKVIEW2/node-v14.15.0-linux-x64/bin/LINKVIEW2 --style classic -k karyotype.txt --svg_width 2000 --svg_height 500 --label_font_size 15 --label_angle 30 --label_pos right --chro_axis --gap_length 0.01 --chro_axis_pos top --hl_min1px -hl $ref_highlight $chr.asm.paf
			mv linkview2_output.svg $chr.asm.svg
			mv linkview2_output.png $chr.asm.png
			cp $chr.asm.svg ../output
		else
			echo "karyotype.txt or $chr.asm.paf is empty!"
		fi
	else
		echo "$chr has no assigned sequences"
	fi
	cd $wk/splitbyChrs
done
