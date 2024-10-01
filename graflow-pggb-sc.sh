## settings -------------------------------------------------------------------

### clean pggb data folder (the docker is run by root)
# sudo rm -rf dock-dat

### pggb version
# pggb_version="202307130714058eaf35"
pggb_version="202409162204183c21d2"

### general variables
fasta_rep="${HOME}/data/nano-assemblies-pansn-2024"
n_threads=44
# fagz_files_sorted=('SGDref-0-genome.fa.gz' \
# 'DBVPG6765-0-genome.fa.gz' \
# 'AAC-1-genome.fa.gz' \
# 'AAC-2-genome.fa.gz' \
# 'AAR-0-genome.fa.gz' \
# 'ABA-0-genome.fa.gz' \
# 'ABH-0-genome.fa.gz' \
# 'AFI-0-genome.fa.gz' \
# 'S288C-0-genome.fa.gz')
fagz_files_sorted=('SGDref-0-genome.fa.gz' \
'AAB-0-genome.fa.gz' \
'AAC-1-genome.fa.gz' \
'AAC-2-genome.fa.gz' \
'AAR-0-genome.fa.gz' \
'ABA-0-genome.fa.gz' \
'ABH-0-genome.fa.gz' \
'ACH-0-genome.fa.gz' \
'ADE-0-genome.fa.gz' \
'ADI-0-genome.fa.gz' \
'ADM-1-genome.fa.gz' \
'ADM-2-genome.fa.gz' \
'ADQ-0-genome.fa.gz' \
'ADS-0-genome.fa.gz' \
'AEG-0-genome.fa.gz' \
'AEH-0-genome.fa.gz' \
'AEL-1-genome.fa.gz' \
'AEL-2-genome.fa.gz' \
'AFH-0-genome.fa.gz' \
'AFI-0-genome.fa.gz' \
'AGA_1a-0-genome.fa.gz' \
'AGK-0-genome.fa.gz' \
'AHG-0-genome.fa.gz' \
'AHL-0-genome.fa.gz' \
'AIC-0-genome.fa.gz' \
'AIE-0-genome.fa.gz' \
'AIF-1-genome.fa.gz' \
'AIF-2-genome.fa.gz' \
'AIG-0-genome.fa.gz' \
'AIS-1-genome.fa.gz' \
'AIS-2-genome.fa.gz' \
'AKH_1a-0-genome.fa.gz' \
'AKR-1-genome.fa.gz' \
'AKR-2-genome.fa.gz' \
'ALI-1-genome.fa.gz' \
'ALI-2-genome.fa.gz' \
'ALS_1a-0-genome.fa.gz' \
'AMH_1a-0-genome.fa.gz' \
'AMM_1a-0-genome.fa.gz' \
'ANE-0-genome.fa.gz' \
'APG-0-genome.fa.gz' \
'ASB-1-genome.fa.gz' \
'ASB-2-genome.fa.gz' \
'ASG-0-genome.fa.gz' \
'ASN-1-genome.fa.gz' \
'ASN-2-genome.fa.gz' \
'ATM_1a-0-genome.fa.gz' \
'AVB-0-genome.fa.gz' \
'AVI_1a-0-genome.fa.gz' \
'BAF-1-genome.fa.gz' \
'BAF-2-genome.fa.gz' \
'BAG_1a-0-genome.fa.gz' \
'BAH-0-genome.fa.gz' \
'BAI_1a-0-genome.fa.gz' \
'BAK_1a-0-genome.fa.gz' \
'BAL_1a-0-genome.fa.gz' \
'BAM-0-genome.fa.gz' \
'BAP_1a-0-genome.fa.gz' \
'BAQ_1a-0-genome.fa.gz' \
'BBF-1-genome.fa.gz' \
'BBF-2-genome.fa.gz' \
'BBL-0-genome.fa.gz' \
'BBM_1a-0-genome.fa.gz' \
'BCN-0-genome.fa.gz' \
'BFH-1-genome.fa.gz' \
'BFH-2-genome.fa.gz' \
'BFP_1a-0-genome.fa.gz' \
'BHH-0-genome.fa.gz' \
'BLD_1a-0-genome.fa.gz' \
'BMC_2a-0-genome.fa.gz' \
'BPG-0-genome.fa.gz' \
'BPK-1-genome.fa.gz' \
'BPK-2-genome.fa.gz' \
'CAS_1a-0-genome.fa.gz' \
'CBM-1-genome.fa.gz' \
'CBM-2-genome.fa.gz' \
'CCC_1a-0-genome.fa.gz' \
'CCQ_1a-0-genome.fa.gz' \
'CCT_1a-0-genome.fa.gz' \
'CDA-0-genome.fa.gz' \
'CDG_1a-0-genome.fa.gz' \
'CDN_1a-0-genome.fa.gz' \
'CEI_1a-0-genome.fa.gz' \
'CEL_1a-0-genome.fa.gz' \
'CEQ_1a-0-genome.fa.gz' \
'CFA-0-genome.fa.gz' \
'CFF-1-genome.fa.gz' \
'CFF-2-genome.fa.gz' \
'CIC-1-genome.fa.gz' \
'CIC-2-genome.fa.gz' \
'CIH-1-genome.fa.gz' \
'CIH-2-genome.fa.gz' \
'CKB-1-genome.fa.gz' \
'CKB-2-genome.fa.gz' \
'CLL-1-genome.fa.gz' \
'CLL-2-genome.fa.gz' \
'CLN-0-genome.fa.gz' \
'CMF-1-genome.fa.gz' \
'CMF-2-genome.fa.gz' \
'CNT-1-genome.fa.gz' \
'CNT-2-genome.fa.gz' \
'CPG_1a-0-genome.fa.gz' \
'CQS_1a-0-genome.fa.gz' \
'CRB_1a-0-genome.fa.gz' \
'DBVPG6044-0-genome.fa.gz' \
'YPS128-0-genome.fa.gz' \
'JXXY161-0-genome.fa.gz')

### the reference for calling built-in variants with vg deconstruct
ref_hyphen_hap="SGDref-0"
max_decomp=1000

### set folder variables
### for debug: cd scr && dir_base=$(dirname $PWD) && dir_full=$(echo $PWD)
dir_full=$(cd $(dirname $0) && pwd)
dir_base=$(dirname "${dir_full}")
dir_run=$(basename "${dir_base}")
dir_time="${dir_base}/scr"
dir_dock="${dir_base}/dock-dat"
dir_gra="${dir_base}/gra"
dir_dss="${dir_base}/dss"
dir_calls="${dir_base}/var/calls"
dir_anno="${HOME}/data/nano-assemblies-latest/annotations"
suff_out="final"
pggb_out_suffix="${suff_out}.gfa"
gfa_graph="gr-final.gfa"
### variables for the external resources
# shared_dir="${HOME}/ems/shared-plato"
# plato_dir="${HOME}/ems/plato/Tmp"
# shortreads_dir="${shared_dir}/phenovar-reads/short"
sr_ext="fq.gz"

pair_div=95
match_filt=47
poa_param="asm5"
poa_target_length=700,900,1100
poa_pad=0.001

### pggb: segment length and number of alignments per segment
n_aps=$(echo ${#fagz_files_sorted[@]})
seg_length=5000

### to avoid vg index (gcsa) construction
### to run out of temporary disk space at /tmp
### we use a Plato folder
# if [[ ! -d "${plato_dir}" ]]; then
#   echo "Plato is not mounted"
#   echo "Before re-runnnig try:"
#   echo "> [[ ! -d ~/ems/plato ]] && mkdir ~/ems/plato"
#   echo "> sshfs plato:/home ~/ems/plato"
#   exit 0
# fi
# export TMP="${plato_dir}"

### check short-read folder
# if [[ ! -d "${shortreads_dir}" ]]; then
#   echo "Plato's shared folder is not mounted"
#   echo "Before re-runnnig try:"
#   echo "> [[ ! -d ~/ems/shared-plato ]] && mkdir ~/ems/shared-plato"
#   echo "> sshfs plato:/Liti_Lab_Shared ~/ems/shared-plato"
#   exit 0
# fi

## clmnt ----------------------------------------------------------------------

### making folders for the multifasta
if [[ -d "${dir_dock}" ]]; then
  rm -rf "${dir_dock}"
fi
mkdir -p "${dir_dock}"

### graph folder
if [[ -d "${dir_gra}" ]]; then
  rm -rf "${dir_gra}"
fi
mkdir -p "${dir_gra}"

### distances folder
if [[ -d "${dir_dss}" ]]; then
  rm -rf "${dir_dss}"
fi
mkdir -p "${dir_dss}"

### variants folder
if [[ -d "${dir_calls}" ]]; then
  rm -rf "${dir_calls}"
fi
mkdir -p "${dir_calls}"

## run pggb -------------------------------------------------------------------

{ time {

### copy and format fastas, make all-genomes-ordered.txt, multifasta for pggb
cd "${dir_dock}"
rm -f all-genomes-ordered.txt
rm -f multi.fa*
for fasta_file in "${fagz_files_sorted[@]}"; do
  mito_id=$(basename "${fasta_file}" | cut -f 1 -d "-") # e.g. ADE
  mito_name="${mito_id}-mt-genome.fa.gz"
  mito_path=$(find "${fasta_rep}/" -name "${mito_name}")
  fasta_path=$(find "${fasta_rep}/" -name "${fasta_file}")
  asse_type=$(basename "${fasta_file}" | cut -f 2 -d "-") # e.g. 0 or 1
  if [[ -f "${mito_path}" ]] && [[ "${asse_type}" != "2" ]]; then
    zcat "${fasta_path}" "${mito_path}" > ./nuc-temp.fa
  else
    zcat "${fasta_path}" > ./nuc-temp.fa
  fi
  echo "${fasta_file}" >> all-genomes-ordered.txt
  ### ref_strain_id is e.g. ADE#1
  ref_strain_id=$(basename "${fasta_file}" | cut -f 1,2 -d "-" | sed 's|-|#|')
  awk '/^>/ { print ">'${ref_strain_id}'#"substr($1,2); } \
  $0 !~ />/ {print toupper($0)}' < ./nuc-temp.fa >> multi.fa
done

rm -f nuc-temp.fa
rm -f *genome.fa.gz
bgzip multi.fa
samtools faidx multi.fa.gz

### run pggb and call built-in variants
### with "vg deconstruct" against the paths of reference strain
### WARNING: do not use any double quote (aka ") in the docker command line
### with the version built locally
### (${USER}/pggb:${pggb_version} is the name given to the image 
### by docker build)

### run pggb: for the "locally built" docker: version 202408220711079c2a8b
ref_strain_id=$(echo "${ref_hyphen_hap}" | cut -d "-" -f 1) ### probably to be changed with versions 2024*
cd "${HOME}/tools/pggb/pggb-${pggb_version}"
docker run -it -v ${dir_dock}:/data ${USER}/pggb:${pggb_version} \
pggb -i /data/multi.fa.gz -s ${seg_length} \
-p ${pair_div} -n ${n_aps} -t ${n_threads} \
-k ${match_filt} -P ${poa_param} -O ${poa_pad} \
-G ${poa_target_length} -o /data -V ${ref_strain_id}:${max_decomp}
cd "${dir_dock}"

### run pggb: for the "locally built" docker: version 2023*
# ref_strain_id=$(echo "${ref_hyphen_hap}" | cut -d "-" -f 1)
# cd "${HOME}/tools/pggb/pggb-${pggb_version}"
# docker run -it -v ${dir_dock}:/data ${USER}/pggb:${pggb_version} \
# pggb -i /data/multi.fa.gz -s ${seg_length} \
# -p ${pair_div} -n ${n_aps} -t ${n_threads} \
# -k ${match_filt} -P ${poa_param} -O ${poa_pad} \
# -G ${poa_target_length} -o /data -V ${ref_strain_id}:#

### run pggb: for the "docker pull" version 202307130714058eaf35
# ref_strain_id=$(echo "${ref_hyphen_hap}" | cut -d "-" -f 1)
# docker run -it -v ${dir_dock}:/data ghcr.io/pangenome/pggb:${pggb_version} \
# pggb -i /data/multi.fa.gz -s ${seg_length} \
# -p ${pair_div} -n ${n_aps} -t ${n_threads} \
# -k ${match_filt} -P ${poa_param} -O ${poa_pad} \
# -G ${poa_target_length} -o /data -V ${ref_strain_id}:#

} } 2> "${dir_time}/time-pggb-nt${n_threads}.txt"

### transfer pggb output (the owner of the folder is root!)
pggb_graph=$(find . -name "*${pggb_out_suffix}")
cp "${pggb_graph}" "${dir_gra}/${gfa_graph}"
vcf_pref=$(echo "${pggb_graph}" | sed 's|.gfa||')
### move vcf files to the var folder
cp "${vcf_pref}.${ref_strain_id}.vcf" "${dir_calls}/bin-${ref_hyphen_hap}.vcf"

## graph manipulation ---------------------------------------------------------

cd "${dir_gra}"
### set variables
graph_name=$(echo "${gfa_graph}" | sed 's|\.gfa$||')
gra_base=$(echo "${gfa_graph}" | sed 's|\.gfa$|-chp|')
gra_graph=$(echo "${gfa_graph}" | sed 's|gfa$|vg|')
gra_chopped=$(echo "${gfa_graph}" | sed 's|\.gfa$|-chp.vg|')
gra_pruned=$(echo "${gfa_graph}" | sed 's|\.gfa$|-chp-prn.vg|')
xg_index=$(echo "${gra_chopped}" | sed 's|vg$|xg|')
gcsa_index=$(echo "${gra_chopped}" | sed 's|vg$|gcsa|')
# og_graph=$(echo "${gfa_graph}" | sed 's|\.gfa$|.og|')
og_graph=$(echo "${gfa_graph%%.*}.og")
### variable for graph annotation
gra_snarls=$(echo "${gfa_graph}" | sed 's|gfa$|snarls|')
gra_chp=$(echo "${gra_graph}" | sed 's|\.vg$|-chp.vg|')
gra_ann=$(echo "${gfa_graph}" | sed 's|\.gfa$|-chp-ann.vg|')
ann_gam=$(echo "${gfa_graph}" | sed 's|\.gfa$|-ann.gam|')

## odgi graph and similarity matrix -------------------------------------------

### gfa to og
odgi build -g "${gfa_graph}" -o "${og_graph}"

### similarity sample by contig
cnt_file="cnt-dss-${dir_run}.txt"
odgi similarity -i gr-final.og --distances > "${dir_dss}/${cnt_file}"

### similarity sample by sample
smp_file="smp-dss-${dir_run}.txt"
odgi similarity -i gr-final.og --distances -D "#" -p 2 \
> "${dir_dss}/${smp_file}"

### clean pack files and call files e.g. call-ACA-CFC-0.sh
# cd dir_scr
# rm -f *pack *sh

### clean augmented graphs
# cd dir_aug
# rm -f *-flt-aug.vg

echo "( : the graph is ready : )"
