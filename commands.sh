set -e
wd=$(pwd)
bd=$(dirname $(readlink -f $0))

# Help message
usage="
Usage: $(basename $0) -i illumina -1 f1 -2 f2 -n nanopore -o out

options:
  -i illumina   Illumina-only assembly file
  -1 f1         Illumina forward reads file
  -2 f2         Illumina reverse reads file
  -n nanopore   Nanopore reads file
  -d adb        antiSMASH database directory
  -o out        output directory
  -p cpu        maximum number of CPUs to use (default: all available)"

# Function that logs a message to stderr
function log {
  >&2 echo "$(basename $0): $(date +'%D %T') - $1"
}

# Print help messages when script finishes execution unexpectedly
function finish {
  exit_code=$?
  if [[ ${exit_code} -eq 2 ]]; then
    >&2 echo "$usage"
  else
    log "script execution finished"
  fi
  return
}
trap finish EXIT

# Parse arguments
cpu=$(nproc)
while [[ "$#" > 1 ]]; do
  case $1 in
    -i) illumina="$2"; shift;;
    -1) f1="$2"; shift;;
    -2) f2="$2"; shift;;
    -n) nanopore="$2"; shift;;
    -o) out="$2"; shift;;
    -d) adb="$2"; shift;;
    -p) if [ "$2" -lt 1 ]; then cpu=1
        elif [ "$2" -lt "$cpu" ]; then cpu="$2"
        fi; shift;;
  esac
  shift
done

# Verify parameters
if [ -z "$illumina" ]; then
  log "missing Illumina-only assembly file (set it with -i)"
  exit 2
fi
if [ -z "$f1" ]; then
  log "missing Illumina forward reads file (set it with -1)"
  exit 2
fi
if [ -z "$f2" ]; then
  log "missing Illumina reverse reads file (set it with -2)"
  exit 2
fi
if [ -z "$nanopore" ]; then
  log "missing Nanopore reads file (set it with -n)"
  exit 2
fi
if [ -z "$adb" ]; then
  log "missing antiSMASH database directory (set it with -d)"
  exit 2
fi
if [ -z "$out" ]; then
  log "missing output directory (set it with -o)"
  exit 2
fi

# Verify that files exist
files=( "$illumina" "$f1" "$f2" "$nanopore" )
for file in "${files[@]}"; do
  if [ ! -f "$file" ]; then
    log "file $file does not exist";
    exit 1
  fi
done

# Check that the antiSMASH database directory exists
if [ ! -f "$adb/tigrfam/TIGRFam.hmm" ]; then
  log "antiSMASH database not found in $adb, run download-antismash-databases --help"
  exit 1
fi

# Check that all required programs exist
programs=( 
  "trim_galore" "porechop" "masurca" "statswrapper.sh" "bwa"
  "datasets" "antismash" "agat_sp_keep_longest_isoform.pl"
)
for prog in "${programs[@]}"; do
  if ! command -v $prog &> /dev/null; then
    log "$prog command not found, is the pancluster environment activated?"
    exit 1
  fi
done

# Create output directory
mkdir -p "$out"

# Illumina trimming
if [ ! -f "$out/reads/illumina_1.fq.gz" ]; then
  log "trimming Illumina reads"
  trim_galore \
    -length 40 \
    --basename illumina \
    -o "$out"/reads \
    --cores 4 \
    --paired \
    --no_report_file \
    "$f1" "$f2" >> "$out/log.txt" 2>&1
  mv -v "$out/reads/illumina_val_1.fq.gz" "$out/reads/illumina_1.fq.gz" \
    >> "$out/log.txt" 2>&1
  mv -v "$out/reads/illumina_val_2.fq.gz" "$out/reads/illumina_2.fq.gz" \
    >> "$out/log.txt" 2>&1
else
  log "skipping Illumina read trimming"
fi

# Nanopore trimming
if [ ! -f "$out/reads/nanopore.fq.gz" ]; then
  log "trimming Nanopore reads"
  porechop \
    --format fastq.gz \
    -t 16 \
    -i "$nanopore" \
    -o "$out/reads/nanopore.fq.gz" >> "$out/log.txt" 2>&1
else
  log "skipping Nanopore read trimming"
fi

# MaSuRCA assembly
if [ ! -f "$out/assembly/hybrid.fa" ]; then
  log "performing hybrid assembly"
  mkdir -p "$out/assembly/masurca"
  cd "$out/assembly/masurca"
  masurca \
    -t 16 \
    -i '../../reads/illumina_1.fq.gz,../../reads/illumina_2.fq.gz' \
    -r ../../reads/nanopore.fq.gz >> "../../log.txt" 2>&1
  cd ..
  ln -s $(readlink -f masurca/*/primary.genome.scf.fasta) hybrid.fa
  cd "$wd"
else
  log "skipping hybrid assembly"
fi

# Compare assemblies
if [ ! -f "$out/assembly/stats.tsv" ]; then
  log "comparing assemblies"
  statswrapper.sh in="$illumina,$out/assembly/hybrid.fa" \
    2>> $out/log.txt 1> "$out/assembly/stats.tsv"
else
  log "skipping assembly comparison"
fi

# Download all reference genomes of Fusarium
mkdir -p "$out"/genomes/{fna,gff}
downloaded=$(
  ls $out/genomes/fna/ 2> /dev/null | while read line; do
    echo ${line%%.fna}
  done
)
missing=$(comm -23 <(sort "$bd/ref/accessions.txt") <(echo "$downloaded"))
if [ -z "$missing" ]; then
  log "skipping reference genome download"
else
  log "downloading reference genomes"
  datasets download genome accession \
    --no-progressbar \
    --include 'genome,gff3' \
    --filename "$out/ncbi_dataset.zip" $missing >> "$out/log.txt" 2>&1
  zipinfo -1 "$out/ncbi_dataset.zip" | while read file; do
    if [[ "$file" == *genomic.??? ]]; then
      name="$(basename "$(dirname $file)")"
      ext=${file##*genomic.}
      unzip -p "$out/ncbi_dataset.zip" "$file" > "$out/genomes/$ext/$name.$ext"
    fi
  done
  rm "$out/ncbi_dataset.zip"
fi

# Create symlink to F. oxysporum for reference alignment
if [[ ! -L "$out/alignment/reference.fa" ]]; then
  mkdir -p "$out/alignment"
  cd "$out/alignment"
  ln -s "$(readlink -f ../genomes/fna/GCF_013085055.1.fna)" reference.fa
  cd "$wd"
fi

# Index reference genome
if [ ! -f "$out/alignment/reference.fa.sa" ]; then
  log "indexing reference genome"
  bwa index "$out/alignment/reference.fa" >> "$out/log.txt" 2>&1
else
  log "skipping reference genome indexing"
fi

# Align assembly to reference genome
if [ ! -f "$out/alignment/alignment.sam" ]; then
  log "aligning assembly to reference genome"
  bwa mem "$out/alignment/reference.fa" "$out/assembly/hybrid.fa" \
    2>> $out/log.txt 1> "$out/alignment/alignment.sam"
else
  log "skipping reference genome alignment"
fi

# Filter gff files to keep only longest isoforms
mkdir -p "$out"/genomes/gff_filtered/
missing=$(comm -23 <(ls "$out"/genomes/gff) <(ls "$out"/genomes/gff_filtered))
if [ -z "$missing" ]; then
  log "skipping gff filtering"
else
  log "filtering gff files"
  for name in ${missing[@]}; do
    agat_sp_keep_longest_isoform.pl -gff "$out"/genomes/gff/$name \
      -o "$out"/genomes/gff_filtered/$name >> "$out/log.txt" 2>&1
    rm ${name%%.fna}.agat.log
  done
fi

# Annotate hybrid assembly with antiSMASH
mkdir -p "$out"/bgc/
if [ ! -d "$out"/bgc/hybrid/ ]; then
  log "annotating hybrid assembly with antiSMASH"
  antismash \
    -c $cpu \
    --taxon fungi \
    --output-dir "$out"/bgc/hybrid/ \
    --output-basename hybrid \
    --genefinding-tool glimmerhmm \
    --databases $adb \
    "$out"/assembly/hybrid.fa
else
  log "skipping hybrid assembly antiSMASH annotation"
fi

# Annotate all reference genomes with antiSMASH
all=$(
  ls "$out"/genomes/gff_filtered/ | while read line; do
    echo ${line%%.gff}
  done
)
finished=$(
  ls -d "$out"/bgc/GCF_* 2> /dev/null | while read line; do
    echo $(basename ${line})
  done
)
missing=$(comm -23 <(echo "$all") <(echo "$finished"))

if [ -z "$missing" ]; then
  log "skipping reference genome annotation with antiSMASH"
else
  log "annotating reference genomes with antiSMASH"
  for name in ${missing[@]}; do
    antismash \
      -c $cpu \
      --taxon fungi \
      --output-dir "$out"/bgc/$name \
      --output-basename $name \
      --genefinding-tool none \
      --databases $adb \
      --genefinding-gff3 "$out"/genomes/gff_filtered/$name.gff \
      "$out"/genomes/fna/$name.fna
  done
fi
