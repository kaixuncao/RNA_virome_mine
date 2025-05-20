#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -p kshcnormal
#AUTHOR: CAO Kaixun


# === CONFIGURATION PARAMETERS ===
# Database paths (modify these as needed)
MMSEQS_DB=/public/home/ac99cxsnzt/RdRp/RdRpdb/mmseqs_db/Cell_plam_mmseqsdb
DIAMOND_DB=/public/home/ac99cxsnzt/RdRp/RdRpdb/Cell_plam_db_diamond/Cell_plam_db
BLAST_DB=/public/home/ac99cxsnzt/RdRp/RdRpdb/Cell_plam_blastdb/Cell_plam_blastdb
HMM_DB=/public/home/ac99cxsnzt/RdRp/RdRpdb/HMMER_DB/Cell_plam_hmmdb.hmm

# Processing parameters
KMER_MIN=27
KMER_MAX=127
KMER_STEP=10
MIN_COUNT=1
MIN_ORF_LENGTH=30
THREADS=32

# === MAIN PIPELINE ===
# Process each R1 FASTQ file in the current directory
for file in *_R1.fastq.gz; do
    # Extract base name by removing R1 suffix and extension
    name=${file%_R1.fastq.gz}
    
    # ASSEMBLY STAGE: De novo transcriptome assembly using MEGAHIT
    megahit \
        -1 $file \
        -2 ${name}_R2.fastq.gz \
        --k-min $KMER_MIN \
        --k-max $KMER_MAX \
        --k-step $KMER_STEP \
        --min-count $MIN_COUNT \
        -o ./megahit
    
    # ORF PREDICTION STAGE: Using TransDecoder to identify coding regions
    TransDecoder.LongOrfs \
        -t ./megahit/final.contigs.fa \
        -O ./transdecoder_virus \
        -m $MIN_ORF_LENGTH
    
    TransDecoder.Predict \
        -t ./megahit/final.contigs.fa \
        -O ./transdecoder_virus
    
    # Rename output files with sample-specific identifiers
    mv ./megahit/final.contigs.fa ${name}_final.contigs.fa
    mv ./transdecoder_virus/final.contigs.fa ${name}_contigs.fa
    mv ./transdecoder_virus/final.contigs.fa.transdecoder.bed ${name}_contigs.fa.transdecoder.bed
    mv ./transdecoder_virus/final.contigs.fa.transdecoder.cds ${name}_contigs.fa.transdecoder.cds
    mv ./transdecoder_virus/final.contigs.fa.transdecoder.gff3 ${name}_contigs.fa.transdecoder.gff3
    mv ./transdecoder_virus/final.contigs.fa.transdecoder.pep ${name}_contigs.fa.transdecoder.pep
    
    # VIRAL DETECTION STAGE: Multiple search strategies for viral proteins
    mkdir ${name}_virus_find
    
    # 1. MMseqs2 search: Fast sequence search against reference database
    mmseqs easy-search \
        ${name}_contigs.fa.transdecoder.pep \
        $MMSEQS_DB \
        ${name}_virus_find/${name}_mmseqs_note.m8 tmp

    # 2. DIAMOND BLASTP: Accelerated protein alignment
    diamond blastp \
        --db $DIAMOND_DB \
        --query ${name}_contigs.fa.transdecoder.pep \
        --out ${name}_virus_find/${name}_diamond_note.txt

    # 3. PSI-BLAST: Iterative search for distantly related sequences
    psiblast \
        -query ${name}_contigs.fa.transdecoder.pep \
        -out ${name}_virus_find/${name}_psiblast.out.txt \
        -db $BLAST_DB \
        -outfmt 6 \
        -num_threads $THREADS \
        -word_size 2 \
        -evalue 0.5

    # 4. HMMER search: Profile-based search against HMM database
    hmmsearch \
        --noali \
        --cpu $THREADS \
        --incE 0.5 \
        --domtblout ${name}_virus_find/${name}_hmmsearch_note.tsv \
        $HMM_DB \
        ${name}_contigs.fa.transdecoder.pep

    # RESULT INTEGRATION STAGE: Consolidate hits from multiple methods
    
    # Extract hit IDs from each tool's output
    csvtk cut -t -f 1 ${name}_virus_find/${name}_mmseqs_note.m8 > ${name}_virus_find/${name}_mmseqs_list.txt
    csvtk cut -t -f 1 ${name}_virus_find/${name}_diamond_note.txt > ${name}_virus_find/${name}_diamond_list.txt
    csvtk cut -t -f 1 ${name}_virus_find/${name}_psiblast.out.txt > ${name}_virus_find/${name}_psiblast_list.txt
    csvtk cut -t -f 1 ${name}_virus_find/${name}_hmmsearch_note.tsv > ${name}_virus_find/${name}_hmmsearch_list.txt

    # Combine all hit lists and remove duplicates
    cat ${name}_virus_find/${name}_mmseqs_list.txt \
        ${name}_virus_find/${name}_diamond_list.txt \
        ${name}_virus_find/${name}_psiblast_list.txt \
        ${name}_virus_find/${name}_hmmsearch_list.txt > \
        ${name}_virus_find/${name}_potential_virus_list.txt

    # Extract matching sequences from the protein database
    seqkit grep \
        --pattern-file ${name}_virus_find/${name}_potential_virus_list.txt \
        ${name}_contigs.fa.transdecoder.pep > \
        ${name}_virus_find/${name}_potential_virus.fasta

    # Remove duplicate sequences
    seqkit rmdup -s \
        ${name}_virus_find/${name}_potential_virus.fasta > \
        ${name}_virus_find/${name}_potential_virus_rmdup.fasta

    # Final validation with BLASTP
    blastp \
        -query ${name}_virus_find/${name}_potential_virus_rmdup.fasta \
        -out ${name}_virus_find/${name}_RdRp_list.txt \
        -db $BLAST_DB \
        -outfmt 6 \
        -num_threads $THREADS \
        -max_target_seqs 1
done