

module unload python
module load python/3.6.1
module load unicycler

#Porechop version 0.2.4
PC=/N/u/wrshoema/Carbonate/Porechop/porechop-runner.py
IN4=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled/run4
BARCODE_DIR4=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_pcbc/run4

$PC -i $IN4 --format fastq \
      -b $BARCODE_DIR4 \
      --barcode_threshold 75 \
      --barcode_diff 5 \
      --require_two_barcodes \
      --discard_middle

PMB12_nano=/N/dc2/projects/muri2/Task2/LTDE/data/PMB12/PMB12.fastq
PMB12_nano_clean=/N/dc2/projects/muri2/Task2/LTDE/data/PMB12/PMB12_clean.fastq

cp /N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_pcbc/run4/BC07.fastq $PMB12_nano

#NanoFilt version 2.3.0
NanoFilt=/N/u/wrshoema/Carbonate/.local/lib/python3.6/site-packages/nanofilt/NanoFilt.py
python $NanoFilt -q 10 -l 1000 --headcrop 100 $PMB12_nano > $PMB12_nano_clean

# get illumina reads from SRA using accession number SRR548025
/N/u/wrshoema/Carbonate/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump -I --split-files -O /N/dc2/projects/muri2/Task2/LTDE/data/PMB12 SRR548025
# per-base quality looks good, just go ahead and assemble

# run unicycler
#unicycler -s /N/dc2/projects/muri2/Task2/LTDE/data/PMB12/SRR548025_1.fastq \
#          -l $PMB12_nano_clean \
#          -o /N/dc2/projects/muri2/Task2/LTDE/data/PMB12/unicycler_assembly


# try spades
module load spades
spades.py --careful --cov-cutoff auto \
      --nanopore $PMB12_nano_clean \
      -s /N/dc2/projects/muri2/Task2/LTDE/data/PMB12/SRR548025_1.fastq \
      -o /N/dc2/projects/muri2/Task2/LTDE/data/PMB12/PMB12_spades_assembly
