# TEProf3 (modified) for TE-derived transcript assembly

Modifications to the existing TEProf3 codebase have been made
to accommodate the Ensembl annotations and alternatively stranded
or mixed stranded protocols for joint assembly (https://github.com/Yonghao-Holden/TEProf3).

## Modified usage info

```
usage: teprof3 [-h] [--test] [--teprof2 TEPROF2] [--gdcjson GDCJSON] [--gtexjson GTEXJSON] [--split SPLIT] [--blastp BLASTP]
               [--blastpshort BLASTPSHORT] [--blastpdatabase BLASTPDATABASE] [--taxiddb TAXIDDB] [--blastjobnum BLASTJOBNUM] [--classify CLASSIFY]
               [--detectedprotein] [--parsenetmhcpan PARSENETMHCPAN] [--neoantigen] [--repeatmasker REPEATMASKER] [--geneannotation GENEANNOTATION]
               [--geneannotationfortranslation] [--hervannotation HERVANNOTATION] [--geneprotein] [-f MANIFEST] [-ki] [-s SAMPLENUMBER] [-g GUIDED]
               [-rs] [-v] [-am ASSEMBLEMODE] [-at ASSEMBLETHREAD] [-al ASSEMBLELENGTH] [-as ASSEMBLESAMPLENUMBER] [-ast ASSEMBLESTRAND]
               [-aj ASSEMBLEJUNCTIONREAD] [-pt PROCESSTPM] [-ps PROCESSSAMPLENUMBER] [-ptn PROCESSTRANSCRIPTNUMBER] [-fm FILTERMODE]
               [-fs FILTERSAMPLENUMBER] [-ft FILTERTHREAD] [-fa] [-fi FILTERINTRONRETENTION] [-fmo] [-fmot FILTERMONOEXONTPM]
               [-fdm FILTERDOWNSTREAMMATE] [-fr FILTERRATIO] [-fncf] [-fljt LONGREADSJTOLERANCE] [-flmt LONGREADMTOLERANCE] [-tt TACOTHREAD]
               [-tto TACOTOLERANCE] [-qm QUANMODE] [-qs QUANSAMPLENUMBER] [-qsc QUANSAMPLENUMBERCON] [-qnp] [-qng] [-ql QUANREADLENGTH]
               [-ti TRANSLATION] [-tm TRANSLATIONMODE] [-tl TRANSLATIONLENGTH] [-tg TRANSLATIONGENOME]

TEProf3 takes aligned files (bam files) and/or assembled data (gtf files) and provides you a list of TE-derived transcripts with expression across
samples.

options:
  -h, --help            show this help message and exit
  --test                testing mode
  --teprof2 TEPROF2     directly use teprof2 output for translation
  --gdcjson GDCJSON     json file downloaded from GDC for splice junction tsv files to generate sample_manifest.txt file
  --gtexjson GTEXJSON   json file downloaded from GTEx for splice junction tsv files to generate sample_manifest.txt file
  --split SPLIT         run blastp on sbatch cluster like htcf, put the fasta file here
  --blastp BLASTP       run blastp on sbatch cluster like htcf, put the fasta file here
  --blastpshort BLASTPSHORT
                        run blastp-short on sbatch cluster like htcf, put the fasta file here
  --blastpdatabase BLASTPDATABASE
                        database that blastp will run against (e.g., nr or gene)
  --taxiddb TAXIDDB     folder to taxid file
  --blastjobnum BLASTJOBNUM
                        number of blast jobs run simultaneously
  --classify CLASSIFY   path to teprof3_protein_information.tsv table from teprof3 --translation output
  --detectedprotein     process BLAST and BLAT result to get a list of detected TE-derived proteins
  --parsenetmhcpan PARSENETMHCPAN
                        parse output from netMHCpan
  --neoantigen          get list of neoantigen candidates
  --repeatmasker REPEATMASKER
                        repeatmasker file, only provide it if you are making a new reference
  --geneannotation GENEANNOTATION
                        gene annotation file, only provide it if you are making a new reference. GENCODE gtf is recommended
  --geneannotationfortranslation
                        gene annotation file for translation, only provide it if you are making a new reference. GENCODE gtf is recommended
  --hervannotation HERVANNOTATION
                        herv annotation file
  --geneprotein         generate protein sequence fasta file from gene annotation
  -f MANIFEST, --manifest MANIFEST
                        Dataset manifest file: sample name, short/long read, bam/gtf file (ended with bam or gtf). (tab-delimited). Optionally add
                        a 4th column for strand code (0,1,2).
  -ki, --keepintermediate
                        Keep intermediate files
  -s SAMPLENUMBER, --samplenumber SAMPLENUMBER
                        number of samples to be processed together at the same time (default 10)
  -g GUIDED, --guided GUIDED
                        run teprof3 in guided mode...
  -rs, --reset          reset the environment from a previous run
  -v, --version         Print version of teprof3
  -am ASSEMBLEMODE, --assemblemode ASSEMBLEMODE
                        how to run transcript de novo assembly: 0 = none (use provided GTF), 1 = short-read only, 2 = long-read only, 3 =
                        short+long in a single hybrid run (--mix), 4 = short & long assembled separately. (default 0)
  -at ASSEMBLETHREAD, --assemblethread ASSEMBLETHREAD
                        number of threads used for assembly (default 4)
  -al ASSEMBLELENGTH, --assemblelength ASSEMBLELENGTH
                        minimum transcript length for stringtie assembly (default 200)
  -as ASSEMBLESAMPLENUMBER, --assemblesamplenumber ASSEMBLESAMPLENUMBER
                        number of samples to be processed in parallel (default 10)
  -ast ASSEMBLESTRAND, --assemblestrand ASSEMBLESTRAND
                        default strandedness (0=unstranded, 1=--rf, 2=--fr). (default 0)
  -aj ASSEMBLEJUNCTIONREAD, --assemblejunctionread ASSEMBLEJUNCTIONREAD
                        minimum junction coverage (default 1)
  -pt PROCESSTPM, --processtpm PROCESSTPM
                        TPM cutoff to filter te-derived transcripts right after assembly (default 0.5)
  -ps PROCESSSAMPLENUMBER, --processsamplenumber PROCESSSAMPLENUMBER
                        number of samples to process at the same time (default 10)
  -ptn PROCESSTRANSCRIPTNUMBER, --processtranscriptnumber PROCESSTRANSCRIPTNUMBER
                        only include samples with > X TE-derived transcripts in sample. (default 100)
  -fm FILTERMODE, --filtermode FILTERMODE
                        how to filter TE-derived transcripts (1=only short read, 2=short+long) (default 1)
  -fs FILTERSAMPLENUMBER, --filtersamplenumber FILTERSAMPLENUMBER
                        number of samples to process at once (default 10)
  -ft FILTERTHREAD, --filterthread FILTERTHREAD
                        threads for filtering transcripts (default 10)
  -fa, --filterannotated
                        keep annotated transcripts. Default is to exclude them.
  -fi FILTERINTRONRETENTION, --filterintronretention FILTERINTRONRETENTION
                        cutoff for intron retention filtering. (default 3)
  -fmo, --filtermonoexon
                        exclude mono-exonic transcripts (default false).
  -fmot FILTERMONOEXONTPM, --filtermonoexontpm FILTERMONOEXONTPM
                        TPM cutoff for mono-exonic transcripts (default 1).
  -fdm FILTERDOWNSTREAMMATE, --filterdownstreammate FILTERDOWNSTREAMMATE
                        # reads capturing splicing junction to downstream exon (default 2).
  -fr FILTERRATIO, --filterratio FILTERRATIO
                        ratio used in chimeric mate filtering (default 0.5).
  -fncf, --filternochimericfilter
                        disable chimeric mate filter
  -fljt LONGREADSJTOLERANCE, --longreadsjtolerance LONGREADSJTOLERANCE
                        SJ tolerance for LR vs stringtie (default 3)
  -flmt LONGREADMTOLERANCE, --longreadmtolerance LONGREADMTOLERANCE
                        tolerance for mono-exonic check in LR data (default 50)
  -tt TACOTHREAD, --tacothread TACOTHREAD
                        number of threads for TACO (default 10)
  -tto TACOTOLERANCE, --tacotolerance TACOTOLERANCE
                        edge tolerance for TACO (default 3)
  -qm QUANMODE, --quanmode QUANMODE
                        1=short read, 2=SJ only (default 1)
  -qs QUANSAMPLENUMBER, --quansamplenumber QUANSAMPLENUMBER
                        samples in parallel for stringtie quant (default 10)
  -qsc QUANSAMPLENUMBERCON, --quansamplenumbercon QUANSAMPLENUMBERCON
                        samples in parallel for concatenating stringtie quant output (default 50)
  -qnp, --quannoprepde  skip prepDE.py for raw counts (default False)
  -qng, --quannogene    collect TPM for canonical genes (default True)
  -ql QUANREADLENGTH, --quanreadlength QUANREADLENGTH
                        read length of libraries (default 75)
  -ti TRANSLATION, --translation TRANSLATION
                        provide GTF or table for translation
  -tm TRANSLATIONMODE, --translationmode TRANSLATIONMODE
                        translation mode (1 or 2, default 1)
  -tl TRANSLATIONLENGTH, --translationlength TRANSLATIONLENGTH
                        minimum length of translated protein (default 20)
  -tg TRANSLATIONGENOME, --translationgenome TRANSLATIONGENOME
                        reference genome name for translation (default hg38)
```
