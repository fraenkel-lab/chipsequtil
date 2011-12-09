Illumina pipeline script reference
==================================

The following is the output of the scripts provided by this package when invoked
on the command line with *-h*.

.. _top:

Scripts:
  - :ref:`chipseq_pipeline.py <chipseq_pipeline.py>`
  - :ref:`create_pipeline_script.py <create_pipeline_script.py>`
  - :ref:`extract_promoters.py <extract_promoters.py>`
  - :ref:`filter_bed_by_position_count.py <filter_bed_by_position_count.py>`
  - :ref:`filter_macs_peaks.py <filter_macs_peaks.py>`
  - :ref:`filter_gps_peaks.py <filter_gps_peaks.py>`
  - :ref:`filter_mapped_known_genes.py <filter_mapped_known_genes.py>`
  - :ref:`gerald_stats.py <gerald_stats.py>`
  - :ref:`gerald_to_bed.py <gerald_to_bed.py>`
  - :ref:`join_mapped_known_genes.py <join_mapped_known_genes.py>`
  - :ref:`map_intervals.py <map_intervals.py>`
  - :ref:`map_peaks_to_genes.py <map_peaks_to_genes.py>`
  - :ref:`map_peaks_to_known_genes.py <map_peaks_to_known_genes.py>`
  - :ref:`motif_scan.py <motif_scan.py>`
  - :ref:`nibFrag.py <nibFrag.py>`
  - :ref:`org_settings.py <org_settings.py>`
  - :ref:`peaks_to_fasta.py <peaks_to_fasta.py>`
  - :ref:`plot_pos_vs_neg_peaks.py <plot_pos_vs_neg_peaks.py>`
  - :ref:`plot_peak_loc_dist.py <plot_peak_loc_dist.py>`
  - :ref:`rejection_sample_fasta.py <rejection_sample_fasta.py>`
  - :ref:`sort_bed.py <sort_bed.py>`
  - :ref:`wait_for_jobid.py <wait_for_jobid.py>`
  - :ref:`wqsub.py <wqsub.py>`
  - :ref:`wqsub_drmaa.py <wqsub_drmaa.py>`


.. _chipseq_pipeline.py:

chipseq_pipeline.py::

   Usage: chipseq_pipeline.py [options] <organism> <experiment alignment filename> [<control alignment filename>]
   
   1st generation ChIPSeq analysis pipeline:
   
     - runs MACS to find peaks and sorts peaks by p-value
     - sorts peaks by pvalue and isolates top *n*
     - maps peaks to genes
     - extracts fasta files for gene peaks in experiments
     - constructs background sequences matching foreground distribution
     - runs THEME.py on input sequences w/ refinement
     - builds an infosite with stats from this analysis
   
   Control input file is optional.  *organism* argument is passed to the
   *org_settings.py* command to specify organism specific parameters, ensure
   that the following commands return valid paths:
   
   If running MACS:
    - org_settings.py <organism> genome_size
    - org_settings.py <organism> genome_dir
    - org_settings.py <organsim> refgene_anno_path
   
   If running THEME:
    - org_settings.py <organism> theme_hypotheses
    - org_settings.py <organism> theme_markov
   
   
   
   Options:
     -h, --help            show this help message and exit
     --auto                run all steps non-interactively (for batch mode, e.g.)
     --steplist=STEPLIST   with --auto, run specific steps
     --exp-name=EXP_NAME   name for the experiment/pipeline, used for convenience
                           [default: current directory name]
     --bed-args=BED_ARGS   double quote wrapped arguments for gerald_to_bed.py
                           [default: --stdout --chromo-strip=.fa]
     --macs-exec=MACS_EXEC
                           the executable to use for MACS, if not an absolute
                           path it needs to be on your shell environment path
                           [default: macs14]
     --macs-args=MACS_ARGS
                           double quote wrapped arguments for macs, only changing
                           --mfold, --tsize, --bw, and --pvalue recommended
                           [default: --pvalue=1e-5]
     --map-args=MAP_ARGS   double quote wrapped arguments for mapping peaks to
                           genes [default: --tss --upstream-window=10000
                           --downstream-window=10000]
     --filter-peaks-args=FILTER_PEAKS_ARGS
                           double quote wrapped arguments for
                           filter_macs_peaks.py [default: --sort-by=pvalue
                           --top=1000 -f 'tags>20']
     --filter-neg-peaks-args=FILTER_NEG_PEAKS_ARGS
                           double quote wrapped arguments for
                           filter_macs_peaks.py applied to negative peaks
                           [default: -f 'tags>20']
     --peaks-to-fa-args=PEAKS_TO_FA_ARGS
                           double quote wrapped arguments for peaks_to_fasta.py
                           [default: --fixed-peak-width=200]
     --bg-exec=BG_EXEC     the executable to use for generating background
                           sequences for THEME, if not an absolute path it needs
                           to be on your shell environment path [default:
                           rejection_sample_fasta.py]
     --bg-args=BG_ARGS     double quote wrapped arguments for background sequence
                           generation utility [default: --num-seq=2.1x]
     --theme-args=THEME_ARGS
                           double quote wrapped arguments for THEME.py [default:
                           --beta=0.7 --cv=5 --trials=25]
     --motif-pval-cutoff=MOTIF_PVAL
                           the p-value cutoff for sending non-refined enrichmed
                           motifs to THEME for refinement
     --parallelize         parallelize portions of the pipeline using qsub, only
                           works from SGE execution hosts
     --ucsc                perform tasks for automated integration with UCSC
                           genome browser [default:False]
     --build-infosite-args=INFOSITE_ARGS
                           arguments to pass to build_chipseq_infosite.py
                           [default: None]
   
     UCSC Integration Options (with --ucsc):
       --stage-dir=STAGE_DIR
                           root directory where UCSC integration files should be
                           made available [default: ./]
       --stage-url=STAGE_URL
                           URL where UCSC integration files will be made
                           available over the web [default: http://localhost/]
   
   Note: it is advised to leave the --*-args arguments unchanged
   unless you really know what you're doing.
      

:ref:`top <top>`

.. _create_pipeline_script.py:

create_pipeline_script.py::

   This is an interactive script that creates an executable script to use
   for ChIPSeq analyses. When prompted for experiment and control files,
   tab completion is available a la bash or tcsh shells. Press Ctrl-C at
   any time to quit.
   Usage: create_pipeline_script.py
   
   Script for creating a custom run script for ChIPSeq/DNAse hypersensitivity
   experiments.  User is asked for paths and settings required for ChIPSeq
   analysis using the *chipseq_pipeline.py* utility and produces an executable
   run script with helpful information on how to run it.  Also creates a JSON
   formatted file containing all the parameters for this pipeline run.
   
   Options:
     -h, --help  show this help message and exit
   
   Note: this script only works in Unix-style environments
      
   ================= ChIPSeq Experiment Pipeline Script Generator =================
   

:ref:`top <top>`

.. _extract_promoters.py:

extract_promoters.py::

   Usage: extract_promoters.py [options] <organism>
   
   Extract the promoter sequences in FASTA format from all genes
   or a list of genes specified in an input file.  Gene annotation is RefGene
   corresponding to the organism passed in, paths returned by:
   
   $> org_settings.py <organism> refgene_anno_path
   $> org_settings.py <organism> genome_dir
   
   must be valid.
   
   Options:
     -h, --help            show this help message and exit
     -u UPSTREAM, --upstream=UPSTREAM
                           upstream window from TSS to extract [default: 3000]
     -d DOWNSTREAM, --downstream=DOWNSTREAM
                           downstream window from TSS to extract [default: 1000]
     -l GENE_LIST, --gene-list=GENE_LIST
                           file containing a list of gene identifiers to extract,
                           one per line [default: none]
     -t GENE_TYPE, --gene-type=GENE_TYPE
                           type of gene identifier in gene list, choose from
                           ['symbol', 'refgene'] [default: symbol]
     -o OUTPUT, --output=OUTPUT
                           file to write fasta records to [default: stdout]
      

:ref:`top <top>`

.. _filter_bed_by_position_count.py:

filter_bed_by_position_count.py::

   Usage: filter_bed_by_position_count.py [options] <bed file>
   
   Analyze BED file and filter out alignments above some threshold that align to
   a single genomic position.
   
   Options:
     -h, --help            show this help message and exit
     -n MAX_COUNT, --max-count=MAX_COUNT
                           max tag count at a given position, filter above
                           [default: 5]
     --output=OUTPUT       write output to file
   
   Note: only works if BED file is sorted!
      

:ref:`top <top>`

.. _filter_macs_peaks.py:

filter_macs_peaks.py::

   Usage: filter_macs_peaks.py [options] <MACS peak file>
   
   Filter MACS peaks by supplied criteria.  Available filter features are:
   
   length
   tags
   pvalue
   fold_enrichment
   fdr
   
   Filters are provided as expressions using the [-f |--filter] option, e.g. the
   command
   
   filter_macs_peaks.py -f "tags>100" --filter="pvalue<=1e-9"
   --filter="100<length<=200" <MACS peak file>
   
   finds only peaks with more than 100 tags, a pvalue of less than 1e9, and a
   length between 100, exclusive, and 200, inclusive.  Any number of filters may
   be provided, and only peaks that match *all* filters pass.  User is warned if
   filters result in zero results.  Only inequality operators are valid.
   Invoking with no filter arguments returns all peaks.  To sort, use the --sort-
   by option, e.g.
   
   filter_macs_peaks.py -f "pvalue<=1e-9" --sort-by=pvalue <MACS peak file>
   
   sorts peaks with a pvalue smaller than 1e-9 with the smallest pvalue peaks.
   All fields are sorted ascending by default.  Output is prepended with comments
   describing what the file contains, i.e. which filters are applied, how many
   records there are, etc.
   
   Note: MACS -10*log10(pvalue) values are converted to normal pvalues
   
   
   Options:
     -h, --help            show this help message and exit
     -f FILTERS, --filter=FILTERS
                           add filter expression
     --sort-by=SORT_BY     comma delimited list of features to sort by, filtered
                           peaks are not sorted by default, if provided peaks are
                           sorted ascending by default
     --sort-dir=SORT_DIR   direction to sort [default: ASCEND]
     --top=TOP             accepts an integer, output at most this many peaks
                           [default: all]
     --output=OUTPUT       filename to output filtered peaks to [default: stdout]
     --encode-filters      write out records to a file <MACS peaks
                           file>_<filters>.xls (incompatible with --output
                           option)
     --summary             only print out summary information for the filter
     --no-header           do not print out header or metadata info
     --shuffle             shuffle order of filtered records, useful for
                           selecting random peaks
     --print-encoded-fn    print out the filename that would be created by
                           --encode-filters
      

:ref:`top <top>`

.. _filter_gps_peaks.py:

filter_gps_peaks.py::

   Usage: filter_gps_peaks.py [options] <GPS peak file>
   
   Filter GPS peaks by supplied criteria.  Available filter features are:
   
   IP
   Control
   Fold
   qvalue
   pvalue
   IPvsEMP
   IPvsCTR
   
   Filters are provided as expressions using the [-f |--filter] option, e.g. the
   command
   
   filter_gps_peaks.py -f "IP>100" --filter="pvalue<=1e-9" <GPS peak file>
   
   finds only peaks with more than 100 tags and a pvalue of less than 1e9.  Any
   number of filters may be provided, and only peaks that match *all* filters
   pass. User is warned if filters result in zero results.  Only inequality
   operators are valid.  Invoking with no filter arguments returns all peaks.  To
   sort, use the --sort-by option, e.g.
   
   filter_gps_peaks.py -f "pvalue<=1e-9" --sort-by=pvalue <GPS peak file>
   
   sorts peaks with a pvalue smaller than 1e-9 with the smallest pvalue peaks.
   All fields are sorted ascending by default.  Output is prepended with comments
   describing what the file contains, i.e. which filters are applied, how many
   records there are, etc.
   
   Note: GPS P_-log10 and Q_-log10 values are converted to normal pvalues and
   qvalues
   
   
   Options:
     -h, --help            show this help message and exit
     -f FILTERS, --filter=FILTERS
                           add filter expression
     --sort-by=SORT_BY     comma delimited list of features to sort by, filtered
                           peaks are not sorted by default, if provided peaks are
                           sorted ascending by default
     --sort-dir=SORT_DIR   direction to sort [default: ASCEND]
     --top=TOP             accepts an integer, output at most this many peaks
                           [default: all]
     --output=OUTPUT       filename to output filtered peaks to [default: stdout]
     --encode-filters      write out records to a file <GPS peaks
                           file>_<filters>.xls (incompatible with --output
                           option)
     --summary             only print out summary information for the filter
     --no-header           do not print out header or metadata info
     --shuffle             shuffle order of filtered records, useful for
                           selecting random peaks
     --print-encoded-fn    print out the filename that would be created by
                           --encode-filters
      

:ref:`top <top>`

.. _filter_mapped_known_genes.py:

filter_mapped_known_genes.py::

   Usage: filter_mapped_known_genes.py [options] <mapped known genes file>
   
   Filter columns and rows from *join_mapped_known_genes.py* output which was
   invoked with *--binary-plus* and *--field-types* flags.  Specify full column
   names for either binding or expression data with the *--bind-cols* and
   *--affy-cols* arguments, respectively. The special fieldname *MAPPED* from
   *join_mapped_known_genes.py* is used to determine whether a file contains a
   mapping for each gene.  To filter genes by their associated binding or
   expression data, specify *--bind-filter* or *--affy-filter* as follows:
   
     - *any* - report gene if at least one input file maps to the gene
     - *all* - report gene if every input file maps to the gene
     - *absent* - report gene if no input file maps to the gene
     - *none* - do not filter genes at all (default)
   
   Results of binding and expression filters are 'and'ed together, e.g.:
   
   --bind-filter=all --affy-filter=absent
   
   returns only genes for which all binding files and none of the expression
   files map.
   
   
   Options:
     -h, --help            show this help message and exit
     --bind-cols=BIND_COLS
                           comma delimited list of binding data column names to
                           include, [default: all]
     --affy-cols=AFFY_COLS
                           comma delimited list of expression data column names
                           to include, [default: all]
     --bind-filter=BIND_FILT
                           gene set to include based on binding data [default:
                           none]
     --affy-filter=AFFY_FILT
                           gene set to include based on expression data [default:
                           none]
     --output=OUTPUT       write output to file
   
   Note: when specifying column names, be sure to escape characters like
   (,),&,*,etc... that shells interpret with a \, e.g. --bind-
   cols=-10\*log10\(pvalue\)
      

:ref:`top <top>`

.. _gerald_stats.py:

gerald_stats.py::

   Usage: gerald_stats.py [options] <filename> [<filename>...]
   
   Outputs various stats about the GERALD formatted file(s) input. If multiple
   files are provided statistics are aggregated according to the specified output
   format.  Output formats available via --format=X :
   
     # *python* - print an eval()'able python dictionary w/ counts
     # *rst* - print statistics in a reStructured text table (default)
     # *tab* - print statistics in a tab delimited form w/ header names
   
   Except for *python* format, each input file has its own output line.  *python*
   summarizes all alignments.
   
   
   Options:
     -h, --help       show this help message and exit
     --output=OUTPUT  write output to file [default: stdout]
     --format=FORMAT  format to print out stats [default: rst]
      

:ref:`top <top>`

.. _gerald_to_bed.py:

gerald_to_bed.py::

   Usage: gerald_to_bed.py [options] <GERALD file> [<GERALD file>...]
   
   Convert the GERALD alignment formatted files into BED format.  Input file
   named <path>/<filename>.<ext> is translated into <path>/<filename>.bed unless
   --output or --stdout is specified, in which case formatted lines are written
   to file or standard output, respectively.  If multiple input files are
   supplied with the --output or --stdout option all formatted lines are
   concatenated together. Formatting only occurs for GERALD input lines that have
   a valid Match Position field (i.e. successfully aligned somewhere).
   
   Options:
     -h, --help            show this help message and exit
     --output=OUTPUT       write all records to file
     --stdout              write out all formatted lines to stdout
     --min-fields          only format the first three fields
     --pass-only           only format lines with Y in the Pass Filtering field
     --chromo-strip=CHROMO_STRIP
                           pattern to remove from chromo field in BED output
                           (e.g. --chromo-strip=.fa to remve .fa from chrX.fa)
                           [default: .fa]
      

:ref:`top <top>`

.. _join_mapped_known_genes.py:

join_mapped_known_genes.py::

   Usage: join_mapped_known_genes.py -b <mapped DNA binding file>|-a <mapped microarray file> [-b <mapped DNA binding file> ...] [-a <mapped microarray file> ...]
   
   Join all files on the first column, concatenating records with matching
   entries onto one line per entry.  Understands DNA binding data as mapped with
   *map_peaks_to_known_genes.py* utility microarray data as mapped by
   *probeset_to_known_genes.py* utility, passed to program using *-b* and *-a*
   options respectively.  If a file contains more than one mapping to a gene
   additional columns are added. At least one file of either type is required.
   Field names are written as <filename>.<original field name>.<map number>
   
   Options:
     -h, --help            show this help message and exit
     -a AFFY_FILE, --affy-file=AFFY_FILE
                           add a mapped microarray file
     -b BIND_FILE, --bind-file=BIND_FILE
                           add a mapped DNA binding file (e.g. MACS, BED)
     -m MACS_FILE, --macs-file=MACS_FILE
                           DEPRECATED: use -b instead, add a mapped default MACS
                           formatted peaks (*.xls) file
     --output=OUTPUT       file to output joined records to [default: stdout]
     --first-only          only output the first mapping to a gene from each file
     --binary              output only one column per file with a 0 or 1 to
                           indicate whether a mapping exists in that file
     --binary-plus         output one column per file with a 0 or 1 to indicate
                           whether a mapping exists in that file in addition to
                           all other columns
     --field-types         prepend BIND or AFFY to the beginning of all
                           appropriate columns
   
   Note: microarray files should have been created by bioconductor, and all files
   should have a row of fieldnames as the first line
      

:ref:`top <top>`

.. _map_intervals.py:

map_intervals.py::

   Usage: map_intervals.py [options] <from> <to>
   
   Find records in <to> interval file that map to records in <from> interval
   file.  Files should be tab delimited and are expected to have a chromosome
   column, a start column, and an end column.  The indices of these columns can
   be specified on the command line but by default are the first three columns,
   respectively.  Prints out to stdout by default one new line separated row per
   row in <from> with a line from <to> where there is a mapping. If no mapping is
   found (e.g. when specifying a maximum margin to search within) the word None
   is printed.  By default only prints nearest record, with ties settled by
   smallest line number in <to>.
   
   Options:
     -h, --help            show this help message and exit
     -w WINDOW, --window=WINDOW
                           window as <int upstream> <int downstream> to search
                           for intervals [default: (1000000000.0, 1000000000.0)]
     -f FROM_IND, --from=FROM_IND
                           coordinates of chromosome, start, stop in <from> file
     -i, --skip-from-header
                           <from> has a header that should be skipped
     -t TO_IND, --to=TO_IND
                           coordinates of chromosome, start, stop in <to> file
     -j, --skip-to-header  <to> has a header that should be skipped
      

:ref:`top <top>`

.. _map_peaks_to_genes.py:

map_peaks_to_genes.py::

   Usage: map_peaks_to_genes.py [options] <refGene file> <peaks file>
   
    Map the peaks in <peaks file> to genes in <refGene file>.  <refGene file> is
   format is as specified in
   http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.sql. <peaks
   file> format is as produced by MACS.
   
   Options:
     -h, --help            show this help message and exit
     --upstream-window=UPST_WIN
                           window width in base pairs to consider promoter region
                           [default: 5500]
     --downstream-window=DNST_WIN
                           window width in base pairs to consider downstream
                           region [default: 2500]
     --map-output=PEAK_OUTPUT
                           filename to output mapped peaks in BED format to
                           [default: stdout]
     --stats-output=STATS_OUTPUT
                           filename to output summary stats in conversion
                           [default: stderr]
     --peaks-format=PEAKS_FMT
                           format of peaks input file [default: MACS]
      

:ref:`top <top>`

.. _map_peaks_to_known_genes.py:

map_peaks_to_known_genes.py::

   Usage: map_peaks_to_known_genes.py [options] <knownGene file> <knownGene xRef file> <peaks file>
   
   
   Map the peaks in <peaks file> to genes in <knownGene file>.  <knownGene file>
   isformat is as specified in
   http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.sql.<peaks
   file> format is as produced by MACS.  If *auto* is chosen (default) file
   extension is examined for *.xls* for default MACS format or *.bed* for BED
   format.  If the --detailoption is provided, the following extra fields are
   appended to each row:
   
   peak loc, dist from feature, score, map type, map subtype
   
   
   Options:
     -h, --help            show this help message and exit
     --upstream-window=UPST_WIN
                           window width in base pairs to consider promoter region
                           [default: 5500]
     --downstream-window=DNST_WIN
                           window width in base pairs to consider downstream
                           region [default: 2500]
     --tss                 calculate downstream window from transcription start
                           site instead of transcription end site
     --map-output=PEAK_OUTPUT
                           filename to output mapped peaks to [default: stdout]
     --stats-output=STATS_OUTPUT
                           filename to output summary stats in conversion
                           [default: stderr]
     --peaks-format=PEAKS_FMT
                           format of peaks input file [default: auto]
     --detail              add extra fields to output, see description
     --intergenic          write intergenic peaks to the gene file as well with
                           None as gene ID
      

:ref:`top <top>`

.. _motif_scan.py:

motif_scan.py::

   Usage: motif_scan.py [options] <org> <peaks fn> <TAMO motif fn>
   
   Do some motif scanning stuffs
   
   Options:
     -h, --help            show this help message and exit
     -n TOP_N, --top-n=TOP_N
                           use top n peaks by pvalue for sequence scanning
                           [default: all]
     -i MOTIF_IND, --motif-indices=MOTIF_IND
                           which indices from <TAMO motif fn> to use [default:
                           all]
     -d DIR, --dir=DIR     write all results into this directory
     --fixed-peak-width=FIXED_W
                           use only a fixed peak window around the summit instead
                           of whole peak
      

:ref:`top <top>`

.. _nibFrag.py:

nibFrag.py::

   Usage: nibFrag.py [options] file.nib start end strand [outfile]
     -- or --
   nibFrag.py [options] --batch file.nib batchfile [batchfile ...]
   
   A python implementation of Jim Kent's nibFrag utility that allows outputting
   to stdout.  Otherwise the functionality is identical for the non-batch usage.
   Batch mode accepts one or more files containing sets of coordinates to extract
   from the nib file.  Only BED formatting is accepted at the moment. All
   sequences are concatenated together in FASTA format.  To retrieve the entire
   sequence, use END as the end argument.
   
   Options:
     -h, --help            show this help message and exit
     --no-header           only output sequence (no fasta header)
     --wrap-width=WRAP_WIDTH
                           wrap output sequence at this number of bases, 0
                           indicates no wrap (sequence ends up on single line)
                           [default: 50]
     --batch               run in batch mode, interpret arguments after nib file
                           as queries
     --batch-format=BATCH_FORMAT
                           format to interpret batch files [default: BED]
   
     Original nibFrag options:
       --masked            use lower case characters for bases meant to be masked
                           out
       --hardMasked        use upper case for non masked-out and 'N' characters
                           for masked-out bases
       --upper             use upper case characters for all bases
       --name=NAME         Use given name after '>' in output sequence
       --dbHeader=DBHEADER
                           Add full database info to the header, with or without
                           -name option
       --tbaHeader=TBAHEADER
                           Format header for compatibility with tba, takes
                           database name as argument
   
   Note: When specifying --name optionin batch mode, also specify --dbHeader to
   ensure unique FASTA headers.
      

:ref:`top <top>`

.. _org_settings.py:

org_settings.py::

   Usage: org_settings.py [options] [<org key> [<org setting>]]
   
   Tool for retrieving sets of organism-specific settings and paths. Original
   paths are set at install time, and can be overridden in the file ~/.org
   settings.cfg. Allows output of settings in a variety of shell environment
   syntaxes.  The tool attempts to guess which shell environment is being used by
   examining the SHELL environment variable unless explicitly set.  When run
   without an argument, returns a listing of all settings available.
   
   Options:
     -h, --help            show this help message and exit
     -s SYNTAX, --syntax=SYNTAX
                           syntax flavor                   of output to produce
                           [default: %auto]
     -l, --list            print                   all available settings for
                           human consumption
      

:ref:`top <top>`

.. _peaks_to_fasta.py:

peaks_to_fasta.py::

   Usage: peaks_to_fasta.py [options] <organism> <peak file> [<peak file> ...]
   
   Extract sequences for peaks in provided peak file(s).  Can interpret MACS or
   BED output, determined automatically by .xls or .bed extensions respectively
   (force explicit format with --peak-format option).  Outputs fasta sequences
   for the peaks in all files extracted from the reference genome specified by
   the output of *org_settings.py <organism> genome_dir* to stdout by
   default.Chromosome names in peak files must match nib filenames without
   extension (e.g. peak line: chr1 0  100 searches *genome_dir*/chr1.nib).  Fasta
   records have the following format:
   
   ><chromosome>:<start>-<end>;fn=<name of file>:<line number>;db_fn=<db
   filename>;fmt=<format>;<source alignment info>
   <sequence...>
   
   <db filename> is the filename where the sequence was extracted, <format> is
   the format of the input file (MACS or BED), and <source alignment info>
   contains all the fields from the originating alignment according to the source
   format.
   
   Options:
     -h, --help            show this help message and exit
     --min-header          only store <chromosome>:<start>-<end> in header
     --peak-format=PEAK_FORMAT
                           peak file format, 'auto' determines format by
                           extension, choices: MACS, BED, auto [default: auto]
     --output=OUTPUT       filename to output fasta records to [default: stdout]
     --fixed-peak-width=FIXED_PEAK_WIDTH
                           return a fixed number of bases flanking peak summit
                           (*summit* field in MACS, (end-start)/2 in BED),
                           ignoring start/stop coords [default: None]
     --wrap-width=WRAP_WIDTH
                           wrap fasta sequences to specified width. -1 indicates
                           no wrap [default: 70]
      

:ref:`top <top>`

.. _plot_pos_vs_neg_peaks.py:

plot_pos_vs_neg_peaks.py::

   Usage: plot_pos_vs_neg_peaks.py [options] <pos peaks fn> <neg peaks fn>
   
   Options:
     -h, --help            show this help message and exit
     -o OUT_FN, --output=OUT_FN
                           filename of output image
      

:ref:`top <top>`

.. _plot_peak_loc_dist.py:

plot_peak_loc_dist.py::

   Usage: plot_peak_loc_dist.py [options] <peaks fn> <gene list fn>
   
   Produce a pie chart of the locations of peaks in different bins (promoter,
   gene, exon, intron, etc.) and, optionally, save the different records to their
   own files for subsequent analysis.  Also produce a histogram of distance from
   feature values in mapping file. Peaks file is expected to be as output by
   MACS, or alternately as a BED file but then the -b plot is not available.
   Gene list file is expected to be in the format as output by
   peaks_to_known_genes.py script.
   
   Options:
     -h, --help            show this help message and exit
     -b BAR_FN, --bar-fn=BAR_FN
                           filename for pvalue stacked bar chart
     -g GENE_PIE_FN, --gene-pie-fn=GENE_PIE_FN
                           filename for pie chart image
     -p PEAK_PIE_FN, --peak-pie-fn=PEAK_PIE_FN
                           filename for pie chart image
     -f DIST_FN, --dist-fn=DIST_FN
                           filename for distance from feature image
     -s, --save            write out files containing peaks for each category
     -d OUT_DIR, --output-dir=OUT_DIR
                           output files created by --save option to this
                           directory
     --no-plot             dont show (but save) the figure produced
     --peaks-format=PEAK_FMT
                           format of peaks file, either MACS or BED [default:
                           MACS]
      

:ref:`top <top>`

.. _rejection_sample_fasta.py:

rejection_sample_fasta.py::

   Usage: rejection_sample_fasta.py [options] <organism> <fasta file> [<fasta file> ... ]
   
   Use rejection sampling to generate a set of background/random
   sequences matching the distance to nearest transcription start site, sequence
   length, and GC content distributions of the input fasta file(s).  Generated
   sequences are genomic sequences sampled based on these distributions. All
   sequences
   from all files are used to generate the background sequences. The following
   command must output a path to a nib genomic sequence directory and refGene
   annotation, respectively :
   
   $> org_settings.py <organism> genome_dir
   $> org_settings.py <organism> refgene_anno_path
   
   Utility prints out generated fasta records to stdout by default.  Input
   sequences
   from chr20 are mapped to chrX, chr21 are mapped to chrY, and sequences from
   chrM
   are not used.
   
   
   Options:
     -h, --help            show this help message and exit
     -n NUM_SEQS, --num-seqs=NUM_SEQS
                           number of sequences to generate, either absolute
                           number or factor of # input sequences, e.g. 2.5x for
                           2.5 times the # of input sequences [default: 1x]
     --output=OUTPUT       file to output fasta records to [default: stdout]
     --bed                 also produce a BED formatted file representing sampled
                           sequences
     --bed-output=BED_OUTPUT
                           with --bed, file to output BED records to [default:
                           output.bed]
     -v, --verbose         print out debug information
      

:ref:`top <top>`

.. _sort_bed.py:

sort_bed.py::

   Usage: sort_bed.py [options] <BED file> [<BED file> <BED file>...]
   
   Sort the BED formatted files first by chromosome (field 1) and then by start
   coordinate (field 2).  Lines from all files submitted are concatenated and
   sorted in the final output.
   
   Options:
     -h, --help       show this help message and exit
     --output=OUTPUT  filename to write the sorted BED lines [default: stdout]
      

:ref:`top <top>`

.. _wait_for_jobid.py:

wait_for_jobid.py::

   Usage: wait_for_jobid.py [options] <job id> [<job id>...]
   
   Poll qstat and wait until all <job id>s are finished
   
   Options:
     -h, --help  show this help message and exit
      

:ref:`top <top>`

.. _wqsub.py:

wqsub.py::

   Usage: [wqsub.py] [options] command
   
   Wrap the specified command into a qsub script and submit it for execution.
   Script captures both stdout and stderr to the current directory. By default,
   all of the user's environment variables are put into the script (compatible
   with SGE only ATM).
   
   Options:
     -h, --help            show this help message and exit
     --wqsub-name=WQSUB_NAME
                           job name to submit as <--wqsub-name>_<first non-
                           whitespace chars in command> [default: wqsub]
     --wqsub-ext=WQSUB_EXT
                           file extension to use for stdout files
     --wqsub-keep-script   do not delete qsub script generated after job
                           submission
     --wqsub-no-env        do not include any local environment variables in the
                           script
     --wqsub-no-submit     create script but do not submit job (useful for
                           generating scripts)
     --wqsub-drm=DRM       the DRM to generate scripts for [default: SGE]
     --wqsub-drm-arg=DRM_ARGS
                           arguments to pass as parameters in the job script
                           specific to the DRM, use multiple option flags to
                           specify multiple parameters
     --wqsub-wait          poll the DRM and do not return control until job is
                           finished (only works for TORQUE)
   
   Note: this script only works in Unix-style environments.
      

:ref:`top <top>`

.. _wqsub_drmaa.py:

wqsub_drmaa.py::

      Traceback (most recent call last):
     File "../scripts/wqsub_drmaa.py", line 9, in <module>
       import drmaa
   ImportError: No module named drmaa
   

:ref:`top <top>`



