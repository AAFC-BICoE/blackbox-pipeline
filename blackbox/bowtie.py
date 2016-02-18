#!/usr/bin/env python
from Bio.Application import _Option, AbstractCommandline, _Switch
import re

__author__ = 'mike knowles'
__doc__ = 'Wrapper for bowtie2'


class MakeBlastDB(AbstractCommandline):
    """Base makeblastdb wrapper"""

    def __init__(self, cmd='bowtie2', **kwargs):
        assert cmd is not None
        extra_parameters = [
            # Core:
            _Switch(["-h", "h"],
                    "Print USAGE and DESCRIPTION;  ignore other arguments."),
            _Switch(["-help", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description; "
                    "ignore other arguments."),
            _Switch(["-version", "version"],
                    "Print version number;  ignore other arguments."),

            # Input Options
            _Switch(["-q", "fastq"],
                    "Reads (specified with <m1>, <m2>, <s>) are FASTQ files. FASTQ files usually have "
                    "at. See also: --solexa-quals and --int-quals."),
            _Switch(["--qseq", "qseq"],
                    "Reads (specified with <m1>, <m2>, <s>) are QSEQ files. QSEQ files usually end in s."),
            _Switch(["-f", "fasta"],
                    "Reads (specified with <m1>, <m2>, <s>) are FASTA files. FASTA files usually have "
                    "ore-quals is also set."),
            _Switch(["-r", "unformated"],
                    "Reads (specified with <m1>, <m2>, <s>) are unformated files. With one input sequence per "
                    "if --ignore-quals is also set."),
            _Switch(["-c", "csv"],
                    "The read sequences are given on command line. I.e. <m1>, <m2> and <singles> are CSV files of "
                    "reads rather than lists of or qualities, so -c also implies --ignore-quals."),
            _Switch(["--phred33", "phred33"],
                    "Input qualities are ASCII chars equal to the Phred quality plus 33. This is also called "
                    "the Phred+33 encoding, which is used by the very latest Illumina pipelines"),

            _Switch(["--phred64", "phred64"],
                    "Input qualities are ASCII chars equal to the Phred quality plus 64. This is also called "
                    "the Phred+64 encoding"),

            _Switch(["--solexa-quals", "solexa_quals"],
                    "Convert input qualities from Solexa (which can be negative) to Phred (which can't). This "
                    "scheme was used in older Illumina GA Pipeline versions (prior to 1.3). Default: off"),

            _Switch(["--int-quals", "int_quals"],
                    "Quality values are represented in the read input file as space-separated ASCII integers, "
                    "e.g., 40 40 30 40..., rather than ASCII characters, e.g., II?I.... Integers are treated as "
                    "being on the Phred quality scale unless --solexa-quals is also specified. Default: off"),

            # Preset options in --end-to-end mode
            _Switch(["--very-fast", "very_fast"],
                    "Same as: -D 5 -R 1 -N 0 -L 22 -i S,0,2.50"),

            _Switch(["--fast", "fast"],
                    "Same as: -D 10 -R 2 -N 0 -L 22 -i S,0,2.50"),

            _Switch(["--sensitive", "sensitive"],
                    "Same as: -D 15 -R 2 -L 22 -i S,1,1.15 (default in --end-to-end mode)"),

            _Switch(["--very-sensitive", "very_sensitive"],
                    "Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50"),

            # Preset options in --local mode
            _Switch(["--very-fast-local", "very_fast_local"],
                    "Same as: -D 5 -R 1 -N 0 -L 25 -i S,1,2.00"),

            _Switch(["--fast-local", "fast_local"],
                    "Same as: -D 10 -R 2 -N 0 -L 22 -i S,1,1.75"),

            _Switch(["--sensitive-local", "sensitive_local"],
                    "Same as: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default in --local mode)"),

            _Switch(["--very-sensitive-local", "very_sensitive_local"],
                    "Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50"),
            # Input configuration options
            _Option(["--skip", "skip"],
                    "Skip (i.e. do not align) the first <int> reads or "
                    "pairs in the input",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--qupto", "qupto"],
                    "Align the first <int> reads or read pairs from the"
                    " input (after the -s/--skip reads or pairs have been skipped), then stop. Default: no limit",
                    checker_function=lambda value: type(value) is int,
                    equate=False),

            _Option(["--trim5", "trim5"],
                    "Trim <int> bases from 5' (left) end of each read before alignment (default: 0)",
                    checker_function=lambda value: type(value) is int,
                    equate=False),

            _Option(["--trim3", "trim3"],
                    "Trim <int> bases from 3' (right) end of each read before alignment (default: 0)",
                    checker_function=lambda value: type(value) is int,
                    equate=False),

            # Alignment options
            _Option(["-N", "num_mismatches"],
                    "Sets the number of mismatches to allowed in a seed alignment during multiseed "
                    "alignment. Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower) "
                    "but increases sensitivity. Default: 0",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["-L", "seed_length"],
                    "Sets the length of the seed substrings to align during multiseed alignment. "
                    "Smaller values make alignment slower but more senstive. Default: the --sensitive preset is used "
                    "by default, which sets -L to 20 both in --end-to-end mode and in --local mode",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["-i", "i_func"],
                    "Sets a function governing the interval between seed substrings to use during multiseed alignment. "
                    "For instance, if the read has 30 characters, and seed length is 10, and the seed interval is 6, "
                    "the seeds extracted will be: Since it's best to use longer intervals for longer reads, this "
                    "parameter sets the interval as a function of the read length, rather than a single one-size-fits-"
                    "all number. For instance, specifying -i S,1,2.5 sets the interval "
                    "function f to f(x) = 1 + 2.5 * sqrt(x), where x is the read length. "
                    "See also: setting function options. If the function returns a result less than 1, it is rounded up"
                    " to 1. Default: the --sensitive preset is used by default, which sets -i to S,1,1.15 "
                    "in --end-to-end mode to -i S,1,0.75 in --local mode.",
                    checker_function=lambda value: re.match('^[CLSG],[-\d\.],[-\d\.]', value) is not None,
                    equate=False),
            _Option(["--n-ceil", "n_ceil"],
                    "Sets a function governing the maximum number of ambiguous characters (usually Ns and/or .s) "
                    "allowed in a read as a function of read length. For instance, specifying -L,0,0.15 sets the "
                    "N-ceiling function f to f(x) = 0 + 0.15 * x, where x is the read length. See also: setting "
                    "function options. Reads exceeding this ceiling are filtered out. Default: L,0,0.15.",
                    checker_function=lambda value: re.match('^[CLSG],[-\d\.],[-\d\.]', value) is not None,
                    equate=False),
            _Option(["--gbar", "-gbar"],
                    "Disallow gaps within <int> positions of the beginning or end of the read. Default: 4.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--dpad", "-dpad"],
                    "Pads dynamic programming problems by <int> columns on either side to allow gaps. Default: 15.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["--ignore-quals", "ignore-quals"],
                    "When calculating a mismatch penalty, always consider the quality value at the mismatched position "
                    "to be the highest possible, regardless of the actual value. I.e. input is treated as though all "
                    "quality values are high. This is also the default behavior when the input doesn't specify quality "
                    "values (e.g. in -f, -r, or -c modes)"),
            _Switch(["--nofw", "nofw"],
                    "If --nofw is specified, bowtie2 will not attempt to align unpaired reads to the forward (Watson) "
                    "reference strand. In paired-end mode, --nofw and --norc pertain to the fragments; i.e. specifying "
                    "--nofw causes bowtie2 to explore only those paired-end configurations corresponding to fragments "
                    "from the reverse-complement (Crick) strand. Default: both strands enabled"),
            _Switch(["--norc", "norc"],
                    "If --norc is specified, bowtie2 will not attempt to align unpaired reads against the reverse-"
                    "complement Crick reference strand. In paired-end mode, --nofw and --norc pertain to the fragments;"
                    " i.e. specifying --nofw causes bowtie2 to explore only those paired-end configurations "
                    "corresponding to fragments from the reverse-complement (Crick) strand. Default: both strands"),
            _Switch(["--no-1mm-upfront", "no-1mm-upfront"],
                    "By default, Bowtie 2 will attempt to find either an exact or a 1-mismatch end-to-end alignment"
                    " for the read before trying the multiseed heuristic. Such alignments can be found very quickly,"
                    " and many short read alignments have exact or near-exact end-to-end alignments. However, this can "
                    "lead to unexpected alignments when the user also sets options governing the multiseed heuristic, "
                    "like -L and -N. For instance, if the user specifies -N 0 and -L equal to the length of the read, "
                    "the user will be surprised to find 1-mismatch alignments reported. This option prevents Bowtie 2 "
                    "from searching for 1-mismatch end-to-end alignments before using the multiseed heuristic, which "
                    "leads to the expected behavior when combined with options such as -L and -N. This comes at the "
                    "expense of speed"),
            _Switch(["--end-to-end", "end-to-end"],
                    "In this mode, Bowtie 2 requires that the entire read align from one end to the other, without any "
                    "trimming (or soft clipping) of characters from either end. The match bonus --ma always equals 0 in"
                    " this mode, so all alignment scores are less than or equal to 0, and the greatest possible "
                    "alignment score is 0. This is mutually exclusive with --local. --end-to-end is the default mode"),
            _Switch(["--local", "local"],
                    "In this mode, Bowtie 2 does not require that the entire read align from one end to the other. "
                    "Rather, some characters may be omitted (soft clipped) from the ends in order to achieve the "
                    "greatest possible alignment score. The match bonus --ma is used in this mode, and the best "
                    "possible alignment score is equal to the match bonus (--ma) times the length of the read. "
                    "Specifying --local and one of the presets (e.g. --local --very-fast) is equivalent to specifying "
                    "the local version of the preset (--very-fast-local). This is mutually exclusive with --end-to-end."
                    " --end-to-end is the default mode"),

            # Scoring Options
            _Option(["--score-min", "score_min"],
                    "Sets a function governing the minimum alignment score needed for an alignment to be considered "
                    "valid (i.e. good enough to report). This is a function of read length. For instance, specifying "
                    "L,0,-0.6 sets the minimum-score function f to f(x) = 0 + -0.6 * x, where x is the read length."
                    " See also: setting function options. The default in --end-to-end mode is L,-0.6,-0.6 "
                    "and the default in --local mode is G,20,8.",
                    checker_function=lambda value: re.match('^[CLSG],[-\d\.],[-\d\.]', value) is not None,
                    equate=False),
            _Option(["--ma", "ma"],
                    "Sets the match bonus. In --local mode <int> is added to the alignment score for each "
                    "position where a read character aligns to a reference character and the characters match. "
                    "Not used in --end-to-end mode. Default: 2.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--np", "np"],
                    "Sets penalty for positions where the read, reference, or both, contain an ambiguous "
                    "character such as N. Default: 1.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--rdg", "rdg"],
                    "Sets the read gap open (<int1>) and extend (<int2>) penalties. A read gap of length N gets"
                    " a penalty of <int1> + N * <int2>. Default: 5, 3.",
                    checker_function=lambda value: re.match('[-d.],[-d.]', value) is not None,
                    equate=False),
            _Option(["--rfg", "rfg"],
                    "Sets the reference gap open (<int1>) and extend (<int2>) penalties. A reference gap of "
                    "length N gets a penalty of <int1> + N * <int2>. Default: 5, 3.",
                    checker_function=lambda value: re.match('[-d.],[-d.]', value) is not None,
                    equate=False),
            _Option(["--mp", "mp"],
                    "Sets the maximum (MX) and minimum (MN) mismatch penalties, both integers. A number less "
                    "than or equal to MX and greater than or equal to MN is subtracted from the alignment score for "
                    "each position where a read character aligns to a reference character, the characters do not match,"
                    " and neither is an N. If --ignore-quals is specified, the number subtracted quals MX. "
                    "Otherwise, the number subtracted is MN + floor( (MX-MN)(MIN(Q, 40.0)/40.0) ) "
                    "where Q is the Phred quality value. Default: MX = 6, MN = 2.",
                    checker_function=lambda value: re.match('[-d.],[-d.]', value) is not None,
                    equate=False),

            # Reporting Options
            _Option(["-k", "k"],
                    "By default, bowtie2 searches for distinct, valid alignments for each read. When it finds a"
                    " valid alignment, it continues looking for alignments that are nearly as good or better. The best "
                    "alignment found is reported (randomly selected from among best if tied). Information about the "
                    "best alignments is used to estimate mapping quality and to set SAM optional fields, such as "
                    "AS:i and XS:i.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["-a", "a"],
                    "Like -k but with no upper limit on number of alignments to search for. "
                    "-a is mutually exclusive with -k."),

            # Effort Options
            _Option(["-D", "D"],
                    "Up to <int> consecutive seed extension attempts can fail before Bowtie 2 moves on, using"
                    " the alignments found so far. A seed extension fails if it does not yield a new best or a new "
                    "second-best alignment. This limit is automatically adjusted up when -k or -a are specified. "
                    "Default: 15.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["-R", "R"],
                    "<int> is the maximum number of times Bowtie 2 will re-seed reads with repetitive seeds. "
                    "When re-seeding, Bowtie 2 simply chooses a new set of reads (same length, same number of "
                    "mismatches allowed) at different offsets and searches for more alignments. A read is considered "
                    "to have repetitive seeds if the total number of seed hits divided by the number of seeds that "
                    "aligned at least once is greater than 300. Default: 2.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),

            # Paired-end options
            _Option(["--minins", "minins"],
                    "The minimum fragment length for valid paired-end alignments. E.g. if -I 60 is specified "
                    "and a paired-end alignment consists of two 20-bp alignments in the appropriate orientation with "
                    "a 20-bp gap between them, that alignment is considered valid (as long as -X is also satisfied). "
                    "A 19-bp gap would not be valid in that case. If trimming options -3 or -5 are also used, "
                    "the -I constraint is applied with respect to the untrimmed mates. The larger the difference "
                    "between -I and -X, the slower Bowtie 2 will run. This is because larger differences bewteen -I "
                    "and -X require that Bowtie 2 scan a larger window to determine if a concordant alignment exists. "
                    "For typical fragment length ranges (200 to 400 nucleotides), Bowtie 2 is very efficient. "
                    "Default: 0 (essentially imposing no minimum)",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--maxins", "maxins"],
                    "The maximum fragment length for valid paired-end alignments. E.g. if -X 100 is specified "
                    "and a paired-end alignment consists of two 20-bp alignments in the proper orientation with a "
                    "60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied). "
                    "A 61-bp gap would not be valid in that case. If trimming options -3 or -5 are also used, the "
                    "-X constraint is applied with respect to the untrimmed mates, not the trimmed mates. The larger "
                    "the difference between -I and -X, the slower Bowtie 2 will run. This is because larger differences "
                    "bewteen -I and -X require that Bowtie 2 scan a larger window to determine if a concordant "
                    "alignment exists. For typical fragment length ranges (200 to 400 nucleotides), "
                    "Bowtie 2 is very efficient. Default: 500",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["--fr", "fr"],
                    "The upstream/downstream mate orientations for a valid paired-end alignment against the "
                    "forward reference strand. E.g., if --fr is specified and there is a candidate paired-end "
                    "alignment where mate 1 appears upstream of the reverse complement of mate 2 and the fragment "
                    "length constraints (-I and -X) are met, that alignment is valid. Also, if mate 2 appears "
                    "upstream of the reverse complement of mate 1 and all other constraints are met, "
                    "that too is valid. --rf likewise requires that an upstream mate1 be reverse-complemented "
                    "and a downstream mate2 be forward-oriented. --ff requires both an upstream mate 1 and a "
                    "downstream mate 2 to be forward-oriented. "
                    "Default: --fr (appropriate for Illumina's Paired-end Sequencing Assay)."),

            _Switch(["--rf", "rf"],
                    "The upstream/downstream mate orientations for a valid paired-end alignment against the "
                    "forward reference strand. E.g., if --fr is specified and there is a candidate paired-end "
                    "alignment where mate 1 appears upstream of the reverse complement of mate 2 and the fragment "
                    "length constraints (-I and -X) are met, that alignment is valid. Also, if mate 2 appears "
                    "upstream of the reverse complement of mate 1 and all other constraints are met, "
                    "that too is valid. --rf likewise requires that an upstream mate1 be reverse-complemented "
                    "and a downstream mate2 be forward-oriented. --ff requires both an upstream mate 1 and a "
                    "downstream mate 2 to be forward-oriented. "
                    "Default: --fr (appropriate for Illumina's Paired-end Sequencing Assay)."),
            _Switch(["--ff", "ff"],
                    "The upstream/downstream mate orientations for a valid paired-end alignment against the "
                    "forward reference strand. E.g., if --fr is specified and there is a candidate paired-end "
                    "alignment where mate 1 appears upstream of the reverse complement of mate 2 and the fragment "
                    "length constraints (-I and -X) are met, that alignment is valid. Also, if mate 2 appears "
                    "upstream of the reverse complement of mate 1 and all other constraints are met, "
                    "that too is valid. --rf likewise requires that an upstream mate1 be reverse-complemented "
                    "and a downstream mate2 be forward-oriented. --ff requires both an upstream mate 1 and a "
                    "downstream mate 2 to be forward-oriented. "
                    "Default: --fr (appropriate for Illumina's Paired-end Sequencing Assay)."),
            _Switch(["--no-mixed", "no-mixed"],
                    "By default, when bowtie2 cannot find a concordant or discordant alignment for a pair, it "
                    "then tries to find alignments for the individual mates. This option disables that behavior."),

            _Switch(["--no-discordant", "no-discordant"],
                    "By default, bowtie2 looks for discordant alignments if it cannot find any concordant "
                    "alignments. A discordant alignment is an alignment where both mates align uniquely, "
                    "but that does not satisfy the paired-end constraints (--fr/--rf/--ff, -I, -X). "
                    "This option disables that behavior."),

            _Switch(["--dovetail", "dovetail"],
                    "If the mates dovetail, that is if one mate alignment extends past the beginning of the "
                    "other such that the wrong mate begins upstream, consider that to be concordant. See also: "
                    "Mates can overlap, contain or dovetail each other. Default: mates cannot dovetail "
                    "in a concordant alignment."),

            _Switch(["--no-contain", "no-contain"],
                    "If one mate alignment contains the other, consider that to be non-concordant. See also: "
                    "Mates can overlap, contain or dovetail each other. Default: a mate can contain "
                    "the other in a concordant alignment."),

            _Switch(["--no-overlap", "no-overlap"],
                    "If one mate alignment overlaps the other at all, consider that to be non-concordant. See "
                    "also: Mates can overlap, contain or dovetail each other. Default: mates can overlap in "
                    "a concordant alignment."),

            # SAM options
            _Switch(["--no-unal", "no-unal"],
                    "Suppress SAM records for reads that failed to align"),
            _Switch(["--no-hd", "no-hd"],
                    "Suppress SAM header lines (starting with"),
            _Switch(["--no-sq", "no-sq"],
                    "Suppress @SQ SAM header lines"),
            _Switch(["--omit-sec-seq", "omit-sec-seq"],
                    "When printing secondary alignments, Bowtie 2 by default will write out the SEQ and QUAL strings. "
                    "Specifying this option causes Bowtie 2 to print an asterix in those fields instead."),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)


if __name__ == '__main__':
    x = 0
    type
    pass
