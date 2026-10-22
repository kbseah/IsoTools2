from isotools.transcriptome import Transcriptome
from isotools._transcriptome_io import _add_sample_transcript
from pysam import FastaFile
import json


def test_add_orfs_lowercase_fasta():
    # regression test for #26: genome FASTA files commonly soft-mask repeats
    # with lowercase bases. The start/stop codons here are uppercase (so
    # find_orfs locates the ORF regardless), but the surrounding sequence
    # used for the kozak score lookup is lowercase -- this crashed with an
    # uncaught KeyError before add_orfs uppercased the sequence.
    isoseq = Transcriptome.from_reference("tests/data/lowercase_orf_example.gff3")
    gene = isoseq["GENE1"]
    with FastaFile("tests/data/lowercase_orf_example.fa") as genome_fh:
        gene.add_orfs(genome_fh, reference=True, minlen=0, max_5utr_len=12)

    orf = gene.ref_transcripts[0].get("ORF")
    assert (
        orf is not None
    ), "expected an ORF to be found despite lowercase genome sequence"
    genome_start, genome_end, orf_info = orf
    assert (genome_start, genome_end) == (12, 33), "unexpected ORF coordinates"
    assert orf_info["start_codon"] == "ATG"
    assert orf_info["stop_codon"] == "TAA"
    assert isinstance(
        orf_info["kozak"], float
    ), "kozak score should compute without error"


def test_add_orfs_prefer_annotated_init_false():
    # add_orfs prefers existing annotated start codon in ref_transcripts if
    # available. If `prefer_annotated_init=False`, find_orfs crashed with
    # AttributeError because of wrong default. Check that orfs are indeed
    # predicted when instructed to ignore reference
    isoseq = Transcriptome.from_reference("tests/data/zero_cds_example.ref.gtf")
    gene = isoseq["AAEL003046"]
    with open("tests/data/zero_cds_example.annot.json", "r") as fh:
        transcripts = json.load(fh)
    for tr in transcripts:
        _ = _add_sample_transcript(
            isoseq, tr, gene.chrom, fuzzy_junction=0, min_exonic_ref_coverage=0.25
        )
    with FastaFile("tests/data/zero_cds_example.fasta") as genome_fh:
        # max_5utr_len is default 0 in add_orfs but default 500 in add_orf_predictions
        gene.add_orfs(
            genome_fh, reference=False, prefer_annotated_init=False, max_5utr_len=500
        )
    assert "ORF" in gene.transcripts[0]


def test_add_orfs_prefer_annotated_init_true():
    # add_orfs added zero-length ORFs when `prefer_annotated_init=True`
    # because upstream AUGs (uAUGs) were also included in find_orfs output.
    isoseq = Transcriptome.from_reference("tests/data/zero_cds_example.ref.gtf")
    gene = isoseq["AAEL003046"]
    with open("tests/data/zero_cds_example.annot.json", "r") as fh:
        transcripts = json.load(fh)
    for tr in transcripts:
        _ = _add_sample_transcript(
            isoseq, tr, gene.chrom, fuzzy_junction=0, min_exonic_ref_coverage=0.25
        )
    with FastaFile("tests/data/zero_cds_example.fasta") as genome_fh:
        # max_5utr_len is ignored when prefer_annotated_init is true
        gene.add_orfs(genome_fh, reference=False, prefer_annotated_init=True)
    assert "ORF" in gene.transcripts[0]
    assert "ORF" not in gene.transcripts[1]


def test_import_ref_transcripts_split_start_codons():
    # Reference transcripts where start codon features were split across two
    # exons were imported with the wrong coordinate. Parser should treat split
    # start codon as a single feature
    isoseq = Transcriptome.from_reference("tests/data/split_start_codon_example.ref.gtf")
    gene = isoseq["AAEL010787"]
    genome_fh = "tests/data/split_start_codon_example.fasta"
    # _seq = gene.get_sequence(genome_fh, [0], reference=True, protein=True)
    # assert _seq is not None
    assert gene.ref_transcripts[0]["CDS"][0] == 1015

