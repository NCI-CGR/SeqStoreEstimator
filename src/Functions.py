
import plotly.graph_objs as go

def estimate_bam_size_from_nreads(n_reads: int,
                                  read_len: int = 150,
                                  bam_compression_ratio: float = 0.15,
                                  cram_compression_ratio: float = 0.3,
                                  output_format: str = "CRAM",
                                  supplementary_alignments: float = 0.1,
                                  percent_mapped: float = 0.9) -> float:
    """
    Estimate the disk usage in bytes of a BAM/CRAM file from the number of reads.
    
    Parameters
    ----------
    n_reads : int
        Number of sequencing reads.
    read_len : int, default=150
        Read length in bases.
    bam_compression_ratio : float, default=0.15
        Compression ratio for BAM format.
    cram_compression_ratio : float, default=0.3
        Additional compression ratio applied when output format is CRAM.
    output_format : str, default="CRAM"
        Output format ("BAM" or "CRAM").
    supplementary_alignments : float, default=0.1
        Proportion of reads with supplementary alignments.
    percent_mapped : float, default=0.9
        Proportion of reads that are mapped.
        
    Returns
    -------
    float
        Estimated file size in bytes.
    """
    bam_fixed_fields = {
        "block_size":  {"type": "int32_t",  "bytes": 4, "description": "BAM records starts with block size"},
        "refID":       {"type": "int32_t",  "bytes": 4, "description": "Reference sequence ID (−1 if unmapped)"},
        "pos":         {"type": "int32_t",  "bytes": 4, "description": "0-based leftmost coordinate (POS−1)"},
        "l_read_name": {"type": "uint8_t",  "bytes": 1, "description": "Length of read name including NUL"},
        "mapq":        {"type": "uint8_t",  "bytes": 1, "description": "Mapping quality (MAPQ)"},
        "bin":         {"type": "uint16_t", "bytes": 2, "description": "BAI index bin"},
        "n_cigar_op":  {"type": "uint16_t", "bytes": 2, "description": "Number of CIGAR operations"},
        "flag":        {"type": "uint16_t", "bytes": 2, "description": "Bitwise SAM FLAG"},
        "l_seq":       {"type": "uint32_t", "bytes": 4, "description": "Length of the sequence"},
        "next_refID":  {"type": "int32_t",  "bytes": 4, "description": "Reference ID of the next segment"},
        "next_pos":    {"type": "int32_t",  "bytes": 4, "description": "0-based leftmost position of the next segment"},
        "tlen":        {"type": "int32_t",  "bytes": 4, "description": "Template length (TLEN)"},
    }

    cigar_mix = {
        "perfect": {"percent": 0.6, "bytes": 4},
        "softclip": {"percent": 0.35, "bytes": 8},
        "indel": {"percent": 0.04, "bytes": 16},
        "complex": {"percent": 0.01, "bytes": 32},
    }

    aux_tag_types = {
        "integer": {"bytes": 7, "description": "Integer tags (e.g. NM:i, AS:i, MQ:i)"},
        "float": { "bytes": 8, "description": "Float tags (e.g. sd:f)"},
        "string": {"bytes": 16, "description": "String tags (e.g. RG:Z, PG:Z, MC:Z, RX:Z)"},
        "character": { "bytes": 4, "description": "Character tags (rare)"},
        "array_short_integer": { "bytes": 8 + 2 * 150, "description": "Short integer arrays (e.g. cd:B:s, ce:B:s)"},
        "array_integer": { "bytes": 8 + 4 * 0, "description": "Integer arrays"},
        "array_float": {"bytes": 8 + 4 * 0, "description": "Float arrays"},
    }

    aux_tags = {
        "mapped_reads": {"bytes" : aux_tag_types["integer"]["bytes"] * 4 * percent_mapped,
                               "description": "PAS, XS, AM, SM tags for primary alignments"},
        "core" : {"bytes": aux_tag_types["integer"]["bytes"] * 6 +  aux_tag_types["string"]["bytes"] * 4,
                  "description": "Core auxiliary tags present in all reads NM, MQ, UQ, PQ, etc"},
    }   


    bam_variable_fields = {
        "cigar":   {"bytes":  sum(v["percent"] * v["bytes"] for v in cigar_mix.values()), "description": "CIGAR operations (array of uint32_t, 4 bytes each), assumed 5 ops"},
        "read_name":     {"bytes": 35, "description": "Read name (l_read_name bytes), assumed 35 chars + NUL"},
        "aux_tags":     {"bytes": sum(v["bytes"] for v in aux_tags.values()), "description": "Auxiliary tags, assumed 45 bytes"},
        "SA_tag":      {"bytes": 120 * supplementary_alignments, "description": "SA auxiliary tag for supplementary alignments"},
    }


    total_fixed_bytes = sum(f["bytes"] for f in bam_fixed_fields.values())
    total_variable_bytes = sum(f["bytes"] for f in bam_variable_fields.values())


    bytes_quality: int = read_len
    bytes_seq: int = (read_len + 1) // 2  # 4-bit encoding
    bytes_per_read = total_fixed_bytes + total_variable_bytes + bytes_quality + bytes_seq
    total_bytes = n_reads * bytes_per_read
    print( "DEBUG: Total bytes before compression:", total_bytes)
    total_bytes *= bam_compression_ratio
    if output_format.upper() == "CRAM":
        total_bytes *= cram_compression_ratio 

    return total_bytes 


def estimate_fastqz_size_from_nreads(
    n_reads: int,
    read_len: int = 150,
    gzip_compression_ratio: float = 0.25,
    pe: bool = True,
    read_name_len: int = 20,
) -> float:
    """
    Estimate the gzipped FASTQ (.fastq.gz) size in bytes from number of reads.
    
    Parameters
    ----------
    n_reads : int
        Number of sequencing fragments (read pairs if pe=True, else single-end reads).
    read_len : int, default=150
        Read length in bases.
    gzip_compression_ratio : float, default=0.25
        Compression ratio for gzip (~0.22–0.28 typical for Illumina data).
    pe : bool, default=True
        If True, counts both R1 and R2 reads.
    read_name_len : int, default=20
        Length of the read name (without '@' or newline).
    """

    fastq_record = {
        "read_name": {
            "bytes": read_name_len + 2,
            "description": "@readname plus newline",
        },
        "sequence": {
            "bytes": read_len + 1,
            "description": "Sequence bases plus newline",
        },
        "separator": {
            "bytes": 2,
            "description": "'+' separator line plus newline",
        },
        "qualities": {
            "bytes": read_len + 1,
            "description": "Quality scores plus newline",
        },
    }

    bytes_per_read = sum(f["bytes"] for f in fastq_record.values())

    if pe:
        bytes_per_read *= 2

    total_bytes = n_reads * bytes_per_read * gzip_compression_ratio
    return total_bytes


def file_size_converter(size_in_bytes: float) -> str:
    """
    Convert file size in bytes to a human-readable string with appropriate units.
    
    Parameters
    ----------
    size_in_bytes : float
        File size in bytes.
        
    Returns
    -------
    str
        Human-readable file size string.
    """
    if size_in_bytes < 0:
        raise ValueError("File size cannot be negative.")

    units = ["B", "KB", "MB", "GB", "TB", "PB"]
    index = 0

    while size_in_bytes >= 1024 and index < len(units) - 1:
        size_in_bytes /= 1024
        index += 1

    return f"{size_in_bytes:.2f} {units[index]}"


def bases_to_reads(num_bases: float, read_length: int) -> int:
    """
    Convert number of bases to number of reads based on read length.
    
    Parameters
    ----------
    num_bases : float
        Total number of bases.
    read_length : int
        Length of each read in bases.
        
    Returns
    -------
    int
        Estimated number of reads.
    """
    if read_length <= 0:
        raise ValueError("Read length must be a positive integer.")
    
    return int(num_bases / read_length)

def reads_to_bases(num_reads: int, read_length: int) -> int:
    """
    Convert number of reads to number of bases based on read length.
    
    Parameters
    ----------
    num_reads : int
        Total number of reads.
    read_length : int
        Length of each read in bases.
        
    Returns
    -------
    int
        Estimated number of bases.
    """
    if read_length <= 0:
        raise ValueError("Read length must be a positive integer.")
    
    return num_reads * read_length

def to_bases(bases: int, unit: str) -> float:
    """
    Convert bases to a specific unit (bp, kb, Mb, Gb).
    """
    unit_multipliers = {
        "bp": 1,
        "Kb": 1e3,
        "Mb": 1e6,
        "Gb": 1e9,
    }

    if unit not in unit_multipliers:
        raise ValueError(f"Unknown unit: {unit}")

    return bases / unit_multipliers[unit]


def plot_cummulative_cost_over_years(monthly_cost: float, years: int = 5) -> dict:
    """
    Generate data for plotting monthly storage cost over a number of years.
    
    Parameters
    ----------
    monthly_cost : float
        Monthly storage cost in dollars.
    years : int, default=5
        Number of years to project.
        
    Returns
    -------
    fig : plotly.graph_objs.Figure
        Plotly figure object for the monthly cost over years.
    """
    months = list(range(1, years * 12 + 1))
    cumulative_cost = [monthly_cost * month for month in months]


    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=months,
        y=cumulative_cost,
        mode='lines+markers',
        name='Cumulative Cost'
    ))
    fig.update_layout(
        title="Cumulative Storage Cost Over 5 Years",
        xaxis_title="Month",
        yaxis_title="Cumulative Cost (USD)",
        template="plotly_white"
    )
    return fig
