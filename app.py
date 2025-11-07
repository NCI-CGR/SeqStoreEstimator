from shiny import reactive
from shiny.express import input, render, ui
import shinyswatch
from shinywidgets import render_plotly
from src.Functions import estimate_bam_size_from_nreads,file_size_converter, plot_cummulative_cost_over_years, reads_to_bases
from faicons import icon_svg



ui.page_opts(title=ui.span(
        icon_svg("calculator"),
        " SeqStoreEstimator"
    ), fillable=True, theme=shinyswatch.theme.simplex,navbar_color="primary")


with ui.sidebar(title=ui.span(icon_svg("gears"), "Parameters")):

    ui.input_numeric(
            id="num_reads",
            label="Number of Reads",
            min=0,
            max=6_000_000_000,
            value=3_771_780_000_000,
            step=1_000_000,
        )

    # with ui.card():
    #     ui.input_numeric(
    #     id="num_bases",
    #     label="Number of Bases",
    #     min=0,
    #     max=6_000_000_000,
    #     value=6_000_000_000,
    #     step=1_000_000,
    #     ),

    #     ui.input_select(
    #         id="unit",
    #         label="Unit",
    #         choices=["bp", "Kb", "Mb", "Gb"],
    #         selected="Gb",
    #     )

    with ui.card():
        ui.input_numeric(
            id="read_length",
            label="Read Length",
            min=0,
            max=800,
            value=150,
            step=1,
        )
        ui.input_numeric(
            id="supplementary_alignments",
            label="Supplementary Alignments",
            min=0,
            max=1,
            value=0.1,
        )
        ui.input_numeric(
            id="mapped",
            label="Mapped Reads",
            min=0,
            max=1,
            value=0.9,
        )
    with ui.card():
        ui.input_numeric(
        id="bam_compression_ratio",
        label="BAM Compression Ratio",
        min=0.0,
        max=1.0,
        value=0.15,
        )

        ui.input_numeric(
            id="cram_compression_ratio",
            label="CRAM Compression Ratio",
            min=0.0,
            max=1.0,
            value=0.30)
    
    ui.input_select(
        id="output_format",
        label="Output Format",
        choices=["BAM", "CRAM"],
        selected="CRAM")

    ui.input_numeric(
        id="cost_per_month_per_gb",
        label="Cost per Month per GB ($)",
        value=0.0064,
        step=0.0001
    )

with ui.layout_column_wrap(fill=False):
    with ui.value_box(theme="warning", showcase=icon_svg("file-zipper", style="solid")):
        @render.text
        def estimated_disk_usage_title():
            return f"Estimated {input.output_format()} File Size"
        @render.text
        def estimated_bam_size():
            bam_size_string = file_size_converter(estimated_bam_size_bytes())
            return f"{bam_size_string}"
    
    with ui.value_box(theme="warning", showcase=icon_svg("sack-dollar", style="solid")):
        "Estimated Monthly Storage Cost"
        @render.text
        def estimated_monthly_cost_display():
            return f"${estimated_monthly_cost():.2f} per month"

with ui.card(full_screen=True):
    @render_plotly
    def cummulative_cost_chart():
        return plot_cummulative_cost_over_years(
            monthly_cost=estimated_monthly_cost(),
            years=5
        )

with ui.layout_column_wrap(fill=False):
    with ui.value_box(showcase_layout="top right", theme="light"):
        "Number of Reads"
        @render.text
        def num_reads_display():
            return f"{input.num_reads()}"


    with ui.value_box(showcase_layout="top right",  theme="light"):
        "Number of Bases"
        @render.text
        def num_bases_display():
            num_bases = reads_to_bases(input.num_reads(), input.read_length())
            return f"{num_bases} bases"

    with ui.value_box(showcase_layout="top right",  theme="light"):
        "Read Length"
        @render.text
        def read_length_display():
            return f"{input.read_length()} nt"


## calculations
@reactive.Calc
def estimated_bam_size_bytes() -> float:
    return estimate_bam_size_from_nreads(n_reads = input.num_reads(), 
                                                read_len = input.read_length(),
                                                bam_compression_ratio = input.bam_compression_ratio(),
                                                supplementary_alignments = input.supplementary_alignments(),
                                                percent_mapped = input.mapped(),
                                                output_format = input.output_format(),
                                                cram_compression_ratio = input.cram_compression_ratio()
                                                )
@reactive.Calc
def estimated_monthly_cost() -> float:
    size_in_gb = estimated_bam_size_bytes() / (1024 ** 3)
    return size_in_gb * input.cost_per_month_per_gb()