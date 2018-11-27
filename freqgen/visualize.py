import math

from bokeh.core.properties import value
from bokeh.io import output_file, save
from bokeh.io import show as _show
from bokeh.models import ColumnDataSource
from bokeh.palettes import Set2_3
from bokeh.plotting import figure
from bokeh.transform import dodge


def visualize(k_mers,
              target_freqs,
              optimized_freqs,
              original_freqs=None,
              title="Freqgen Optimization Results",
              plot_height=400,
              plot_width=1200,
              show=True,
              filepath="freqgen.html",
              codons=False):
    '''Creates a visualization of the results of a Freqgen optimization.

    Note:
        Currently, this function does not support the visualization of codon
        optimization *and* *k*-mer optimization simultaneously.

    Args:
        k_mers (list): A list of the *k*-mers to use as the labels for the *x*-axis.
        target_freqs (list): A list of the target frequencies in the same order as the `k_mers` argument.
        optimized_freqs (list): A list of the resultant frequencies in the same order as the `k_mers` argument.
        original_freqs (list, optional): A list of the original frequencies in the same order as the `k_mers` argument.
        title (str, optional): A title to use for the graph. Defaults to "Freqgen Optimization Results".
        plot_height (int, optional): The height for the graph. Defaults to 400.
        plot_width (int, optional): The width for the graph. Defaults to 1200.
        show (bool, optional): Whether to show the plot or simply return it. Defaults to True.
        filepath (str, optional): The output filepath. Defaults to "freqgen.html".
        codons (bool, optional): Whether codons are included in the input vectors. If they are, the *x*-axis legend will updated accordingly.

    Note:
        Codons must be denoted with a ``*`` in the ``k_mers`` argument. For example, the codon GAG should be passed as ``GAG*``

    Returns:
        bokeh.plotting.figure.Figure: A Bokeh figure containing the bar graph.
    '''

    # validate that all the codons that should be there are there
    if codons and len([k_mer for k_mer in k_mers if k_mer.endswith("*")]) != 64:
        raise ValueError("You appear to be passing an incomplete list of codons.")

    codons_only = False
    if all([k_mer.endswith("*") for k_mer in k_mers]) and codons:
        k_mers = [k_mer[:-1] for k_mer in k_mers]
        codons_only = True

    output_file(filepath)

    categories = ['Original', 'Target', 'Optimized']

    data = {'k_mers': k_mers,
            'Target': target_freqs,
            'Optimized': optimized_freqs}

    # adjust spacing depending on if there's three bars or two
    if not isinstance(original_freqs, type(None)):
        offset = 0.25
        data["Original"] = original_freqs
    else:
        offset = 0.15

    # identify the greatest y value for setting bounds
    y_max = max((max(target_freqs), max(optimized_freqs)))
    if not isinstance(original_freqs, type(None)):
        y_max = max((max(original_freqs), y_max))

    source = ColumnDataSource(data=data)

    p = figure(x_range=k_mers,
               plot_height=plot_height,
               plot_width=plot_width,
               title=title,
               y_range=(0, 1.2 * y_max))

    if not isinstance(original_freqs, type(None)):
        p.vbar(x=dodge('k_mers', -0.25, range=p.x_range), top='Original', width=0.2, source=source,
               color=Set2_3[0], legend=value("Original"))

    p.vbar(x=dodge('k_mers',  0.0 if not isinstance(original_freqs, type(None)) else -offset,  range=p.x_range), top='Target', width=0.2, source=source,
           color=Set2_3[1], legend=value("Target"))

    p.vbar(x=dodge('k_mers',  offset, range=p.x_range), top='Optimized', width=0.2, source=source,
           color=Set2_3[2], legend=value("Optimized"))

    p.x_range.range_padding = 0.05
    p.xgrid.grid_line_color = None
    p.legend.location = "top_right"
    p.legend.orientation = "horizontal"
    p.legend.click_policy = "hide"

    # decide the x-axis label
    if codons_only:
        p.xaxis.axis_label = "codon"
    elif codons:
        p.xaxis.axis_label = "k-mer (* denotes codon)"
    else:
        p.xaxis.axis_label = 'k-mer'

    if len(k_mers) >= 32:
        p.xaxis.major_label_orientation = math.pi / 2

    p.yaxis.axis_label = 'frequency'
    if show:
        _show(p)
    else:
        save(p, filename=filepath)
    return p
