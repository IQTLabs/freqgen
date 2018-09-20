from bokeh.core.properties import value
from bokeh.io import show as _show
from bokeh.io import output_file
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
from bokeh.transform import dodge
from bokeh.palettes import Set2_3

def visualize(k_mers,
              target_freqs,
              optimized_freqs,
              original_freqs=None,
              title="Freqgen Optimization Results",
              plot_height=400,
              plot_width=1200,
              show=True):
    '''Creates a visualization of the results of a Freqgen optimization.

    Note:
        Currently, this function does not support the visualization of codon
        optimization *and* *k*-mer optimization simultaneously.

    Args:
        k_mers (list): A list of the *k*-mers to use as the labels for the *x*-axis.
        target_freqs (list): A list of the target frequencies in the same order as the `k_mers` argument.
        optimized_freqs (list): A list of the resultant frequencies in the same order as the `k_mers` argument.
        original_freqs (list, optional): A list of the original frequencies in the same order as the `k_mers` argument.
        title (str): A title to use for the graph. Defaults to "Freqgen Optimization Results".
        plot_height (int): The height for the graph. Defaults to 400.
        plot_width (int): The width for the graph. Defaults to 1200.
        show (bool): Whether to show the plot or simply return it. Defaults to True.

    Returns:
        bokeh.plotting.figure.Figure: A Bokeh figure containing the bar graph.
    '''

    output_file("freqgen.html")

    categories = ['Original', 'Target', 'Optimized']

    data = {'k_mers': k_mers,
            'Target': target_freqs,
            'Optimized': optimized_freqs}

    if not isinstance(original_freqs, type(None)):
        offset = 0.25
        data["Original"] = original_freqs
    else:
        offset = 0.15

    source = ColumnDataSource(data=data)

    p = figure(x_range=k_mers, plot_height=plot_height, plot_width=plot_width, title=title)

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
    p.xaxis.axis_label = 'k-mer'

    _show(p)
    return p
