"""xsec_figs_style.py -- the style of the released cross-section figures, for the matplotlib figures.

White canvas, black axes and text, Okabe-Ito generator colours and the dashed line styles of
UnfolderNuMI (#0072B2, #009E73, #CC79A7, #D55E00 with styles 2, 7, 9, 3), black markers for the
(fake) data, grey for the fake-data truth. Sequential maps use one hue; ratio maps diverge about a
neutral grey.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

GEN = [('h_gen_GENIE', 'GENIE', '#0072B2', (0, (6, 3))),
       ('h_gen_GiBUU', 'GiBUU', '#009E73', (0, (1, 1.6))),
       ('h_gen_NEUT', 'NEUT', '#CC79A7', (0, (7, 2.5, 1.5, 2.5))),
       ('h_gen_NuWro', 'NuWro', '#D55E00', (0, (3, 2)))]
TRUTH, TUNE, DATA = '#666666', '#56B4E9', 'black'
SEQ = LinearSegmentedColormap.from_list('seq', ['#ffffff', '#cfe3f2', '#8ec4e4', '#3b8fc4', '#0072B2', '#004c78'])
DIV = LinearSegmentedColormap.from_list('div', ['#7f3000', '#D55E00', '#f0c9ae', '#e8e8e8', '#a8cfe6', '#0072B2', '#00405f'])


def apply():
    plt.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'axes.edgecolor': 'black', 'axes.labelcolor': 'black', 'text.color': 'black',
        'xtick.color': 'black', 'ytick.color': 'black', 'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.top': True, 'ytick.right': True, 'axes.grid': False, 'axes.spines.top': True,
        'axes.spines.right': True, 'axes.linewidth': 0.8, 'font.size': 10, 'legend.frameon': False,
        'lines.linewidth': 1.8, 'figure.dpi': 150})
