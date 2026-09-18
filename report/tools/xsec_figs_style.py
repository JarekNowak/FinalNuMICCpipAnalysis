"""xsec_figs_style.py -- the style of the released cross-section figures, for the matplotlib figures.

White canvas, black axes and text, Okabe-Ito generator colours and the dashed line styles of
UnfolderNuMI (#0072B2, #009E73, #CC79A7, #D55E00 with styles 2, 7, 9, 3), black markers for the
(fake) data, grey for the fake-data truth.

Colour maps use ROOT's kBird, which is what UnfolderNuMI sets for every matrix figure of the
released note (gStyle->SetPalette(kBird), UnfolderNuMI.C:368 and :1249), so a map here reads like a
migration or covariance matrix there. BIRD is built from kBird's nine control points; against the
255-colour table dumped from this ROOT installation it agrees to 0.008 in RGB, which is the 8-bit
quantisation of that table.
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
_BIRD = [(0.2082, 0.1664, 0.5293), (0.0592, 0.3599, 0.8684), (0.0780, 0.5041, 0.8385),
         (0.0232, 0.6419, 0.7914), (0.1802, 0.7178, 0.6425), (0.5301, 0.7492, 0.4662),
         (0.8186, 0.7328, 0.3499), (0.9956, 0.7862, 0.1968), (0.9764, 0.9832, 0.0539)]
BIRD = LinearSegmentedColormap.from_list('bird', _BIRD)


def ink(rgba):
    """Readable text on a kBird cell: the palette runs from dark blue to yellow, so the label has to
    follow the background rather than sit at one fixed colour."""
    r, g, b = rgba[:3]
    return 'black' if 0.299 * r + 0.587 * g + 0.114 * b > 0.55 else 'white'


def apply():
    plt.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'axes.edgecolor': 'black', 'axes.labelcolor': 'black', 'text.color': 'black',
        'xtick.color': 'black', 'ytick.color': 'black', 'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.top': True, 'ytick.right': True, 'axes.grid': False, 'axes.spines.top': True,
        'axes.spines.right': True, 'axes.linewidth': 0.8, 'font.size': 10, 'legend.frameon': False,
        'lines.linewidth': 1.8, 'figure.dpi': 150})
