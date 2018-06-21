#!/usr/bin/env python
"""
Example command line interfacer for HIST feasibility
Michael Hirsch
"""
import matplotlib
matplotlib.use('Agg')
print(matplotlib.get_backend())
from matplotlib.pyplot import show  # noqa: E402
import seaborn as sns  # noqa: E402
sns.color_palette("cubehelix")
sns.set(context='paper', style='whitegrid', font_scale=2,
        rc={'image.cmap': 'cubehelix_r'})
#
from histfeas import userinput  # noqa: E402
from histfeas.main_hist import doSim  # noqa: E402


if __name__ == '__main__':
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    P = userinput()
    doSim(P)

    show()
