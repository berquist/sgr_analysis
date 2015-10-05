#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import math
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def norm_pdf(x, mu, sigma):
    """Return the normal probability distribution function."""

    prefac = 1 / (sigma * math.sqrt(2 * math.pi))
    expt = np.exp(-(x - mu)**2 / (2 * sigma**2))

    return prefac * expt

### example normal distributions
### https://en.wikipedia.org/wiki/Normal_distribution

fig, ax = plt.subplots()

x = np.linspace(-5, 5, 300)

ax.plot(x, norm_pdf(x, mu=0, sigma=math.sqrt(0.2)), label=r'$\mu=0,\,\sigma^2=0.2$', color='blue')
ax.plot(x, norm_pdf(x, mu=0, sigma=math.sqrt(1.0)), label=r'$\mu=0,\,\sigma^2=1.0$', color='red')
ax.plot(x, norm_pdf(x, mu=0, sigma=math.sqrt(5.0)), label=r'$\mu=0,\,\sigma^2=5.0$', color='yellow')
ax.plot(x, norm_pdf(x, mu=-2, sigma=math.sqrt(0.5)), label=r'$\mu=-2,\,\sigma^2=0.5$', color='green')

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$\varphi_{\mu,\sigma^2}(x)$')

ax.legend(fancybox=True, loc='best', framealpha=0.50)
fig.savefig('norm_pdf.pdf', bbox_inches='tight')

plt.close(fig)
