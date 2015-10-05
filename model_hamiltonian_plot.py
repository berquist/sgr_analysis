#!/usr/bin/env python

from __future__ import print_function

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

file_old = "model_hamiltonian_frequencies_box.txt"
file_new = "model_hamiltonian_frequencies.txt"

with open(file_old) as fh_old:
    file_contents_old = fh_old.read()
with open(file_new) as fh_new:
    file_contents_new = fh_new.read()

def frequencies_from_file_contents(file_contents):
    return [float(line.split()[3].replace(']', ''))
            for line in file_contents.splitlines()[:-2]]

frequencies_old = frequencies_from_file_contents(file_contents_old)
frequencies_new = frequencies_from_file_contents(file_contents_new)

fig, ax = plt.subplots()
ax.plot(frequencies_old, label='old', color='blue')
ax.set_title('old')
ax.set_xlabel('snapshot #')
ax.set_ylabel(r'$\omega_{\mathrm{antisymmetric}}$ (cm$^{-1}$)')
fig.savefig('model_hamiltonian_frequencies_old.pdf', bbox_inches='tight')
plt.close(fig)

fig, ax = plt.subplots()
ax.plot(frequencies_new, label='new', color='orange')
ax.set_title('new')
ax.set_xlabel('snapshot #')
ax.set_ylabel(r'$\omega_{\mathrm{antisymmetric}}$ (cm$^{-1}$)')
fig.savefig('model_hamiltonian_frequencies_new.pdf', bbox_inches='tight')
plt.close(fig)

fig, ax = plt.subplots()
ax.plot(frequencies_old, label='old', color='blue')
ax.plot(frequencies_new, label='new', color='orange')
ax.set_title('combined')
ax.set_xlabel('snapshot #')
ax.set_ylabel(r'$\omega_{\mathrm{antisymmetric}}$ (cm$^{-1}$)')
ax.legend(loc='best', fancybox=True, framealpha=0.50)
fig.savefig('model_hamiltonian_frequencies_combined.pdf', bbox_inches='tight')
plt.close(fig)
