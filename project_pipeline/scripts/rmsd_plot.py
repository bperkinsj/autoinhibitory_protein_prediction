import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(snakemake.input[0], sep='\t').astype('object')

fig, ax = plt.subplots()
df.plot(kind='scatter', y='2.0_comp', x='2.0_aligned', ax=ax)
ax.set_xlim([0, 25])
ax.set(title='Rmsds of Proteins Aligned at Domain', xlabel='Domain RMSD', ylabel='Tail RMSD')

plt.savefig(snakemake.output[0])
