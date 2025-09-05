 #!/usr/bin/env python
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('font', size=4)
# Load from a file
tree = Phylo.read("RAxML_bestTree.R3860222561-3205-16-B08_core", "newick")
#tree = Phylo.parse('RAxML_bestTree.R3860222561-3205-16-B08_core', 'phyloxml').next()
# Or load from a string
# tree = Phylo.read(io.StringIO("((A:0.2,B:0.3):0.3,C:0.5);"), "newick")

fig = plt.figure(figsize=(10, 8))
axes = fig.add_subplot(1, 1, 1)

Phylo.draw(tree, do_show=False, axes=axes)

# Change the thickness of all lines in the plot
for line in axes.get_lines():
    line.set_linewidth(0.8)  # Set your desired thickness here

plt.savefig("output_tree.png", dpi=300)
