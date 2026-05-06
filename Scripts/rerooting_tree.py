from ete3 import Tree

t = Tree(r"results/concatenated.fa.treefile")

t.set_outgroup("rhabditophanes_kr3021")

t.write(outfile=r"results/rooted.treefile")
