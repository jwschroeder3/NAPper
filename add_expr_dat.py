import pandas
import sys

# add expression data from a file into a data frame

in_table=sys.argv[1]
out_table=sys.argv[2]
expr_file=sys.argv[3]

fullframe=pandas.read_csv(in_table, sep="\t")
gene_names=fullframe["refseq"]
gene_exprs = []

for g in gene_names:

    this_exp = []

    with open(expr_file) as f:
        for l in f:
            for gn in g.split(";"):
                if gn.strip() == "":
                    continue

                if l.count(gn) > 0:
                    linearr=l.rstrip().split()
                    this_exp.append(float(linearr[1]))

    if len(this_exp) == 0:
        print("Warning: couldn't find a hit for %s" % g)
        gene_exprs.append(None)

    if len(this_exp) == 1:
        gene_exprs.append(this_exp[0])

    if len(this_exp) > 1:
        gene_exprs.append(this_exp[0])
        print("Warning: found multiple hits for %s" % g)

fullframe["%s" % expr_file] = gene_exprs
fullframe.to_csv(out_table, sep="\t", index=False)



