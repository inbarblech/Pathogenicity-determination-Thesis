class Mutation:
    def __init__(self, gene_name, gene_id, wt_aa, pos, mut_aa):
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.wt_aa = wt_aa
        self.pos = pos
        self.mut_aa = mut_aa

    def __repr__(self):
        return f"{self.wt_aa}{self.pos}{self.mut_aa}"

    def __str__(self):
        return f"{self.gene_name}_{self.gene_id}_{self.wt_aa}{self.pos}{self.mut_aa}"