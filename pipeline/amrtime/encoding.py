import os
import subprocess
import pandas as pd
import numpy as np
from amrtime import parsers


class Homology():
    """
    Generate a read encoding
    """

    def __init__(self, simulated_reads, reference_db):

        self.reads = simulated_reads
        self.db = reference_db

        if not os.path.exists('training_data/diamond.out6') :
            self.alignment_fp = self.run_alignment(self.reads, self.db)
        else:
            self.alignment_fp = 'training_data/diamond.out6'

    def make_database(self, db):
        """
        Make the diamond search db if it hasn't been made
        """
        if not os.path.exists(db + '.dmnd'):
            subprocess.check_call('diamond makedb --in {0} --db {0}'.format(db),
                    shell=True)

    def run_alignment(self, reads, db):
        if not os.path.exists('training_data/diamond.out6'):
            self.make_database(self.db)
            subprocess.check_call('diamond blastx --db {} --out training_data/diamond.out6 --outfmt 6 --threads 2 --query {} --more-sensitive'.format(db, reads), shell=True)
        return 'training_data/diamond.out6'

    def encode(self, card, dissimilarity=False):

        alignment_fh = open(self.alignment_fp)
        reads_fh = open(self.reads)

        amr_family_field = {family: ix for ix, family in enumerate(set(card.aro_to_gene_family.values()))}
        field_to_family = {v: k for k, v in amr_family_field.items()}

        gene_family_encoding = {}
        for ix, read in enumerate(reads_fh):
            if ix % 4 == 0:
                read_acc = read.strip().replace('@gb', 'gb')
                sim_vector = np.zeros(len(amr_family_field))
                gene_family_encoding.update({read_acc : sim_vector})

        for alignment in alignment_fh:
            alignment = alignment.strip().split('\t')
            query_acc = alignment[0]
            alignment_aro = alignment[1].split('|')[2]

            amr_family = card.aro_to_gene_family[alignment_aro]
            field = amr_family_field[amr_family]
            gene_family_encoding[query_acc][field] += float(alignment[11])

        alignment_fh.close()
        reads_fh.close()

        family_index = [field_to_family[x] for x in range(len(field_to_family))]
        df = pd.DataFrame(gene_family_encoding)
        df = df.transpose()
        df.columns = family_index

        if dissimilarity:
            df = df.applymap(lambda x: 1-x)
            return df
        else:
            df = df.divide(df.max(axis=1), axis=0)
            return df

class TNF():

    def __init__(self, metagenome_fp):
        self.metagenome_fp = metagenome_fp

    def encode(self):
        tnf = ["".join(x) for x in itertools.product(['A', 'T', 'G', 'C', 'N'], repeat=4)]
        tnf_encoder = {v:k for k,v in enumerate(tnf)}

        X = []
        with open(self.metagenome_fp) as fh:
            for ix, line in enumerate(fh):
                if ix % 4 == 1:
                    clean_seq = re.sub("[M,X,R,S,Y,K]", 'N', line.strip())
                    encoded_seq = self.read_encode(clean_seq, tnf_encoder)
                    X.append(encoded_seq)
        return np.vstack(X)

    def read_encode(self, seq, tnf):
        """
        k-mer decomposition might be simplest although might have to be
        done at the file level in order to represent all k-mers
        """
        encoded_vec = np.zeros(len(tnf))
        for tetranucleotide in self.window(4):
            encoded_vec[tnf[tetranucleotide]] += 1
        return encoded_vec

    def window(self, window_size):
        "Returns a sliding window (of width w) over data from the iterable"
        "   s -> (s0,s1,...s[w-1]), (s1,s2,...,sw), ...                   "
        it = iter(self.seq)
        result = tuple(itertools.islice(it, window_size))
        if len(result) == window_size:
            yield "".join(result)
        for elem in it:
            result = result[1:] + (elem,)
            yield "".join(result)

# simple kmer embedding
class KMer():
    pass

# dna2vec kmer embedding
class KMer_embedding():
    pass
