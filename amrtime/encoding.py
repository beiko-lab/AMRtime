import os
from tqdm import tqdm
import sys
import subprocess
import pandas as pd
import numpy as np
from amrtime import parsers
import math
import itertools
import re

from sklearn.preprocessing import normalize

class Homology():
    """
    Generate a read encoding
    """
    def __init__(self, simulated_reads, data_type, card, tool):

        self.reads = simulated_reads
        self.data_type = data_type

        self.db = 'training_data/card_proteins.faa'
        card.write_proteins(self.db)

        if tool == 'DIAMOND':
            if not os.path.exists(f'training_data/{data_type}_diamond.out6') :
                self.alignment_fp = self.run_diamond_alignment(self.reads, self.db)
            else:
                print(f"training_data/{data_type}_diamond.out6 already exists so re-using, use --redo to rebuild")
                self.alignment_fp = f'training_data/{data_type}_diamond.out6'
        elif tool == 'MMSEQS2':
            if not os.path.exists(f'training_data/{data_type}_mmseqs.out6') :
                self.alignment_fp = self.run_alignment(self.reads, self.db)
            else:
                print(f"training_data/{data_type}_mmseqs.out6 already exists so re-using, use --redo to rebuild")
                self.alignment_fp = f'training_data/{data_type}_mmseqs.out6'

    def run_diamond_alignment(self, reads, db):
        """
        Perform a DIAMOND BLASTX search to gather homology data
        """
        # build database if it doesn't exist
        if not os.path.exists(db + '.dmnd'):
            subprocess.check_call('diamond makedb --in {0} --db {0}'.format(db),
        shell=True)

        # run alignment
        subprocess.check_call(f'diamond blastx --db {db} --out training_data/{self.data_type}_diamond.out6 --outfmt 6 --threads 2 --query {reads} --more-sensitive', shell=True)
        return f'training_data/{self.data_type}_diamond.out6'

    def run_mmseqs_alignment(self, reads, db):
        subprocess.check_call(f'mmseqs easy-search {db} {reads} training_data/{self.data_type}_mmseqs.out6 /tmp', shell=True)


    def encode(self, card, metric, norm=False, dissimilarity=False):
        alignment_fh = open(self.alignment_fp)
        reads_fh = open(self.reads)

        # build reference to get the correct row for each AMR family
        if self.data_type == 'family':
            label_to_field = {family: ix for ix, family in enumerate(set(card.aro_to_gene_family.values()))}
            field_to_label = {v: k for k, v in label_to_field.items()}

        # build the same for the aros
        elif self.data_type == 'gene':
            label_to_field = {aro: ix for ix, aro in enumerate(set(card.aro_to_gene_family.keys()))}
            field_to_label = {v: k for k, v in field_to_label.items()}

        # read input fastq and initialise an empty vectors for each read
        # and store in a dictionary of read_acc: vector
        #read_ixs = []
        #for ix, read in enumerate(reads_fh):
        #    if ix % 4 == 0:
        #        read_acc = read.strip().replace('@gb', 'gb')
        #        read_ixs.append(read_acc)
        #        # initalise the gene family and aro similarity vectors
        #        # i.e. x_j for j is the similarity to the gene family
        #        # or aro of interest
        #        family_sim_vector = np.zeros(len(label_to_field))
        #        aro_sim_vector = np.zeros(len(aro_to_field))

        #        gene_family_encoding.update({read_acc : family_sim_vector})
        #        aro_encoding.update({read_acc: aro_sim_vector})


        # read the alignment file and store the top blast score per family
        # for each read
        if metric == 'bitscore':
            out_field = 11
        elif metric == 'evalue':
            out_field = 10
        elif metric == 'pident':
            out_field = 2
        else:
            raise ValueError('metric must be: {bitscore,evalue,pident}')

        scores = {}
        # without this we don't have encodings for anything with no DIAMOND
        # hits i.e. all 0 vectors.
        if self.data_type == 'aro':
            folder = 'training_data/subfamily_training_data'
        elif self.data_type == 'family':
            folder = 'training_data/family_training_data'

        encoding_fp = f'{folder}/{self.data_type}_X.txt'
        read_names_fp = f'{folder}/{self.data_type}_read_names.txt'

        if os.path.exists(encoding_fp) and os.path.exists(read_names_fp):
            print(f'Encoded {self.data_type} already exists so re-using'
                    ', use --redo to rebuild')
            alignment_fh.close()
            reads_fh.close()
        else:
            encoding_fh = open(encoding_fp, 'w')
            read_names_fh = open(read_names_fp, 'w')


            align_iter = itertools.groupby(alignment_fh,
                                         lambda x: x.split('\t')[0])


            for query_acc, query_hits in align_iter:
                # make new vector and assign new read to current
                vector = np.zeros(len(label_to_field))

                for hit in query_hits:
                    alignment = hit.strip().split('\t')
                    alignment_aro = alignment[1].split('|')[2]
                    score = float(alignment[out_field])

                    if self.data_type == 'family':
                        label = card.aro_to_gene_family[alignment_aro]
                    elif self.data_type == 'aro':
                        label = alignment_aro
                    else:
                        raise ValueError("data_type must be {family,aro}")

                    field = label_to_field[label]

                    # if this bitscore is greater than the already highest for that
                    # read and family
                    if metric == 'evalue':
                        score = math.log(score)
                        if score < vector[field]:
                            vector[field] = score
                    else:
                        if score > vector[field]:
                            vector[field] = score


                # moved onto the next read so we can dump the vector and move on
                encoded_vector = "\t".join([str(x) for x in np.nditer(vector)])
                encoding_fh.write(encoded_vector + "\n")
                read_names_fh.write(query_acc+ "\n")


            encoding_fh.close()
            read_names_fh.close()
            alignment_fh.close()

        x_encoding = np.loadtxt(encoding_fp, delimiter='\t')

        print(x_encoding.shape)

        # normalise bitscores
        if self.data_type == 'family':
            card.calculate_maximum_bitscores_per_family()
            max_bitscores = card.max_family_bitscores
        elif self.data_type == 'aro':
            card.calculate_maximum_bitscores_per_aro()
            max_bitscores = card.max_aro_bitscores

        if norm and metric == 'bitscore':
            print("Normalising encoding")
            x_encoding = normalize(x_encoding, axis=1, norm='l1')
            # numpy divide column by max per column
            #normalising = np.zeros(len(label_to_field))
            #for label, field in label_to_field.items():
            #    normalising[field] = card.max_family_bitscores[label]
            #    x_encoding = x_encoding / normalising


        elif norm and metric != 'bitscore':
            print("Can't normalise non bitscore metrics currently, must set metric to bitscore")
            sys.exit(1)

        # normalise and calculate dissimilarity matrices
        #if dissimilarity and norm:
        #    print("Calculating Dissimilarity")
        #    # convert this to new matrix
        #    family_df = family_df.applymap(lambda x: 1-x)
        #    family_df = family_df.fillna(0)
        #    aro_df = aro_df.applymap(lambda x: 1-x)
        #    aro_df = aro_df.fillna(0)
        #elif dissimilarity and not norm:
        #    print("Can't create dissimilarity matrix for non-normalised bitscores, must set dissimilarity and norm to True")
        #    sys.exit(1)

        return x_encoding


class Kmer():
    def __init__(self, metagenome_fp, k):
        self.metagenome_fp = metagenome_fp
        self.k = k

    def encode(self):
        tnf = ["".join(x) for x in itertools.product(['A', 'T', 'G', 'C', 'N'], repeat=self.k)]
        tnf_encoder = {v:k for k,v in enumerate(tnf)}

        X = []
        with open(self.metagenome_fp) as fh:
            for ix, line in enumerate(fh):
                if ix % 4 == 1:
                    clean_seq = re.sub("[M,X,R,S,Y,K]", 'N', line.strip())
                    self.seq = clean_seq
                    encoded_seq = self.read_encode(clean_seq, tnf_encoder)
                    X.append(encoded_seq)
        return np.vstack(X)

    def read_encode(self, seq, tnf):
        """
        k-mer decomposition might be simplest although might have to be
        done at the file level in order to represent all k-mers
        """
        encoded_vec = np.zeros(len(tnf))
        for tetranucleotide in self.window(self.k):
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

# dna2vec kmer embedding
class KMer_embedding():
    pass
