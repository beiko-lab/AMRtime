#!/usr/bin/env python
import json

import itertools
import numpy as np
import re

class CARD():
    def __init__(self, card_json_fp):
        with open(card_json_fp) as fh:
            with open(card_json_fp) as fh:
                self.card = json.load(fh)
            self.version = self.card['_version']

            # to avoid having to except them later when parsing the other
            # entries
            del self.card['_version']
            del self.card['_timestamp']
            del self.card['_comment']

            self.aro_to_gene_family = self.build_aro_to_gene_family()
            self.gene_family_to_aro = self.build_gene_family_to_aro()

            self.accessions, self.sequences = self.get_protein_sequences()

    def build_aro_to_gene_family(self):
        aro_to_gene_family = {}
        for card_item in self.card.values():
            ARO_acc = card_item['ARO_accession']

            # as multiple gene families are possible per ARO
            # although they are relatively rare so maybe I should talk with
            # Andrew about making them unique
            gene_families = []
            for category in card_item['ARO_category'].values():
                if category['category_aro_class_name'] == 'AMR Gene Family':
                    gene_families.append(category['category_aro_name'])

            # fix the situations where there are multiple gene families manually
            # all glycopeptide resistant gene clusters have 2 gene families one
            # indicating that it is a grgc and the other with its class
            # for now we are just using the cluster level and will deal with
            # the specifics at ARO level.
            if "glycopeptide resistance gene cluster" in gene_families:
                gene_families = ['glycopeptide resistance gene cluster']

            # this is a fusion protein so can be assigned to a new class
            if ARO_acc == '3002598':
                gene_families = ["AAC(6')_ANT(3'')"]
            if ARO_acc == '3002597':
                gene_families = ["APH(2'')_AAC(6')"]
            if ARO_acc in ['3002546', '3002600']:
                gene_families = ["AAC(3)_AAC(6')"]

            # also a fusion so assigned a new fusion class
            if 'class C LRA beta-lactamase' in gene_families and \
                    'class D LRA beta-lactamase' in gene_families:
                gene_families = ['class D/class C beta-lactamase fusion']

            # 23S with multiple resistance classes
            if ARO_acc == '3004181':
                gene_families = ['23S rRNA with mutation conferring resistance '
                                 'to macrolide and streptogramins antibiotics']

            # additional self-resistance class to indicate resistance genes made by
            # antibiotic producer removing the self resistant term
            if "fluoroquinolone resistant parC" in gene_families and \
                    "fluoroquinolone self resistant parC" in gene_families:
                gene_families = ['fluoroquinolone resistant parC']

            if "kirromycin self resistant EF-Tu" in gene_families and \
                    'elfamycin resistant EF-Tu' in gene_families:
                gene_families = ['elfamycin resistant EF-Tu']

            if 'aminocoumarin self resistant parY' in gene_families and \
                    'aminocoumarin resistant parY' in gene_families:
                gene_families = ['aminocoumarin resistant parY']

            # efflux components
            if ARO_acc in ['3000263', '3000833', '3003382', '3000832',
                           '3000815', '3003896', '3000823', '3003511',
                           '3003381', '3000817', '3003895', '3000676',
                           '3003383', '3003585', '3004107', '3003820']:
                gene_families = ['efflux regulator']
            if ARO_acc == '3000237':
                gene_families = ['efflux component']


            # They are homologous parts of topo IV and II but it looks like this is actually parC
            if 'fluoroquinolone resistant parC' in gene_families and \
                    'fluoroquinolone resistant gyrA' in gene_families:
                gene_families = ['fluoroquinolone resistant parC']


            # things that need fixed
            # this looks like a mistake and is only UhpA
            if 'UhpT' in gene_families and 'UhpA' in gene_families:
                gene_families = ['UhpA']

            # missing families
            if ARO_acc == "3004450":
                gene_families = ['TRU beta-lactamase']
            if ARO_acc == "3004294":
                gene_families = ['BUT beta-lactamase']

            aro_to_gene_family.update({ARO_acc: gene_families})

        # tidy up and make unique aro:amr family relationship
        mapping_failure = False
        for aro, gene_families in aro_to_gene_family.items():
            if len(gene_families) != 1:
                mapping_failure = True
                print(aro, gene_families)
            else:
                aro_to_gene_family[aro] = gene_families[0]

        if mapping_failure:
            raise ValueError("AROs and gene families don't map 1:1")

        return aro_to_gene_family

    def build_gene_family_to_aro(self):
        gene_family_to_aro = {}
        for key, value in self.aro_to_gene_family.items():
            for gene_family in value:
                if gene_family not in gene_family_to_aro:
                    gene_family_to_aro.update({gene_family: [key]})
                else:
                    gene_family_to_aro[gene_family].append(key)
        return gene_family_to_aro

    def get_protein_sequences(self):
        """
        Gather list of accession, sequence tuples from the card.json
        """
        # from rgi

        accessions = []
        sequences = []

        for card_item in self.card.values():
	    # model_type: protein homolog model
            if card_item['model_type_id'] in ['40292', '40293', '41091']:

                for sequence_data in card_item['model_sequences']['sequence'].values():
                    accessions.append("ARO:{}|{}".format(card_item['ARO_id'],
                                                         card_item['ARO_name']))
                    seq = sequence_data['dna_sequence']['sequence']
                    #if sequence_data['dna_sequence']['strand'] == '-':
                    #    seq = self.base_complement(seq)
                    sequences.append(seqs)

        return accessions, sequences


def read_metagenome(fp):
    """
    Prepare metagenome
    """
    tnf = ["".join(x) for x in itertools.product(['A', 'T', 'G', 'C', 'N'], repeat=4)]
    tnf_encoder = {v:k for k,v in enumerate(tnf)}

    X = []
    with open(fp) as fh:
        for ix, line in enumerate(fh):
            if ix % 4 == 1:
                seq = Read(line.strip())
                encoded_seq = seq.encode(tnf_encoder)
                X.append(encoded_seq)
    return np.vstack(X)


def prepare_labels(fp):
    aros = []
    with open(fp) as fh:
        for line in fh:
            aro = line.split()[2]
            aros.append(aro)
    return aros


class Read():
    """
    Class to keep the reads encoding nice and tidy
    """
    def __init__(self, read_seq):
        self.seq = re.sub("[M,X,R,S,Y,K]", 'N', read_seq)

    def encode(self, tnf):
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

