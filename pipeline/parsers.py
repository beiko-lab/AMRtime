#!/usr/bin/env python
import json


class CARD():
    def __init__(self, card_json_fp):
        with open(card_json_fp) as fh:
            self.card = json.load(card_json_fp)
            self.version = self.card['_version']

            # to avoid having to except them later when parsing the other
            # entries
            del self.card['_version']
            del self.card['_timestamp']
            del self.card['_comment']

            self.aro_to_gene_family = build_aro_to_gene_family()
            self.gene_family_to_aro = build_gene_family_to_aro()

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
            aro_to_gene_family.update({ARO_acc: gene_family})
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
        pass


class Read():
    """
    Class to keep the reads encoding nice and tidy
    """
    def __init__(self, read_seq):
        self.seq = read_seq
        self.encoded_seq = self.encode()

    def encode(self):
        """
        k-mer decomposition might be simplest although might have to be
        done at the file level in order to represent all k-mers
        """
        pass

