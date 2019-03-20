#!/usr/bin/env python
import json
import math
import itertools
import numpy as np
import glob
import re
import os
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62

class CARD():
    """
    Parser for the CARD database
    """
    def __init__(self, card_json_fp, rrna=False):
        with open(card_json_fp) as fh:
            with open(card_json_fp) as fh:
                self.card = json.load(fh)
            self.version = self.card['_version']

            # to avoid having to except them later when parsing the other
            # entries
            del self.card['_version']
            del self.card['_timestamp']
            del self.card['_comment']

            """
            All models currently in card.json:
                {'efflux pump system meta-model',
                'gene cluster meta-model',
                'protein domain meta-model',
                'protein homolog model',
                'protein knockout model',
                'protein overexpression model',
                'protein variant model',
                'rRNA gene variant model'}

            rRNA is being handled by metaRNA approach

            Other metamodels can be added later
            """


            self.supported_models = ['protein homolog model',
                                     'protein variant model',
                                     'protein overexpression model']
            if rrna:
                self.supported_models.append('rRNA gene variant model')

            self.proteins, self.nucleotides = self.get_sequences()
            self.aro_to_gene_family = self.build_aro_to_gene_family()
            self.gene_family_to_aro = self.build_gene_family_to_aro()

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
        aros_without_protein = []
        for aro, gene_families in aro_to_gene_family.items():
            if len(gene_families) != 1:
                mapping_failure = True
                print(aro, gene_families)
            else:
                aro_to_gene_family[aro] = gene_families[0]

            if aro not in self.proteins:
                aros_without_protein.append(aro)

        if mapping_failure:
            raise ValueError("AROs and gene families don't map 1:1")

        # remove any aro without a protein sequence
        for aro in aros_without_protein:
            del aro_to_gene_family[aro]

        return aro_to_gene_family

    def build_gene_family_to_aro(self):
        gene_family_to_aro = {}
        for aro, gene_family in self.aro_to_gene_family.items():
            if gene_family not in gene_family_to_aro:
                gene_family_to_aro.update({gene_family: [aro]})
            else:
                gene_family_to_aro[gene_family].append(aro)

        return gene_family_to_aro

    def get_sequences(self):
        """
        Gather list of accession, prot and nucleotide sequence tuples from
        the card.json
        """
        data = {'protein_sequence': {}, 'dna_sequence': {}}
        for key, card_item in self.card.items():
            if card_item['model_type'] in self.supported_models:
                aro = card_item['ARO_accession']
                aro_name = card_item['ARO_name']
                sequences = card_item['model_sequences']['sequence']
                for seq_ix in sequences:
                    for sequence_type in sequences[seq_ix]:
                        if sequence_type in ['protein_sequence', 'dna_sequence']:
                            sequence = sequences[seq_ix][sequence_type]
                            if aro not in data[sequence_type]:
                                acc = ">gb|{}|{}|{}|".format(sequence['accession'],
                                                             aro,
                                                             aro_name.replace(' ', '_'))
                                data[sequence_type].update({aro: (acc,
                                                                  sequence['sequence'])})
        # what the fuck is happening why are we getting identical sequences
        # temporarily let's just add the first sequence
                            #else:
                            #    data[sequence_type][aro].append((acc,
                            #                                     sequence['sequence']))

        return data['protein_sequence'], data['dna_sequence']

    def write_seqs(self, seq_dict, seq_file_fp):
        if not os.path.exists(seq_file_fp):
            with open(seq_file_fp, 'w') as fh:
                for seq in seq_dict.values():
                    fh.write("{}\n{}\n".format(seq[0], seq[1]))

    def write_proteins(self, seq_file_fp):
        self.write_seqs(self.proteins, seq_file_fp)

    def write_nucleoties(self, seq_file_fp):
        self.write_seqs(self.nucleotides, seq_file_fp)

    def get_nucleotide_per_family(self):
        """
        Write nucleotides sequences to per family fasta files
        """
        if not os.path.exists('family_fasta'):
            os.mkdir('family_fasta')

        for key, card_item in self.card.items():
            if card_item['model_type'] in self.supported_models:
                aro = card_item['ARO_accession']
                aro_name = card_item['ARO_name']
                sequences = card_item['model_sequences']['sequence']
                for seq_ix in sequences:
                    for sequence_type in sequences[seq_ix]:
                        if sequence_type == 'dna_sequence':
                            sequence = sequences[seq_ix][sequence_type]
                            with open(self.convert_amr_family_to_filename(self.aro_to_gene_family[aro]), 'a') as out_fh:
                                acc = ">gb|{}|{}|{}|".format(sequence['accession'],
                                                             aro,
                                                             aro_name.replace(' ', '_'))
                                seq = sequence['sequence']
                                out_fh.write(acc + '\n' + seq + '\n')

    def convert_amr_family_to_filename(self, family):
        fp = os.path.join('family_fasta', family.replace(' ', '_').replace('/', '_'))
        return fp

    def add_prevalence_to_family(self, prevalence_folder):
        for nt_fp in glob.glob(os.path.join(prevalence_folder, 'nucleotide_fasta_*')):
            for record in SeqIO.parse(nt_fp, "fasta"):
                try:
                    aro = record.description.split('|')[2].replace('ARO:', '')
                except:
                    print(record)
                    assert False
                with open(self.convert_amr_family_to_filename(self.aro_to_gene_family[aro]), 'a') as out_fh:
                    SeqIO.write(record, out_fh, "fasta")


if __name__ == '__main__':
    card = CARD('training_data/card-data/card.json')
    card.get_nucleotide_per_family()
    card.add_prevalence_to_family('training_data/card-data/prevalence')

