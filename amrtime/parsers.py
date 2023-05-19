#!/usr/bin/env python
import json
import math
import itertools
import subprocess
import numpy as np
import glob
import re
import os
from Bio import SeqIO
from Bio.Align import substitution_matrices


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
                                     'protein variant model']
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

            ## efflux components
            if len(gene_families) > 1 and any(["pump" in fam for fam in gene_families]):
                gene_families = ['efflux component']


          # also a fusion so assigned a new fusion class
            if 'class C LRA beta-lactamase' in gene_families and \
                    'class D LRA beta-lactamase' in gene_families:
                gene_families = ['class D/class C beta-lactamase fusion']

            # 23S with multiple resistance classes
            if ARO_acc == '3004181':
                gene_families = ['23S rRNA with mutation conferring resistance '
                                 'to macrolide and streptogramins antibiotics']

            if ARO_acc == '3007431':
                gene_families = ['alm glycyl carrier protein']

            if ARO_acc == '3007433':
                gene_families = ['alm glycyltransferase']

            if ARO_acc in ["3002547", "3005112", "3005114", "3005115",
                           "3005116", "3005117", "3005118", "3005119"]:
                gene_families = ["AAC(6')"]

            if ARO_acc == '3004364':
                gene_families = ['lipid A acyltransferase']

            if ARO_acc == "3007203":
                gene_families = ['pmr phosphoethanolamine transferase']

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
        """
        Get reverse dictionary of gene family to aros
        """
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

    def write_nucleotide_families(self, family_fasta_folder):
        """
        Write nucleotides sequences to per family fasta files
        """
        file_locations = {}
        for key, card_item in self.card.items():
            if card_item['model_type'] in self.supported_models:
                aro = card_item['ARO_accession']
                aro_name = card_item['ARO_name']
                sequences = card_item['model_sequences']['sequence']
                for seq_ix in sequences:
                    for sequence_type in sequences[seq_ix]:
                        if sequence_type == 'dna_sequence':
                            family = self.aro_to_gene_family[aro]
                            sequence = sequences[seq_ix][sequence_type]
                            fp = self.convert_amr_family_to_filename(family + '.fn',
                                    folder=family_fasta_folder)
                            with open(fp, 'a') as out_fh:
                                acc = ">gb|{}|{}|{}|".format(sequence['accession'],
                                                             aro,
                                                             aro_name.replace(' ', '_'))
                                seq = sequence['sequence']
                                out_fh.write(acc + '\n' + seq + '\n')
                            if family not in file_locations:
                                file_locations.update({family: fp})
        return file_locations

    def convert_amr_family_to_filename(self, family, folder=None):
        family = family.replace(' ', '_').replace('/', '@').replace("'", "^").replace('"', "+")
        if folder:
            fp = os.path.join(folder, family)
        else:
            fp = family
        return fp

    def convert_amr_family_filename_to_family(self, family):
        family = family.replace('_', ' ').replace('@', '/').replace("^", "'").replace("+", '"')
        return family

    def add_prevalence_to_family(self, prevalence_folder):
        for nt_fp in glob.glob(os.path.join(prevalence_folder, 'nucleotide_fasta_*')):
            for record in SeqIO.parse(nt_fp, "fasta"):
                try:
                    aro = record.description.split('|')[2].replace('ARO:', '')
                except:
                    print(record)
                    assert False
                with open(self.convert_amr_family_to_filename(\
                        self.aro_to_gene_family[aro], 'family_fasta'), 'a') as out_fh:
                    SeqIO.write(record, out_fh, "fasta")

    def calculate_maximum_bitscores_per_aro(self):
        """
        Determine the maximum bitscore per ARO and add a dictionary
        containing it to self
        """
        k_param = 0.711
        lambda_param = 1.37
        sub_matrix = substitution_matrices.load("BLOSUM62")

        aro_bitscores = {}
        for aro in self.proteins:
            raw_score = 0

            for aa in self.proteins[aro][1]:
                raw_score += sub_matrix[(aa, aa)]
            bitscore = (lambda_param * raw_score - math.log(k_param)) / math.log(2)

            aro_bitscores.update({aro: bitscore})

        self.max_aro_bitscores = aro_bitscores


    def calculate_maximum_bitscores_per_family(self):
        """
        Determine the maximum bitscore per AMR family and add a dictionary
        containing it to self
        """
        # k and lambda parameters are needed to normalise raw scores into
        # bit-scores but I can't find the damn things on the BLAST website
        k_param = 0.711
        lambda_param = 1.37
        sub_matrix = substitution_matrices.load("BLOSUM62")

        family_bitscores = {}
        for aro in self.proteins:
            raw_score = 0

            # to skip over the aro without proteins i.e. rRNA
            if len(self.proteins[aro][1]) == 0:
                continue

            for aa in self.proteins[aro][1]:
                raw_score += sub_matrix[(aa, aa)]
            bitscore = (lambda_param * raw_score - math.log(k_param)) / math.log(2)

            gene_family = self.aro_to_gene_family[aro]
            if gene_family not in family_bitscores:
                family_bitscores[gene_family] = [bitscore]
            else:
                family_bitscores[gene_family].append(bitscore)

        max_family_bitscores = {}
        for family, bitscores in family_bitscores.items():
            max_family_bitscores.update({family: max(bitscores)})

        self.max_family_bitscores = max_family_bitscores


    def generate_reads(self, aro, read_length, read_number, family=False):
        """
        Simulate reads using ART for the corresponding ARO,
        Optionally from the whole family

        Now how to handle balance when the diversity is variable
        is a real problem

        i.e. one classes covers more of sequence space than the other.
        """
        if not family:
            seq = {aro: self.nucleotides[aro]}
            tmp_fasta_file = f'training_data/{aro}.fasta'
        if family:
            family = self.aro_to_gene_family[aro]
            family_aros = self.gene_family_to_aro[family]
            seq = {aro: self.nucleotides[aro] for aro in family_aros}
            fp = self.convert_amr_family_to_filename(family)
            tmp_fasta_file = f'training_data/{fp}.fasta'

        self.write_seqs(seq, 'training_data/tmp.fasta')

        # -c is read count per sequence and -f is coverage
        subprocess.check_call(f'art_illumina -q -na -ss MSv3 -i \
                    {tmp_fasta_file} -c {read_number} -l {read_length} \
                    -rs 42 -o training_data/{aro}', shell=True)


    #def find_family_overlaps(self, read_length):
    #    """
    #    Find perfect overlaps of length read_length between different
    #    members of a AMR family

    #    i.e. how many 250bp overlaps are there between all members of NDM
    #    in protein space

    #    Why use protein space? Because that's what we are querying and it
    #    means shorter sequences to compare quickly
    #    """

    #    for family in
    #    if not family:
    #        seq = {aro: self.nucleotides[aro]}
    #        tmp_fasta_file = f'training_data/{aro}.fasta'
    #    if family:
    #        family = self.aro_to_gene_family[aro]
    #        family_aros = self.gene_family_to_aro[family]
    #        seq = {aro: self.nucleotides[aro] for aro in family_aros}



def prepare_labels(fp, card):
    """
    Parse label file and return list of ARO labels and their higher level
    gene family level labels (using the CARD class lookup
    """
    families = []
    aros = []
    with open(fp) as fh:
        for line in fh:
            aro = line.split()[2]
            family = card.aro_to_gene_family[aro]
            families.append(family)
            aros.append(aro)

    return families, aros
