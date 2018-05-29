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
            if "glycopeptide resistance gene cluster" in gene_families:
                gene_families.remove('glycopeptide resistance gene cluster')
            if "fluoroquinolone self resistant parC" in gene_families:
                gene_families.remove('fluoroquinolone self resistant parC')
            if "AAC(6')" in gene_families:
                gene_families = ["AAC(6')"]
            if 'ATP-binding cassette (ABC) antibiotic efflux pump' in gene_families and 'pmr phosphoethanolamine transferase' in gene_families:
                gene_families = ['pmr phosphoethanolamine transferase']
            if 'resistance-nodulation-cell division (RND) antibiotic efflux pump' in gene_families and 'General Bacterial Porin with reduced permeability to beta-lactams' in gene_families:
                gene_families = ['resistance-nodulation-cell division (RND) antibiotic efflux pump']
            if "kirromycin self resistant EF-Tu" in gene_families:
                gene_families.remove('kirromycin self resistant EF-Tu')
            if 'aminocoumarin self resistant parY' in gene_families:
                gene_families.remove('aminocoumarin self resistant parY')
            if 'major facilitator superfamily (MFS) antibiotic efflux pump' in gene_families and 'resistance-nodulation-cell division (RND) antibiotic efflux pump' in gene_families:
                gene_families = ['efflux pump']
            if 'major facilitator superfamily (MFS) antibiotic efflux pump' in gene_families and 'ATP-binding cassette (ABC) antibiotic efflux pump' in gene_families:
                gene_families = ['efflux pump']
            if 'class C LRA beta-lactamase' in gene_families and 'class D LRA beta-lactamase' in gene_families:
                gene_families = ['class D/class C beta-lactamase fusion']
            if 'fluoroquinolone resistant parC' in gene_families and 'fluoroquinolone resistant gyrA' in gene_families:
                gene_families = ['fluoroquinolone resistant parC']
            if 'UhpT' in gene_families and 'UhpA' in gene_families:
                gene_families = ['UhpA']

            if ARO_acc == "3004450":
                gene_families = ['TRU beta-lactamase']
            if ARO_acc == "3004294":
                gene_families = ['BUT beta-lactamase']



            aro_to_gene_family.update({ARO_acc: gene_families})

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
        """
        for i in j:
	    # model_type: protein homolog model
	    if j[i]['model_type_id'] == '40292':
                pass_bit_score = j[i]['model_param']['blastp_bit_score']['param_value']

                for seq in j[i]['model_sequences']['sequence']:
                    fout.write('>%s_%s | model_type_id: 40292 | pass_bitscore: %s | %s\n' % (i, seq, pass_bit_score, j[i]['ARO_name']))
	            fout.write('%s\n' %(j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))

	        	# model_type: protein variant model
            elif j[i]["model_type_id"] == "40293":
                pass_bit_score = j[i]['model_param']['blastp_bit_score']['param_value']
                snpList = [j[i]['model_param']['snp']['param_value'][k] for k in j[i]['model_param']['snp']['param_value']]
                for seq in j[i]['model_sequences']['sequence']:
	            fout.write('>%s_%s | model_type_id: 40293 | pass_bit_score: %s | SNP: %s | %s\n' \
	                            % (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
                    fout.write('%s\n' % (j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))

	        	# model_type: protein overexpression model
            elif j[i]["model_type_id"] == "41091":
                pass_bit_score = j[i]["model_param"]["blastp_bit_score"]["param_value"]
                snpList = [j[i]['model_param']['snp']['param_value'][k] for k in j[i]['model_param']['snp']['param_value']]
                for seq in j[i]['model_sequences']['sequence']:
	            fout.write('>%s_%s | model_type_id: 41091 | pass_bit_score: %s | SNP: %s | %s\n' \
	                        % (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
                    fout.write('%s\n' % (j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))


            elif j[i]["model_type_id"] == "40295":
                pass_bit_score = j[i]['model_param']['blastn_bit_score']['param_value']
                snpList = [j[i]['model_param']['snp']['param_value'][k] for k in j[i]['model_param']['snp']['param_value']]
	        for s in snpList:
	            if "16S" in j[i]['ARO_name']:
	                if s not in snpList_16s:
                            snpList_16s.append(s)
                    if "23S" in j[i]['ARO_name']:
	                if s not in snpList_23s:
                            snpList_23s.append(s)

                    for seq in j[i]['model_sequences']['sequence']:
	                if j[i]['model_sequences']['sequence'][seq]['dna_sequence']['strand'] == "-":
                            basecomplement = self.complementary_strand(j[i]['model_sequences']['sequence'][seq]['dna_sequence']['sequence'])
                            fout.write('>%s_%s | model_type_id: 40295 | pass_bit_score: %s | SNP: %s | %s\n' \
	    	                    % (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
                            fout.write('%s\n' % (basecomplement))

	    	    else:
	    	        fout.write('>%s_%s | model_type_id: 40295 | pass_bit_score: %s | SNP: %s | %s\n' \
	    	                    % (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
                            fout.write('%s\n' % (j[i]['model_sequences']['sequence'][seq]['dna_sequence']['sequence']))
        """

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

