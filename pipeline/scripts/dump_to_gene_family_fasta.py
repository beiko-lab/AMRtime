from amrtime import parsers

if __name__ == '__main__':
    card = parsers.CARD('training_data/card-data/card.json')
    card.get_nucleotide_per_family()
    card.add_prevalence_to_family('training_data/card-data/prevalence')

