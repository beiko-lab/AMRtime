
if __name__=="__main__":

    name_to_aro = {'mel': "3000616", 'pmrA':"3000822" , 'patB': "3000025"}
    strict_perfect = {'500': "strict", '1000': "perfect"}

    expected_bedtool_labels = open('bedtools_labels.tsv', 'w')

    labels = {}
    with open('intersection.tsv') as fh:
        for line in fh:
            line = line.strip().split('\t')
            labels.update({line[3]: [line[3], line[0], name_to_aro[line[-4]],
                            line[-4], strict_perfect[line[-3]], line[-1]]})

    with open('output.fq') as fh:
        for line_num, line in enumerate(fh):
            if line_num % 4 == 0:
                line = line.strip()
                line = line.replace("@", '')
                if line in labels:
                    expected_bedtool_labels.write("\t".join(labels[line])+"\n")
                else:
                    expected_bedtool_labels.write("\t".join(["na"]*6)+"\n")




