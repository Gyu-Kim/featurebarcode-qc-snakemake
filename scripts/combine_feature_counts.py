

import os
import sys


def parse_sample_files(sample_files, names):
    sample_file_dict = {}

    sample_files = sample_files.split(',')
    if names is None:
        names = [ os.path.splitext(os.path.basename(sample_file))[0]
                                for sample_file in sample_files ]
    else:
        names = names.split(',')

    for sample_file, name in zip(sample_files, names):
        if not os.path.isfile(sample_file):
            sys.exit("Error: sample file " + sample_file + " not found")
        sample_file_dict[name] = sample_file

    return sample_file_dict



def collapse_feature_counts(sample_files, out_file, names):


    sample_file_dict = parse_sample_files(sample_files, names)

    o = open(out_file, 'w')
    print("sample", "feature", "cell", "count", sep="\t", file=o)

    for sample,count_file in sample_file_dict.items():
        i = open(count_file, 'r')
        header = i.readline()
        for line in i:
            o.write(sample + '\t' + line)
        i.close()

    o.close()


if __name__ == "__main__":

##############################################################################
    # Initiates argument parser
    import argparse

    parser = argparse.ArgumentParser(description='Combine multiple feature barcode count files')

    parser.add_argument('sample_files', help='Comma-separated list of sample count files', type=str)

    parser.add_argument('out_file', help = 'Where to write the CSV output file', type=str)

    # Optional arguments

    parser.add_argument('--names', help='Comma-separated list of sample names',
                        default=None, type=str)

    args = parser.parse_args()


##############################################################################


    collapse_feature_counts(args.sample_files, args.out_file, args.names)






