# -*- coding: utf-8 -*-

import argparse
import biovec


parser = argparse.ArgumentParser(description='This script is ...')
parser.add_argument('-i', '--input_fasta_file')
parser.add_argument('-o', '--output_file')
args = parser.parse_args()

model = biovec.ProtVec(corpus_fname=args.input_fasta_file)
model.save(args.output_file)
