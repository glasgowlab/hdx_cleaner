import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Rangeslist merger')
parser.add_argument('-r', dest='r', help="input", required=True)
parser.add_argument('-m', dest='m', help="mask", required=True)
parser.add_argument('-o', dest='o', help="output", required=True)
parser.add_argument('-e', '-exclude', dest='exclude', action='store_true', help='exclude rather than include rangeslist')

args = parser.parse_args()

r = pd.read_csv(args.r)
mask = pd.read_csv(args.m).drop_duplicates()

if args.exclude:
    exc = pd.merge(mask, r, how='left', indicator=True)
    exc[exc['_merge'] == 'left_only'].drop('_merge', axis=1).to_csv(args.o, index=False)
else:
    pd.merge(r,mask).to_csv(args.o, index=False)
