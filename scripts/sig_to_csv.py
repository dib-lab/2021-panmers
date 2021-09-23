#!/usr/bin/env python3
"""
Given a sourmash signature, extract the minhash integers and abund info
to a csv file with the basename of the sig file name as the name of the 
abund column.
"""

from sourmash import signature
import pandas as pd
import os
import sys
import argparse

def main():
    p = argparse.ArgumentParser()
    p.add_argument('signature')       # sourmash signature
    p.add_argument('output')          # output csv file name
    args = p.parse_args()

    # load the signature from disk
    sigfp = open(args.signature, 'rt')
    siglist = list(signature.load_signatures(sigfp))
    loaded_sig = siglist[0]

    # Get the minhashes
    mins = loaded_sig.minhash.hashes.keys()
    
    name = os.path.basename(args.signature)
    df = pd.DataFrame(mins, columns=[name])
    
    # write to a csv
    df.to_csv(args.output, index = False)

if __name__ == '__main__':
    sys.exit(main())
