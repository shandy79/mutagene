import argparse
import sys

from mutagene.mutability.mutability import THRESHOLD_DRIVER, THRESHOLD_PASSENGER
import mutagene.core.rank as rank

class RankMenu(object):
    def __init__(self, parser):
        parser.add_argument("--infile", "-i", nargs='*', help="Input file in MAF format", type=argparse.FileType('r'))
        parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
        parser.add_argument('--genome', "-g", help="Location of genome assembly file", type=str)

        # parser.add_argument('--mode', "-m", help="n (nucleotide) or aa (amino acid) mode", type=str, default="aa")

        parser.add_argument('--cohorts-file', type=str, help="Location of tar.gz container or directory for cohorts", default="cohorts.tar.gz")
        parser.add_argument('--cohort', "-c", type=str, help="Name of cohort with observed mutations")

        parser.add_argument('--profile', "-p", help="Override profile to calculate mutability, may also describe cohort size", type=str)
        parser.add_argument('--nsamples', "-n", type=int, help="Override cohort size")
        parser.add_argument('--threshold-driver', "-td", help="BScore threshold between Driver and Pontential Driver mutations", type=float, default=THRESHOLD_DRIVER)
        parser.add_argument('--threshold-passenger', "-tp", help="BScore threshold between Pontential Driver and Passenger mutations", type=float, default=THRESHOLD_PASSENGER)

    def callback(self, args):
        rank.rank_driver(input_files=args.infile, genome=args.genome, cohort=args.cohort, output_file=args.outfile, cohorts_file=args.cohorts_file,
                profile_file=args.profile, nsamples=args.nsamples, threshold_driver=args.threshold_driver, threshold_passenger=args.threshold_passenger)
