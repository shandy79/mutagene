import sys
import logging
from pathlib import Path

# https://www.python.org/dev/peps/pep-0443/
from functools import singledispatch
import argparse, _io, yaml

from mutagene.mutability.mutability import rank, THRESHOLD_DRIVER, THRESHOLD_PASSENGER
from mutagene.io.cohorts import read_cohort_mutations_from_tar
from mutagene.io.cohorts import read_cohort_size_from_profile_file, list_cohorts_in_tar
from mutagene.io.fetch import GENOME_ERROR_MESSAGE
from mutagene.io.profile import read_profile_file
from mutagene.io.protein_mutations_MAF import read_MAF_with_genomic_context
from mutagene.core.core_utils import inputs_to_yml


logger = logging.getLogger(__name__)


_DEFAULTS_DICT = {
    'command': 'rank',
    'genome': None,
    'input_files': None,
    'cohort': None,
    'output_file': sys.stdout,
    'cohorts_file': 'cohorts.tar.gz',
    'profile_file': None,
    'nsamples': None,
    'threshold_driver': THRESHOLD_DRIVER,
    'threshold_passenger': THRESHOLD_PASSENGER
}


def create_input_dict(genome, input_files, cohort, output_file, cohorts_file, profile_file, nsamples, threshold_driver, threshold_passenger):
    input_dict = _DEFAULTS_DICT.copy()

    # Merge inputs w/defaults
    if genome is not None:  input_dict['genome'] = genome
    if input_files is not None:  input_dict['input_files'] = input_files  # and len(input_files) > 0
    if cohort is not None:  input_dict['cohort'] = cohort
    if output_file is not None:  input_dict['output_file'] = output_file
    if cohorts_file is not None:  input_dict['cohorts_file'] = cohorts_file
    if profile_file is not None:  input_dict['profile_file'] = profile_file
    if nsamples is not None:  input_dict['nsamples'] = nsamples  # and nsamples >= 0
    if threshold_driver is not None:  input_dict['threshold_driver'] = threshold_driver  # and threshold_driver > 0
    if threshold_passenger is not None:  input_dict['threshold_passenger'] = threshold_passenger  # and threshold_driver > 0

    return input_dict


@singledispatch
def rank_adapter(genome, input_files, cohort, output_file=sys.stdout, cohorts_file='cohorts.tar.gz', profile_file=None, nsamples=None,
        threshold_driver=THRESHOLD_DRIVER, threshold_passenger=THRESHOLD_PASSENGER):
    input_dict = create_input_dict(genome, input_files, cohort, output_file, cohorts_file, profile_file, nsamples, threshold_driver, threshold_driver)
    rank_driver(input_dict)


@rank_adapter.register(argparse.Namespace)
def rank_adapter_argparse(argparse_dict):
#    print(argparse_dict)
    input_dict = create_input_dict(argparse_dict.genome, argparse_dict.infile, argparse_dict.cohort, argparse_dict.outfile,
            argparse_dict.cohorts_file, argparse_dict.profile, argparse_dict.nsamples, argparse_dict.threshold_driver, argparse_dict.threshold_passenger)
    rank_driver(input_dict)


@rank_adapter.register(_io.TextIOWrapper)
def rank_adapter_yml(yml_file):
    yml_data = yaml.safe_load(yml_file)
    print(yml_data)
    input_dict = create_input_dict(yml_data['genome'] if 'genome' in yml_data else None,
            yml_data['input_files'] if 'input_files' in yml_data else None,
            yml_data['cohort'] if 'cohort' in yml_data else None,
            yml_data['output_file'] if 'output_file' in yml_data else None,
            yml_data['cohorts_file'] if 'cohorts_file' in yml_data else None,
            yml_data['profile_file'] if 'profile_file' in yml_data else None,
            yml_data['nsamples'] if 'nsamples' in yml_data else None,
            yml_data['threshold_driver'] if 'threshold_driver' in yml_data else None,
            yml_data['threshold_passenger'] if 'threshold_passenger' in yml_data else None)
    rank_driver(input_dict)


def rank_driver(input_dict):
    # Perform validation, throw exception or output warning if invalid
    # Question of how many arguments should be allowed to have defaults

    if not input_dict['genome']:
        logger.warning(GENOME_ERROR_MESSAGE)
        return

    if not Path(input_dict['cohorts_file']).is_file():
        logger.warning("Cohorts file missing. Download with \"mutagene fetch_cohorts\"")
        return

    if input_dict['cohort'] and Path(input_dict['cohorts_file']).is_file():
        profile, cohort_size, cohort_aa_mutations, cohort_na_mutations = read_cohort_mutations_from_tar(input_dict['cohorts_file'], input_dict['cohort'])
        logger.info('Profile and cohort size loaded from precalculated cohorts N=' + str(cohort_size))
    else:
        logger.warning('Cohort required')

        if input_dict['cohorts_file']:
            logger.warning("List of available cohorts:\n" + list_cohorts_in_tar(input_dict['cohorts_file']))
        return

    if input_dict['profile_file']:
        profile = read_profile_file(input_dict['profile_file'])
        if profile:
            logger.info('Profile overridden')
        else:
            return
        cohort_size_new = read_cohort_size_from_profile_file(input_dict['profile_file'])
        if cohort_size_new:
            cohort_size = cohort_size_new
            logger.info('Cohort size loaded from profile N=' + str(cohort_size))

    if input_dict['nsamples']:
        cohort_size = input_dict['nsamples']
        logger.info('Cohort size overridden N=' + str(cohort_size))

    if isinstance(input_dict['input_files'], list):
        if len(input_dict['input_files']) > 1:
            logger.info('Multiple input files provided')
        mutations_to_rank = {}
        processing_stats = {'loaded': 0, 'skipped': 0}

        # Updated to accept either a file handle or a file name as a string to better support the API
        for infile in input_dict['input_files']:
            if isinstance(infile, (str, bytes, Path)):
                with open(infile, 'r') as f:
                    mut, stats = read_MAF_with_genomic_context(f, input_dict['genome'])
            else:
                mut, stats = read_MAF_with_genomic_context(infile, input_dict['genome'])

            mutations_to_rank.update(mut)
            processing_stats['loaded'] += stats['loaded']
            processing_stats['skipped'] += stats['skipped']
    else:
        if isinstance(input_dict['input_files'], (str, bytes, Path)):
            with open(input_dict['input_files'], 'r') as f:
                mutations_to_rank, processing_stats = read_MAF_with_genomic_context(f, input_dict['genome'])
        else:
            mutations_to_rank, processing_stats = read_MAF_with_genomic_context(input_dict['input_files'], input_dict['genome'])

    if not len(mutations_to_rank):
        logger.warning('MutaGene rank failed: No mutations to rank. Check that the infile is in MAF format')
        return

    msg = "Loaded {} mutations".format(processing_stats['loaded'])
    if processing_stats['skipped'] > 0:
        msg += " skipped {} mutations".format(processing_stats['skipped'])
    logger.info(msg)

    # mutations_to_rank = list(mutations_to_rank)

    # not_unique = len(mutations_to_rank)
    # mutations_to_rank = list(set(mutations_to_rank))
    # unique = len(mutations_to_rank)
    # if not_unique != unique:
    #     logger.info("Removed {} duplicate mutations".format(not_unique - unique))

    logger.info("THRESHOLD_DRIVER: {}".format(input_dict['threshold_driver']))
    logger.info("THRESHOLD_PASSENGER: {}".format(input_dict['threshold_passenger']))

    # Output run inputs and other information to YML file for reproducibility
    inputs_to_yml(input_dict)

    rank(mutations_to_rank, input_dict['output_file'], profile, cohort_aa_mutations, cohort_size,
            input_dict['threshold_driver'], input_dict['threshold_passenger'])
