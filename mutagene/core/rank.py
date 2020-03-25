import sys
import logging
from pathlib import Path

from mutagene.mutability.mutability import rank, THRESHOLD_DRIVER, THRESHOLD_PASSENGER
from mutagene.io.cohorts import read_cohort_mutations_from_tar
from mutagene.io.cohorts import read_cohort_size_from_profile_file, list_cohorts_in_tar
from mutagene.io.fetch import GENOME_ERROR_MESSAGE
from mutagene.io.profile import read_profile_file
from mutagene.io.protein_mutations_MAF import read_MAF_with_genomic_context


logger = logging.getLogger(__name__)


def rank_driver(genome, input_files, cohort, output_file=sys.stdout, cohorts_file='cohorts.tar.gz', profile_file='', nsamples=0,
        threshold_driver=THRESHOLD_DRIVER, threshold_passenger=THRESHOLD_PASSENGER):

    if not genome:
        logger.warning(GENOME_ERROR_MESSAGE)
        return

    if not Path(cohorts_file).is_file():
        logger.warning("Requires a cohorts file (e.g., 'cohorts.tar.gz').  Either specify with the 'cohorts_file' argument or use 'mutagene fetch cohorts' to download.")
        return

    if cohort and Path(cohorts_file).is_file():
        profile, cohort_size, cohort_aa_mutations, cohort_na_mutations = read_cohort_mutations_from_tar(cohorts_file, cohort)
        logger.info('Profile and cohort size loaded from precalculated cohorts N=' + str(cohort_size))
    else:
        logger.warning('Cohort required')

        if cohorts_file:
            logger.warning("List of available cohorts:\n" + list_cohorts_in_tar(cohorts_file))
        return

    if profile_file:
        profile = read_profile_file(profile_file)
        if profile:
            logger.info('Profile overridden')
        else:
            return
        cohort_size_new = read_cohort_size_from_profile_file(profile_file)
        if cohort_size_new:
            cohort_size = cohort_size_new
            logger.info('Cohort size loaded from profile N=' + str(cohort_size))

    if nsamples:
        cohort_size = nsamples
        logger.info('Cohort size overridden N=' + str(cohort_size))

    if isinstance(input_files, list):
        if len(input_files) > 1:
            logger.info('Multiple input files provided')
        mutations_to_rank = {}
        processing_stats = {'loaded': 0, 'skipped': 0}

        for infile in input_files:
            if isinstance(infile, (str, bytes, Path)):
                with open(infile, 'r') as f:
                    mut, stats = read_MAF_with_genomic_context(f, genome)
            else:
                mut, stats = read_MAF_with_genomic_context(infile, genome)

            mutations_to_rank.update(mut)
            processing_stats['loaded'] += stats['loaded']
            processing_stats['skipped'] += stats['skipped']
    else:
        if isinstance(input_files, (str, bytes, Path)):
            with open(input_files, 'r') as f:
                mutations_to_rank, processing_stats = read_MAF_with_genomic_context(f, genome)
        else:
            mutations_to_rank, processing_stats = read_MAF_with_genomic_context(input_files, genome)

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

    logger.info("THRESHOLD_DRIVER: {}".format(threshold_driver))
    logger.info("THRESHOLD_PASSENGER: {}".format(threshold_passenger))

    rank(mutations_to_rank, output_file, profile, cohort_aa_mutations, cohort_size, threshold_driver, threshold_passenger)
