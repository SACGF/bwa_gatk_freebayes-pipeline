#!/usr/bin/env python

"""

"""

import pysam

def get_sample_list(bam):
    """
    input is a BAM file
    returns a list of samples defined in the Header section of the input BAM
    """
    samples = []
    header = pysam.AlignmentFile(bam, 'rb').header
    rgs = header.get('RG', None)
    if rgs:
        for rg in rgs:
            sm = rg.get('SM', None)
            if sm:
                samples.append(sm)
    return samples

def get_platform_list(bam):
    """
    input is a BAM file
    returns a list of platforms defined in the Header section of the input BAM
    """
    platforms = []
    header = pysam.AlignmentFile(bam, 'rb').header
    rgs = header.get('RG', None)
    if rgs:
        for rg in rgs:
            sm = rg.get('PL', None)
            if sm:
                platforms.append(sm)
    return platforms
