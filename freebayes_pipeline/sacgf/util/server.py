#!/usr/bin/env python

"""

"""

import os
from socket import gethostname

def is_on_server():
    """
    return True if hostname is one of tango, sacgf and bigmem
    """
    hostname = gethostname()
    if hostname in ['sacgf.ersa.edu.au',
                        'bigmem-head-01', 'bigmem512-02', 'bigmem512-03',
                        'bigmem1024-1.tizard.ersa.edu.au',
                        'bigmem512-1.tizard.ersa.edu.au',
                        'tango-head-01.ersa.edu.au']:
        return True
    elif hostname.startswith('tango'):
        #tango's compute nodes' names; just don't name your machine with "tango*", OK?
        return True
    elif hostname.startswith('bigmem'):
        #bigmem's compute nodes' names; just don't name your machine with "bigmem*", OK?
        return True
    else:
        return False

def get_sacgf_dir():
    """
    return '/data/sacgf/' if it is on server,
    else return shell environment variable $SACGF_DIR
    """
    if is_on_server():
        return '/data/sacgf/'
    else:
        return os.environ.get('SACGF_DIR')
