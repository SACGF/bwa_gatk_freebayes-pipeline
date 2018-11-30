"""Core Python library for ACRF Cancer Genomics Facility"""

import HTSeq
import os
import pandas

from monkey_patch import htseq_monkey_patch, pandas_monkey_patch

if os.getenv("MONKEY_PATCH_HTSEQ"):
    htseq_monkey_patch.monkey_patch_HTSeq(HTSeq)

if os.getenv("MONKEY_PATCH_PANDAS"):
    pandas_monkey_patch.monkey_patch_pandas(pandas)
