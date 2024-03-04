#!/usr/bin/env python
"""
TADbit based code. It contains slightly modified quality_plot function from TADbit.
It was modified to produce machine-readable file with stats including percentage of diggested sites.
Also, a generation of a figure in multiple formats was added as option.
"""

import re
import os
import bz2
import gzip
import zipfile
import tarfile
import argparse

from io import open, TextIOWrapper, BufferedReader
from warnings import catch_warnings, simplefilter
from itertools import zip_longest
from subprocess import Popen, PIPE
from collections import OrderedDict
from multiprocessing import cpu_count

from numpy import nanstd, nanmean, linspace, nansum
import matplotlib.pyplot as plt


try:
    basestring
except NameError:
    basestring = str


class RE_dict(dict):
    def __getitem__(self, i):
        try:
            return super(RE_dict, self).__getitem__(i)
        except KeyError:
            for nam in self:
                if nam.lower() == i.lower():
                    return self[nam]
            raise KeyError('Restriction Enzyme %s not found\n' % (i))


RESTRICTION_ENZYMES = RE_dict([('AanI'       , 'TTA|TAA'                     ),
                               ('AarI'       , 'CACCTGC|'                    ),
                               ('AasI'       , 'GACNNNN|NNGTC'               ),
                               ('AatII'      , 'GACGT|C'                     ),
                               ('AbaSI'      , 'C|'                          ),
                               ('AbsI'       , 'CC|TCGAGG'                   ),
                               ('Acc16I'     , 'TGC|GCA'                     ),
                               ('Acc36I'     , 'ACCTGC|'                     ),
                               ('Acc65I'     , 'G|GTACC'                     ),
                               ('AccB1I'     , 'G|GYRCC'                     ),
                               ('AccB7I'     , 'CCANNNN|NTGG'                ),
                               ('AccBSI'     , 'CCG|CTC'                     ),
                               ('AccI'       , 'GT|MKAC'                     ),
                               ('AccII'      , 'CG|CG'                       ),
                               ('AccIII'     , 'T|CCGGA'                     ),
                               ('AceIII'     , 'CAGCTC|'                     ),
                               ('AciI'       , 'C|CGC'                       ),
                               ('AclI'       , 'AA|CGTT'                     ),
                               ('AclWI'      , 'GGATC|'                      ),
                               ('AcoI'       , 'Y|GGCCR'                     ),
                               ('AcsI'       , 'R|AATTY'                     ),
                               ('AcuI'       , 'CTGAAG|'                     ),
                               ('AcvI'       , 'CAC|GTG'                     ),
                               ('AcyI'       , 'GR|CGYC'                     ),
                               ('AdeI'       , 'CACNNN|GTG'                  ),
                               ('AfaI'       , 'GT|AC'                       ),
                               ('AfeI'       , 'AGC|GCT'                     ),
                               ('AfiI'       , 'CCNNNNN|NNGG'                ),
                               ('AflII'      , 'C|TTAAG'                     ),
                               ('AflIII'     , 'A|CRYGT'                     ),
                               ('AgeI'       , 'A|CCGGT'                     ),
                               ('AgsI'       , 'TTS|AA'                      ),
                               ('AhaIII'     , 'TTT|AAA'                     ),
                               ('AhdI'       , 'GACNNN|NNGTC'                ),
                               ('AhlI'       , 'A|CTAGT'                     ),
                               ('AjiI'       , 'CAC|GTC'                     ),
                               ('AjnI'       , '|CCWGG'                      ),
                               ('AjuI'       , 'GAANNNN|NNNTTGG'             ),
                               ('AleI'       , 'CACNN|NNGTG'                 ),
                               ('AlfI'       , 'GC|ANNNNNNTGC'               ),
                               ('AloI'       , 'GAACNN|NNNNTCC'              ),
                               ('AluBI'      , 'AG|CT'                       ),
                               ('AluI'       , 'AG|CT'                       ),
                               ('Alw21I'     , 'GWGCW|C'                     ),
                               ('Alw26I'     , 'GTCTC|'                      ),
                               ('Alw44I'     , 'G|TGCAC'                     ),
                               ('AlwFI'      , 'GAAAYNNNNNRTG|GAAAYNNNNNRTG' ),
                               ('AlwI'       , 'GGATC|'                      ),
                               ('AlwNI'      , 'CAGNNN|CTG'                  ),
                               ('Ama87I'     , 'C|YCGRG'                     ),
                               ('Aor13HI'    , 'T|CCGGA'                     ),
                               ('Aor51HI'    , 'AGC|GCT'                     ),
                               ('AoxI'       , '|GGCC'                       ),
                               ('ApaBI'      , 'GCANNNNN|TGC'                ),
                               ('ApaI'       , 'GGGCC|C'                     ),
                               ('ApaLI'      , 'G|TGCAC'                     ),
                               ('ApeKI'      , 'G|CWGC'                      ),
                               ('ApoI'       , 'R|AATTY'                     ),
                               ('ApyPI'      , 'ATCGAC|'                     ),
                               ('AquII'      , 'GCCGNAC|'                    ),
                               ('AquIII'     , 'GAGGAG|'                     ),
                               ('AquIV'      , 'GRGGAAG|'                    ),
                               ('ArsI'       , 'GACNN|NNNNTTYG'              ),
                               ('AscI'       , 'GG|CGCGCC'                   ),
                               ('AseI'       , 'AT|TAAT'                     ),
                               ('Asi256I'    , 'G|ATC'                       ),
                               ('AsiGI'      , 'A|CCGGT'                     ),
                               ('AsiSI'      , 'GCGAT|CGC'                   ),
                               ('Asp700I'    , 'GAANN|NNTTC'                 ),
                               ('Asp718I'    , 'G|GTACC'                     ),
                               ('AspA2I'     , 'C|CTAGG'                     ),
                               ('AspBHI'     , 'YSCNS|'                      ),
                               ('AspLEI'     , 'GCG|C'                       ),
                               ('AspS9I'     , 'G|GNCC'                      ),
                               ('AssI'       , 'AGT|ACT'                     ),
                               ('AsuC2I'     , 'CC|SGG'                      ),
                               ('AsuHPI'     , 'GGTGA|'                      ),
                               ('AsuI'       , 'G|GNCC'                      ),
                               ('AsuII'      , 'TT|CGAA'                     ),
                               ('AsuNHI'     , 'G|CTAGC'                     ),
                               ('AvaI'       , 'C|YCGRG'                     ),
                               ('AvaII'      , 'G|GWCC'                      ),
                               ('AvaIII'     , 'ATGCAT|ATGCAT'               ),
                               ('AvrII'      , 'C|CTAGG'                     ),
                               ('AxyI'       , 'CC|TNAGG'                    ),
                               ('BaeGI'      , 'GKGCM|C'                     ),
                               ('BaeI'       , 'A|CNNNNGTAYC'                ),
                               ('BalI'       , 'TGG|CCA'                     ),
                               ('BamHI'      , 'G|GATCC'                     ),
                               ('BanI'       , 'G|GYRCC'                     ),
                               ('BanII'      , 'GRGCY|C'                     ),
                               ('BarI'       , 'GAAGNN|NNNNTAC'              ),
                               ('BasI'       , 'CCANNNN|NTGG'                ),
                               ('BauI'       , 'C|ACGAG'                     ),
                               ('Bbr7I'      , 'GAAGAC|'                     ),
                               ('BbrPI'      , 'CAC|GTG'                     ),
                               ('BbsI'       , 'GAAGAC|'                     ),
                               ('Bbv12I'     , 'GWGCW|C'                     ),
                               ('BbvCI'      , 'CC|TCAGC'                    ),
                               ('BbvI'       , 'GCAGC|'                      ),
                               ('BbvII'      , 'GAAGAC|'                     ),
                               ('BccI'       , 'CCATC|'                      ),
                               ('Bce83I'     , 'CTTGAG|'                     ),
                               ('BceAI'      , 'ACGGC|'                      ),
                               ('BcefI'      , 'ACGGC|'                      ),
                               ('BcgI'       , 'CG|ANNNNNNTGC'               ),
                               ('BciT130I'   , 'CC|WGG'                      ),
                               ('BciVI'      , 'GTATCC|'                     ),
                               ('BclI'       , 'T|GATCA'                     ),
                               ('BcnI'       , 'CC|SGG'                      ),
                               ('BcoDI'      , 'GTCTC|'                      ),
                               ('BcuI'       , 'A|CTAGT'                     ),
                               ('BdaI'       , 'TG|ANNNNNNTCA'               ),
                               ('BetI'       , 'W|CCGGW'                     ),
                               ('BfaI'       , 'C|TAG'                       ),
                               ('BfiI'       , 'ACTGGG|'                     ),
                               ('BfmI'       , 'C|TRYAG'                     ),
                               ('BfoI'       , 'RGCGC|Y'                     ),
                               ('BfrI'       , 'C|TTAAG'                     ),
                               ('BfuAI'      , 'ACCTGC|'                     ),
                               ('BfuCI'      , '|GATC'                       ),
                               ('BfuI'       , 'GTATCC|'                     ),
                               ('BglI'       , 'GCCNNNN|NGGC'                ),
                               ('BglII'      , 'A|GATCT'                     ),
                               ('BinI'       , 'GGATC|'                      ),
                               ('BisI'       , 'GC|NGC'                      ),
                               ('BlnI'       , 'C|CTAGG'                     ),
                               ('BlpI'       , 'GC|TNAGC'                    ),
                               ('BlsI'       , 'GCN|GC'                      ),
                               ('BmcAI'      , 'AGT|ACT'                     ),
                               ('Bme1390I'   , 'CC|NGG'                      ),
                               ('Bme18I'     , 'G|GWCC'                      ),
                               ('BmeDI'      , 'C|'                          ),
                               ('BmeRI'      , 'GACNNN|NNGTC'                ),
                               ('BmeT110I'   , 'C|YCGRG'                     ),
                               ('BmgBI'      , 'CAC|GTC'                     ),
                               ('BmgI'       , 'GKGCCC|GKGCCC'               ),
                               ('BmgT120I'   , 'G|GNCC'                      ),
                               ('BmiI'       , 'GGN|NCC'                     ),
                               ('BmrFI'      , 'CC|NGG'                      ),
                               ('BmrI'       , 'ACTGGG|'                     ),
                               ('BmsI'       , 'GCATC|'                      ),
                               ('BmtI'       , 'GCTAG|C'                     ),
                               ('BmuI'       , 'ACTGGG|'                     ),
                               ('BoxI'       , 'GACNN|NNGTC'                 ),
                               ('BpiI'       , 'GAAGAC|'                     ),
                               ('BplI'       , 'GAG|NNNNNCTC'                ),
                               ('BpmI'       , 'CTGGAG|'                     ),
                               ('Bpu10I'     , 'CC|TNAGC'                    ),
                               ('Bpu1102I'   , 'GC|TNAGC'                    ),
                               ('Bpu14I'     , 'TT|CGAA'                     ),
                               ('BpuEI'      , 'CTTGAG|'                     ),
                               ('BpuMI'      , 'CC|SGG'                      ),
                               ('BpvUI'      , 'CGAT|CG'                     ),
                               ('Bsa29I'     , 'AT|CGAT'                     ),
                               ('BsaAI'      , 'YAC|GTR'                     ),
                               ('BsaBI'      , 'GATNN|NNATC'                 ),
                               ('BsaHI'      , 'GR|CGYC'                     ),
                               ('BsaI'       , 'GGTCTC|'                     ),
                               ('BsaJI'      , 'C|CNNGG'                     ),
                               ('BsaWI'      , 'W|CCGGW'                     ),
                               ('BsaXI'      , 'AC|NNNNNCTCC'                ),
                               ('BsbI'       , 'CAACAC|'                     ),
                               ('Bsc4I'      , 'CCNNNNN|NNGG'                ),
                               ('BscAI'      , 'GCATC|'                      ),
                               ('BscGI'      , 'CCCGT|CCCGT'                 ),
                               ('Bse118I'    , 'R|CCGGY'                     ),
                               ('Bse1I'      , 'ACTGG|'                      ),
                               ('Bse21I'     , 'CC|TNAGG'                    ),
                               ('Bse3DI'     , 'GCAATG|'                     ),
                               ('Bse8I'      , 'GATNN|NNATC'                 ),
                               ('BseAI'      , 'T|CCGGA'                     ),
                               ('BseBI'      , 'CC|WGG'                      ),
                               ('BseCI'      , 'AT|CGAT'                     ),
                               ('BseDI'      , 'C|CNNGG'                     ),
                               ('BseGI'      , 'GGATG|'                      ),
                               ('BseJI'      , 'GATNN|NNATC'                 ),
                               ('BseLI'      , 'CCNNNNN|NNGG'                ),
                               ('BseMI'      , 'GCAATG|'                     ),
                               ('BseMII'     , 'CTCAG|'                      ),
                               ('BseNI'      , 'ACTGG|'                      ),
                               ('BsePI'      , 'G|CGCGC'                     ),
                               ('BseRI'      , 'GAGGAG|'                     ),
                               ('BseSI'      , 'GKGCM|C'                     ),
                               ('BseX3I'     , 'C|GGCCG'                     ),
                               ('BseXI'      , 'GCAGC|'                      ),
                               ('BseYI'      , 'C|CCAGC'                     ),
                               ('BsgI'       , 'GTGCAG|'                     ),
                               ('Bsh1236I'   , 'CG|CG'                       ),
                               ('Bsh1285I'   , 'CGRY|CG'                     ),
                               ('BshFI'      , 'GG|CC'                       ),
                               ('BshNI'      , 'G|GYRCC'                     ),
                               ('BshTI'      , 'A|CCGGT'                     ),
                               ('BshVI'      , 'AT|CGAT'                     ),
                               ('BsiEI'      , 'CGRY|CG'                     ),
                               ('BsiHKAI'    , 'GWGCW|C'                     ),
                               ('BsiHKCI'    , 'C|YCGRG'                     ),
                               ('BsiI'       , 'C|ACGAG'                     ),
                               ('BsiSI'      , 'C|CGG'                       ),
                               ('BsiWI'      , 'C|GTACG'                     ),
                               ('BsiYI'      , 'CCNNNNN|NNGG'                ),
                               ('BslFI'      , 'GGGAC|'                      ),
                               ('BslI'       , 'CCNNNNN|NNGG'                ),
                               ('BsmAI'      , 'GTCTC|'                      ),
                               ('BsmBI'      , 'CGTCTC|'                     ),
                               ('BsmFI'      , 'GGGAC|'                      ),
                               ('BsmI'       , 'GAATGC|'                     ),
                               ('BsnI'       , 'GG|CC'                       ),
                               ('Bso31I'     , 'GGTCTC|'                     ),
                               ('BsoBI'      , 'C|YCGRG'                     ),
                               ('Bsp119I'    , 'TT|CGAA'                     ),
                               ('Bsp120I'    , 'G|GGCCC'                     ),
                               ('Bsp1286I'   , 'GDGCH|C'                     ),
                               ('Bsp13I'     , 'T|CCGGA'                     ),
                               ('Bsp1407I'   , 'T|GTACA'                     ),
                               ('Bsp143I'    , '|GATC'                       ),
                               ('Bsp1720I'   , 'GC|TNAGC'                    ),
                               ('Bsp19I'     , 'C|CATGG'                     ),
                               ('Bsp24I'     , 'GACN|NNNNNTGG'               ),
                               ('Bsp68I'     , 'TCG|CGA'                     ),
                               ('BspACI'     , 'C|CGC'                       ),
                               ('BspCNI'     , 'CTCAG|'                      ),
                               ('BspD6I'     , 'GACTC|'                      ),
                               ('BspDI'      , 'AT|CGAT'                     ),
                               ('BspEI'      , 'T|CCGGA'                     ),
                               ('BspFNI'     , 'CG|CG'                       ),
                               ('BspGI'      , 'CTGGAC|CTGGAC'               ),
                               ('BspHI'      , 'T|CATGA'                     ),
                               ('BspLI'      , 'GGN|NCC'                     ),
                               ('BspLU11I'   , 'A|CATGT'                     ),
                               ('BspMI'      , 'ACCTGC|'                     ),
                               ('BspMII'     , 'T|CCGGA'                     ),
                               ('BspNCI'     , 'CCAGA|CCAGA'                 ),
                               ('BspOI'      , 'GCTAG|C'                     ),
                               ('BspPI'      , 'GGATC|'                      ),
                               ('BspQI'      , 'GCTCTTC|'                    ),
                               ('BspT104I'   , 'TT|CGAA'                     ),
                               ('BspT107I'   , 'G|GYRCC'                     ),
                               ('BspTI'      , 'C|TTAAG'                     ),
                               ('BsrBI'      , 'CCG|CTC'                     ),
                               ('BsrDI'      , 'GCAATG|'                     ),
                               ('BsrFI'      , 'R|CCGGY'                     ),
                               ('BsrGI'      , 'T|GTACA'                     ),
                               ('BsrI'       , 'ACTGG|'                      ),
                               ('BsrSI'      , 'ACTGG|'                      ),
                               ('BssAI'      , 'R|CCGGY'                     ),
                               ('BssECI'     , 'C|CNNGG'                     ),
                               ('BssHII'     , 'G|CGCGC'                     ),
                               ('BssKI'      , '|CCNGG'                      ),
                               ('BssMI'      , '|GATC'                       ),
                               ('BssNAI'     , 'GTA|TAC'                     ),
                               ('BssNI'      , 'GR|CGYC'                     ),
                               ('BssSI'      , 'C|ACGAG'                     ),
                               ('BssT1I'     , 'C|CWWGG'                     ),
                               ('Bst1107I'   , 'GTA|TAC'                     ),
                               ('Bst2BI'     , 'C|ACGAG'                     ),
                               ('Bst2UI'     , 'CC|WGG'                      ),
                               ('Bst4CI'     , 'ACN|GT'                      ),
                               ('Bst6I'      , 'CTCTTC|'                     ),
                               ('BstACI'     , 'GR|CGYC'                     ),
                               ('BstAFI'     , 'C|TTAAG'                     ),
                               ('BstAPI'     , 'GCANNNN|NTGC'                ),
                               ('BstAUI'     , 'T|GTACA'                     ),
                               ('BstBAI'     , 'YAC|GTR'                     ),
                               ('BstBI'      , 'TT|CGAA'                     ),
                               ('BstC8I'     , 'GCN|NGC'                     ),
                               ('BstDEI'     , 'C|TNAG'                      ),
                               ('BstDSI'     , 'C|CRYGG'                     ),
                               ('BstEII'     , 'G|GTNACC'                    ),
                               ('BstENI'     , 'CCTNN|NNNAGG'                ),
                               ('BstF5I'     , 'GGATG|'                      ),
                               ('BstFNI'     , 'CG|CG'                       ),
                               ('BstH2I'     , 'RGCGC|Y'                     ),
                               ('BstHHI'     , 'GCG|C'                       ),
                               ('BstKTI'     , 'GAT|C'                       ),
                               ('BstMAI'     , 'GTCTC|'                      ),
                               ('BstMBI'     , '|GATC'                       ),
                               ('BstMCI'     , 'CGRY|CG'                     ),
                               ('BstMWI'     , 'GCNNNNN|NNGC'                ),
                               ('BstNI'      , 'CC|WGG'                      ),
                               ('BstNSI'     , 'RCATG|Y'                     ),
                               ('BstOI'      , 'CC|WGG'                      ),
                               ('BstPAI'     , 'GACNN|NNGTC'                 ),
                               ('BstPI'      , 'G|GTNACC'                    ),
                               ('BstSCI'     , '|CCNGG'                      ),
                               ('BstSFI'     , 'C|TRYAG'                     ),
                               ('BstSLI'     , 'GKGCM|C'                     ),
                               ('BstSNI'     , 'TAC|GTA'                     ),
                               ('BstUI'      , 'CG|CG'                       ),
                               ('BstV1I'     , 'GCAGC|'                      ),
                               ('BstV2I'     , 'GAAGAC|'                     ),
                               ('BstX2I'     , 'R|GATCY'                     ),
                               ('BstXI'      , 'CCANNNNN|NTGG'               ),
                               ('BstYI'      , 'R|GATCY'                     ),
                               ('BstZ17I'    , 'GTA|TAC'                     ),
                               ('BstZI'      , 'C|GGCCG'                     ),
                               ('Bsu15I'     , 'AT|CGAT'                     ),
                               ('Bsu36I'     , 'CC|TNAGG'                    ),
                               ('BsuI'       , 'GTATCC|'                     ),
                               ('BsuRI'      , 'GG|CC'                       ),
                               ('BtgI'       , 'C|CRYGG'                     ),
                               ('BtgZI'      , 'GCGATG|'                     ),
                               ('BthCI'      , 'GCNG|C'                      ),
                               ('BtrI'       , 'CAC|GTC'                     ),
                               ('BtsCI'      , 'GGATG|'                      ),
                               ('BtsI'       , 'GCAGTG|'                     ),
                               ('BtsIMutI'   , 'CAGTG|'                      ),
                               ('BtuMI'      , 'TCG|CGA'                     ),
                               ('BveI'       , 'ACCTGC|'                     ),
                               ('Cac8I'      , 'GCN|NGC'                     ),
                               ('CaiI'       , 'CAGNNN|CTG'                  ),
                               ('CauII'      , 'CC|SGG'                      ),
                               ('CchII'      , 'GGARGA|'                     ),
                               ('CchIII'     , 'CCCAAG|'                     ),
                               ('CciI'       , 'T|CATGA'                     ),
                               ('CciNI'      , 'GC|GGCCGC'                   ),
                               ('Cdi630V'    , 'CAAAAA|CAAAAA'               ),
                               ('CdiI'       , 'CATC|G'                      ),
                               ('CdpI'       , 'GCGGAG|'                     ),
                               ('CfoI'       , 'GCG|C'                       ),
                               ('Cfr10I'     , 'R|CCGGY'                     ),
                               ('Cfr13I'     , 'G|GNCC'                      ),
                               ('Cfr42I'     , 'CCGC|GG'                     ),
                               ('Cfr9I'      , 'C|CCGGG'                     ),
                               ('CfrI'       , 'Y|GGCCR'                     ),
                               ('Cgl13032I'  , 'GGCGCA|GGCGCA'               ),
                               ('Cgl13032II' , 'ACGABGG|ACGABGG'             ),
                               ('ChaI'       , 'GATC|'                       ),
                               ('CjeFIII'    , 'GCAAGG|GCAAGG'               ),
                               ('CjeFV'      , 'GGRCA|GGRCA'                 ),
                               ('CjeI'       , 'CCA|NNNNNNGT'                ),
                               ('CjeNII'     , 'GAGNNNNNGT|GAGNNNNNGT'       ),
                               ('CjeNIII'    , 'GKAAYG|'                     ),
                               ('CjeP659IV'  , 'CACNNNNNNNGAA|CACNNNNNNNGAA' ),
                               ('CjePI'      , 'CCANN|NNNNNTC'               ),
                               ('CjuI'       , 'CAYNNNNNRTG|CAYNNNNNRTG'     ),
                               ('CjuII'      , 'CAYNNNNNCTC|CAYNNNNNCTC'     ),
                               ('ClaI'       , 'AT|CGAT'                     ),
                               ('CpoI'       , 'CG|GWCCG'                    ),
                               ('CseI'       , 'GACGC|'                      ),
                               ('CsiI'       , 'A|CCWGGT'                    ),
                               ('Csp6I'      , 'G|TAC'                       ),
                               ('CspAI'      , 'A|CCGGT'                     ),
                               ('CspCI'      , 'C|AANNNNNGTGG'               ),
                               ('CspI'       , 'CG|GWCCG'                    ),
                               ('CstMI'      , 'AAGGAG|'                     ),
                               ('CviAII'     , 'C|ATG'                       ),
                               ('CviJI'      , 'RG|CY'                       ),
                               ('CviKI_1'    , 'RG|CY'                       ),
                               ('CviQI'      , 'G|TAC'                       ),
                               ('CviRI'      , 'TG|CA'                       ),
                               ('DdeI'       , 'C|TNAG'                      ),
                               ('DinI'       , 'GGC|GCC'                     ),
                               ('DpnI'       , 'GA|TC'                       ),
                               ('DpnII'      , '|GATC'                       ),
                               ('DraI'       , 'TTT|AAA'                     ),
                               ('DraII'      , 'RG|GNCCY'                    ),
                               ('DraIII'     , 'CACNNN|GTG'                  ),
                               ('DraRI'      , 'CAAGNAC|'                    ),
                               ('DrdI'       , 'GACNNNN|NNGTC'               ),
                               ('DrdII'      , 'GAACCA|GAACCA'               ),
                               ('DriI'       , 'GACNNN|NNGTC'                ),
                               ('DsaI'       , 'C|CRYGG'                     ),
                               ('DseDI'      , 'GACNNNN|NNGTC'               ),
                               ('EaeI'       , 'Y|GGCCR'                     ),
                               ('EagI'       , 'C|GGCCG'                     ),
                               ('Eam1104I'   , 'CTCTTC|'                     ),
                               ('Eam1105I'   , 'GACNNN|NNGTC'                ),
                               ('EarI'       , 'CTCTTC|'                     ),
                               ('EciI'       , 'GGCGGA|'                     ),
                               ('Ecl136II'   , 'GAG|CTC'                     ),
                               ('EclXI'      , 'C|GGCCG'                     ),
                               ('Eco105I'    , 'TAC|GTA'                     ),
                               ('Eco130I'    , 'C|CWWGG'                     ),
                               ('Eco147I'    , 'AGG|CCT'                     ),
                               ('Eco24I'     , 'GRGCY|C'                     ),
                               ('Eco31I'     , 'GGTCTC|'                     ),
                               ('Eco32I'     , 'GAT|ATC'                     ),
                               ('Eco47I'     , 'G|GWCC'                      ),
                               ('Eco47III'   , 'AGC|GCT'                     ),
                               ('Eco52I'     , 'C|GGCCG'                     ),
                               ('Eco53kI'    , 'GAG|CTC'                     ),
                               ('Eco57I'     , 'CTGAAG|'                     ),
                               ('Eco57MI'    , 'CTGRAG|'                     ),
                               ('Eco72I'     , 'CAC|GTG'                     ),
                               ('Eco81I'     , 'CC|TNAGG'                    ),
                               ('Eco88I'     , 'C|YCGRG'                     ),
                               ('Eco91I'     , 'G|GTNACC'                    ),
                               ('EcoHI'      , '|CCSGG'                      ),
                               ('EcoICRI'    , 'GAG|CTC'                     ),
                               ('EcoNI'      , 'CCTNN|NNNAGG'                ),
                               ('EcoO109I'   , 'RG|GNCCY'                    ),
                               ('EcoO65I'    , 'G|GTNACC'                    ),
                               ('EcoRI'      , 'G|AATTC'                     ),
                               ('EcoRII'     , '|CCWGG'                      ),
                               ('EcoRV'      , 'GAT|ATC'                     ),
                               ('EcoT14I'    , 'C|CWWGG'                     ),
                               ('EcoT22I'    , 'ATGCA|T'                     ),
                               ('EcoT38I'    , 'GRGCY|C'                     ),
                               ('EgeI'       , 'GGC|GCC'                     ),
                               ('EheI'       , 'GGC|GCC'                     ),
                               ('ErhI'       , 'C|CWWGG'                     ),
                               ('EsaBC3I'    , 'TC|GA'                       ),
                               ('EsaSSI'     , 'GACCAC|GACCAC'               ),
                               ('Esp3I'      , 'CGTCTC|'                     ),
                               ('EspI'       , 'GC|TNAGC'                    ),
                               ('FaeI'       , 'CATG|'                       ),
                               ('FaiI'       , 'YA|TR'                       ),
                               ('FalI'       , 'AAG|NNNNNCTT'                ),
                               ('FaqI'       , 'GGGAC|'                      ),
                               ('FatI'       , '|CATG'                       ),
                               ('FauI'       , 'CCCGC|'                      ),
                               ('FauNDI'     , 'CA|TATG'                     ),
                               ('FbaI'       , 'T|GATCA'                     ),
                               ('FblI'       , 'GT|MKAC'                     ),
                               ('FinI'       , 'GGGAC|GGGAC'                 ),
                               ('FmuI'       , 'GGNC|C'                      ),
                               ('Fnu4HI'     , 'GC|NGC'                      ),
                               ('FnuDII'     , 'CG|CG'                       ),
                               ('FokI'       , 'GGATG|'                      ),
                               ('FriOI'      , 'GRGCY|C'                     ),
                               ('FseI'       , 'GGCCGG|CC'                   ),
                               ('Fsp4HI'     , 'GC|NGC'                      ),
                               ('FspAI'      , 'RTGC|GCAY'                   ),
                               ('FspBI'      , 'C|TAG'                       ),
                               ('FspEI'      , 'CC|'                         ),
                               ('FspI'       , 'TGC|GCA'                     ),
                               ('GauT27I'    , 'CGCGCAGG|CGCGCAGG'           ),
                               ('GdiII'      , 'C|GGCCR'                     ),
                               ('GlaI'       , 'GC|GC'                       ),
                               ('GluI'       , 'GC|NGC'                      ),
                               ('GsaI'       , 'CCCAG|C'                     ),
                               ('GsuI'       , 'CTGGAG|'                     ),
                               ('HaeI'       , 'WGG|CCW'                     ),
                               ('HaeII'      , 'RGCGC|Y'                     ),
                               ('HaeIII'     , 'GG|CC'                       ),
                               ('HapII'      , 'C|CGG'                       ),
                               ('HauII'      , 'TGGCCA|'                     ),
                               ('HgaI'       , 'GACGC|'                      ),
                               ('HgiAI'      , 'GWGCW|C'                     ),
                               ('HgiCI'      , 'G|GYRCC'                     ),
                               ('HgiEII'     , 'ACCNNNNNNGGT|ACCNNNNNNGGT'   ),
                               ('HgiJII'     , 'GRGCY|C'                     ),
                               ('HhaI'       , 'GCG|C'                       ),
                               ('Hin1I'      , 'GR|CGYC'                     ),
                               ('Hin1II'     , 'CATG|'                       ),
                               ('Hin4I'      , 'GAY|NNNNNVTC'                ),
                               ('Hin4II'     , 'CCTTC|'                      ),
                               ('Hin6I'      , 'G|CGC'                       ),
                               ('HinP1I'     , 'G|CGC'                       ),
                               ('HincII'     , 'GTY|RAC'                     ),
                               ('HindII'     , 'GTY|RAC'                     ),
                               ('HindIII'    , 'A|AGCTT'                     ),
                               ('HinfI'      , 'G|ANTC'                      ),
                               ('HpaI'       , 'GTT|AAC'                     ),
                               ('HpaII'      , 'C|CGG'                       ),
                               ('HphI'       , 'GGTGA|'                      ),
                               ('Hpy166II'   , 'GTN|NAC'                     ),
                               ('Hpy178III'  , 'TC|NNGA'                     ),
                               ('Hpy188I'    , 'TCN|GA'                      ),
                               ('Hpy188III'  , 'TC|NNGA'                     ),
                               ('Hpy8I'      , 'GTN|NAC'                     ),
                               ('Hpy99I'     , 'CGWCG|'                      ),
                               ('Hpy99XIII'  , 'GCCTA|GCCTA'                 ),
                               ('Hpy99XIV'   , 'GGWTAA|GGWTAA'               ),
                               ('HpyAV'      , 'CCTTC|'                      ),
                               ('HpyCH4III'  , 'ACN|GT'                      ),
                               ('HpyCH4IV'   , 'A|CGT'                       ),
                               ('HpyCH4V'    , 'TG|CA'                       ),
                               ('HpyF10VI'   , 'GCNNNNN|NNGC'                ),
                               ('HpyF3I'     , 'C|TNAG'                      ),
                               ('HpySE526I'  , 'A|CGT'                       ),
                               ('Hsp92I'     , 'GR|CGYC'                     ),
                               ('Hsp92II'    , 'CATG|'                       ),
                               ('HspAI'      , 'G|CGC'                       ),
                               ('Jma19592I'  , 'GTATNAC|GTATNAC'             ),
                               ('KasI'       , 'G|GCGCC'                     ),
                               ('KflI'       , 'GG|GWCCC'                    ),
                               ('Kpn2I'      , 'T|CCGGA'                     ),
                               ('KpnI'       , 'GGTAC|C'                     ),
                               ('KroI'       , 'G|CCGGC'                     ),
                               ('Ksp22I'     , 'T|GATCA'                     ),
                               ('Ksp632I'    , 'CTCTTC|'                     ),
                               ('KspAI'      , 'GTT|AAC'                     ),
                               ('KspI'       , 'CCGC|GG'                     ),
                               ('Kzo9I'      , '|GATC'                       ),
                               ('LguI'       , 'GCTCTTC|'                    ),
                               ('LpnI'       , 'RGC|GCY'                     ),
                               ('LpnPI'      , 'CCDG|'                       ),
                               ('Lsp1109I'   , 'GCAGC|'                      ),
                               ('LweI'       , 'GCATC|'                      ),
                               ('MabI'       , 'A|CCWGGT'                    ),
                               ('MaeI'       , 'C|TAG'                       ),
                               ('MaeII'      , 'A|CGT'                       ),
                               ('MaeIII'     , '|GTNAC'                      ),
                               ('MalI'       , 'GA|TC'                       ),
                               ('MaqI'       , 'CRTTGAC|'                    ),
                               ('MauBI'      , 'CG|CGCGCG'                   ),
                               ('MbiI'       , 'CCG|CTC'                     ),
                               ('MboI'       , '|GATC'                       ),
                               ('MboII'      , 'GAAGA|'                      ),
                               ('McaTI'      , 'GCGC|GC'                     ),
                               ('McrI'       , 'CGRY|CG'                     ),
                               ('MfeI'       , 'C|AATTG'                     ),
                               ('MflI'       , 'R|GATCY'                     ),
                               ('MhlI'       , 'GDGCH|C'                     ),
                               ('MjaIV'      , 'GTNNAC|GTNNAC'               ),
                               ('MkaDII'     , 'GAGAYGT|GAGAYGT'             ),
                               ('MlsI'       , 'TGG|CCA'                     ),
                               ('MluCI'      , '|AATT'                       ),
                               ('MluI'       , 'A|CGCGT'                     ),
                               ('MluNI'      , 'TGG|CCA'                     ),
                               ('Mly113I'    , 'GG|CGCC'                     ),
                               ('MlyI'       , 'GAGTC|'                      ),
                               ('MmeI'       , 'TCCRAC|'                     ),
                               ('MnlI'       , 'CCTC|'                       ),
                               ('Mph1103I'   , 'ATGCA|T'                     ),
                               ('MreI'       , 'CG|CCGGCG'                   ),
                               ('MroI'       , 'T|CCGGA'                     ),
                               ('MroNI'      , 'G|CCGGC'                     ),
                               ('MroXI'      , 'GAANN|NNTTC'                 ),
                               ('MscI'       , 'TGG|CCA'                     ),
                               ('MseI'       , 'T|TAA'                       ),
                               ('MslI'       , 'CAYNN|NNRTG'                 ),
                               ('Msp20I'     , 'TGG|CCA'                     ),
                               ('MspA1I'     , 'CMG|CKG'                     ),
                               ('MspCI'      , 'C|TTAAG'                     ),
                               ('MspI'       , 'C|CGG'                       ),
                               ('MspJI'      , 'CNNR|'                       ),
                               ('MspR9I'     , 'CC|NGG'                      ),
                               ('MssI'       , 'GTTT|AAAC'                   ),
                               ('MstI'       , 'TGC|GCA'                     ),
                               ('MunI'       , 'C|AATTG'                     ),
                               ('Mva1269I'   , 'GAATGC|'                     ),
                               ('MvaI'       , 'CC|WGG'                      ),
                               ('MvnI'       , 'CG|CG'                       ),
                               ('MvrI'       , 'CGAT|CG'                     ),
                               ('MwoI'       , 'GCNNNNN|NNGC'                ),
                               ('NaeI'       , 'GCC|GGC'                     ),
                               ('NarI'       , 'GG|CGCC'                     ),
                               ('NciI'       , 'CC|SGG'                      ),
                               ('NcoI'       , 'C|CATGG'                     ),
                               ('NdeI'       , 'CA|TATG'                     ),
                               ('NdeII'      , '|GATC'                       ),
                               ('NgoAVIII'   , '|GACNNNNNTGA'                ),
                               ('NgoMIV'     , 'G|CCGGC'                     ),
                               ('NhaXI'      , 'CAAGRAG|CAAGRAG'             ),
                               ('NheI'       , 'G|CTAGC'                     ),
                               ('NlaCI'      , 'CATCAC|'                     ),
                               ('NlaIII'     , 'CATG|'                       ),
                               ('NlaIV'      , 'GGN|NCC'                     ),
                               ('Nli3877I'   , 'CYCGR|G'                     ),
                               ('NmeAIII'    , 'GCCGAG|'                     ),
                               ('NmeDI'      , '|RCCGGY'                     ),
                               ('NmuCI'      , '|GTSAC'                      ),
                               ('NotI'       , 'GC|GGCCGC'                   ),
                               ('NruI'       , 'TCG|CGA'                     ),
                               ('NsbI'       , 'TGC|GCA'                     ),
                               ('NsiI'       , 'ATGCA|T'                     ),
                               ('NspBII'     , 'CMG|CKG'                     ),
                               ('NspI'       , 'RCATG|Y'                     ),
                               ('NspV'       , 'TT|CGAA'                     ),
                               ('OliI'       , 'CACNN|NNGTG'                 ),
                               ('PabI'       , 'GTA|C'                       ),
                               ('PacI'       , 'TTAAT|TAA'                   ),
                               ('PaeI'       , 'GCATG|C'                     ),
                               ('PaeR7I'     , 'C|TCGAG'                     ),
                               ('PagI'       , 'T|CATGA'                     ),
                               ('PalAI'      , 'GG|CGCGCC'                   ),
                               ('PasI'       , 'CC|CWGGG'                    ),
                               ('PauI'       , 'G|CGCGC'                     ),
                               ('PceI'       , 'AGG|CCT'                     ),
                               ('PciI'       , 'A|CATGT'                     ),
                               ('PciSI'      , 'GCTCTTC|'                    ),
                               ('PcsI'       , 'WCGNNNN|NNNCGW'              ),
                               ('PctI'       , 'GAATGC|'                     ),
                               ('PdiI'       , 'GCC|GGC'                     ),
                               ('PdmI'       , 'GAANN|NNTTC'                 ),
                               ('PenI'       , 'GCAGT|GCAGT'                 ),
                               ('PfeI'       , 'G|AWTC'                      ),
                               ('Pfl1108I'   , 'TCGTAG|TCGTAG'               ),
                               ('Pfl23II'    , 'C|GTACG'                     ),
                               ('PflFI'      , 'GACN|NNGTC'                  ),
                               ('PflMI'      , 'CCANNNN|NTGG'                ),
                               ('PfoI'       , 'T|CCNGGA'                    ),
                               ('PinAI'      , 'A|CCGGT'                     ),
                               ('PlaDI'      , 'CATCAG|'                     ),
                               ('Ple19I'     , 'CGAT|CG'                     ),
                               ('PleI'       , 'GAGTC|'                      ),
                               ('PluTI'      , 'GGCGC|C'                     ),
                               ('PmaCI'      , 'CAC|GTG'                     ),
                               ('PmeI'       , 'GTTT|AAAC'                   ),
                               ('PmlI'       , 'CAC|GTG'                     ),
                               ('PpiI'       , 'GAACN|NNNNCTC'               ),
                               ('PpsI'       , 'GAGTC|'                      ),
                               ('Ppu10I'     , 'A|TGCAT'                     ),
                               ('Ppu21I'     , 'YAC|GTR'                     ),
                               ('PpuMI'      , 'RG|GWCCY'                    ),
                               ('PscI'       , 'A|CATGT'                     ),
                               ('PshAI'      , 'GACNN|NNGTC'                 ),
                               ('PshBI'      , 'AT|TAAT'                     ),
                               ('PsiI'       , 'TTA|TAA'                     ),
                               ('Psp03I'     , 'GGWC|C'                      ),
                               ('Psp124BI'   , 'GAGCT|C'                     ),
                               ('Psp1406I'   , 'AA|CGTT'                     ),
                               ('Psp5II'     , 'RG|GWCCY'                    ),
                               ('Psp6I'      , '|CCWGG'                      ),
                               ('PspCI'      , 'CAC|GTG'                     ),
                               ('PspEI'      , 'G|GTNACC'                    ),
                               ('PspGI'      , '|CCWGG'                      ),
                               ('PspLI'      , 'C|GTACG'                     ),
                               ('PspN4I'     , 'GGN|NCC'                     ),
                               ('PspOMI'     , 'G|GGCCC'                     ),
                               ('PspOMII'    , 'CGCCCAR|'                    ),
                               ('PspPI'      , 'G|GNCC'                      ),
                               ('PspPPI'     , 'RG|GWCCY'                    ),
                               ('PspPRI'     , 'CCYCAG|'                     ),
                               ('PspXI'      , 'VC|TCGAGB'                   ),
                               ('PsrI'       , 'GAACNN|NNNNTAC'              ),
                               ('PssI'       , 'RGGNC|CY'                    ),
                               ('PstI'       , 'CTGCA|G'                     ),
                               ('PstNI'      , 'CAGNNN|CTG'                  ),
                               ('PsuI'       , 'R|GATCY'                     ),
                               ('PsyI'       , 'GACN|NNGTC'                  ),
                               ('PteI'       , 'G|CGCGC'                     ),
                               ('PvuI'       , 'CGAT|CG'                     ),
                               ('PvuII'      , 'CAG|CTG'                     ),
                               ('R2_BceSIV'  , '|GCAGC'                      ),
                               ('RceI'       , 'CATCGAC|'                    ),
                               ('RdeGBI'     , 'CCGCAG|CCGCAG'               ),
                               ('RdeGBII'    , 'ACCCAG|'                     ),
                               ('RdeGBIII'   , '|TGRYCA'                     ),
                               ('RflFIII'    , 'CGCCAG|CGCCAG'               ),
                               ('RgaI'       , 'GCGAT|CGC'                   ),
                               ('RigI'       , 'GGCCGG|CC'                   ),
                               ('RlaI'       , 'VCW|VCW'                     ),
                               ('RleAI'      , 'CCCACA|'                     ),
                               ('RpaB5I'     , 'CGRGGAC|'                    ),
                               ('RpaBI'      , 'CCCGCAG|'                    ),
                               ('RpaI'       , 'GTYGGAG|'                    ),
                               ('RpaTI'      , 'GRTGGAG|GRTGGAG'             ),
                               ('RruI'       , 'TCG|CGA'                     ),
                               ('RsaI'       , 'GT|AC'                       ),
                               ('RsaNI'      , 'G|TAC'                       ),
                               ('RseI'       , 'CAYNN|NNRTG'                 ),
                               ('Rsr2I'      , 'CG|GWCCG'                    ),
                               ('RsrII'      , 'CG|GWCCG'                    ),
                               ('SacI'       , 'GAGCT|C'                     ),
                               ('SacII'      , 'CCGC|GG'                     ),
                               ('SalI'       , 'G|TCGAC'                     ),
                               ('SanDI'      , 'GG|GWCCC'                    ),
                               ('SapI'       , 'GCTCTTC|'                    ),
                               ('SaqAI'      , 'T|TAA'                       ),
                               ('SatI'       , 'GC|NGC'                      ),
                               ('Sau3AI'     , '|GATC'                       ),
                               ('Sau96I'     , 'G|GNCC'                      ),
                               ('SauI'       , 'CC|TNAGG'                    ),
                               ('SbfI'       , 'CCTGCA|GG'                   ),
                               ('ScaI'       , 'AGT|ACT'                     ),
                               ('SchI'       , 'GAGTC|'                      ),
                               ('SciI'       , 'CTC|GAG'                     ),
                               ('ScrFI'      , 'CC|NGG'                      ),
                               ('SdaI'       , 'CCTGCA|GG'                   ),
                               ('SdeAI'      , 'CAGRAG|'                     ),
                               ('SdeOSI'     , '|GACNNNNRTGA'                ),
                               ('SduI'       , 'GDGCH|C'                     ),
                               ('SecI'       , 'C|CNNGG'                     ),
                               ('SelI'       , '|CGCG'                       ),
                               ('SetI'       , 'ASST|'                       ),
                               ('SexAI'      , 'A|CCWGGT'                    ),
                               ('SfaAI'      , 'GCGAT|CGC'                   ),
                               ('SfaNI'      , 'GCATC|'                      ),
                               ('SfcI'       , 'C|TRYAG'                     ),
                               ('SfeI'       , 'C|TRYAG'                     ),
                               ('SfiI'       , 'GGCCNNNN|NGGCC'              ),
                               ('SfoI'       , 'GGC|GCC'                     ),
                               ('Sfr274I'    , 'C|TCGAG'                     ),
                               ('Sfr303I'    , 'CCGC|GG'                     ),
                               ('SfuI'       , 'TT|CGAA'                     ),
                               ('SgeI'       , 'CNNG|'                       ),
                               ('SgfI'       , 'GCGAT|CGC'                   ),
                               ('SgrAI'      , 'CR|CCGGYG'                   ),
                               ('SgrBI'      , 'CCGC|GG'                     ),
                               ('SgrDI'      , 'CG|TCGACG'                   ),
                               ('SgrTI'      , 'CCDS|'                       ),
                               ('SgsI'       , 'GG|CGCGCC'                   ),
                               ('SimI'       , 'GG|GTC'                      ),
                               ('SlaI'       , 'C|TCGAG'                     ),
                               ('SmaI'       , 'CCC|GGG'                     ),
                               ('SmiI'       , 'ATTT|AAAT'                   ),
                               ('SmiMI'      , 'CAYNN|NNRTG'                 ),
                               ('SmlI'       , 'C|TYRAG'                     ),
                               ('SmoI'       , 'C|TYRAG'                     ),
                               ('SnaBI'      , 'TAC|GTA'                     ),
                               ('SnaI'       , 'GTATAC|GTATAC'               ),
                               ('Sno506I'    , 'GGCCGAG|GGCCGAG'             ),
                               ('SpeI'       , 'A|CTAGT'                     ),
                               ('SphI'       , 'GCATG|C'                     ),
                               ('SplI'       , 'C|GTACG'                     ),
                               ('SpoDI'      , 'GCGGRAG|GCGGRAG'             ),
                               ('SrfI'       , 'GCCC|GGGC'                   ),
                               ('Sse232I'    , 'CG|CCGGCG'                   ),
                               ('Sse8387I'   , 'CCTGCA|GG'                   ),
                               ('Sse8647I'   , 'AG|GWCCT'                    ),
                               ('Sse9I'      , '|AATT'                       ),
                               ('SseBI'      , 'AGG|CCT'                     ),
                               ('SsiI'       , 'C|CGC'                       ),
                               ('SspD5I'     , 'GGTGA|'                      ),
                               ('SspDI'      , 'G|GCGCC'                     ),
                               ('SspI'       , 'AAT|ATT'                     ),
                               ('SstE37I'    , 'CGAAGAC|'                    ),
                               ('SstI'       , 'GAGCT|C'                     ),
                               ('Sth132I'    , 'CCCG|'                       ),
                               ('Sth302II'   , 'CC|GG'                       ),
                               ('StrI'       , 'C|TCGAG'                     ),
                               ('StsI'       , 'GGATG|'                      ),
                               ('StuI'       , 'AGG|CCT'                     ),
                               ('StyD4I'     , '|CCNGG'                      ),
                               ('StyI'       , 'C|CWWGG'                     ),
                               ('SwaI'       , 'ATTT|AAAT'                   ),
                               ('TaaI'       , 'ACN|GT'                      ),
                               ('TaiI'       , 'ACGT|'                       ),
                               ('TaqI'       , 'T|CGA'                       ),
                               ('TaqII'      , 'GACCGA|'                     ),
                               ('TasI'       , '|AATT'                       ),
                               ('TatI'       , 'W|GTACW'                     ),
                               ('TauI'       , 'GCSG|C'                      ),
                               ('TfiI'       , 'G|AWTC'                      ),
                               ('Tru1I'      , 'T|TAA'                       ),
                               ('Tru9I'      , 'T|TAA'                       ),
                               ('TscAI'      , 'CASTG|'                      ),
                               ('TseFI'      , '|GTSAC'                      ),
                               ('TseI'       , 'G|CWGC'                      ),
                               ('TsoI'       , 'TARCCA|'                     ),
                               ('Tsp45I'     , '|GTSAC'                      ),
                               ('Tsp4CI'     , 'ACN|GT'                      ),
                               ('TspDTI'     , 'ATGAA|'                      ),
                               ('TspEI'      , '|AATT'                       ),
                               ('TspGWI'     , 'ACGGA|'                      ),
                               ('TspMI'      , 'C|CCGGG'                     ),
                               ('TspRI'      , 'CASTG|'                      ),
                               ('TssI'       , 'GAGNNNCTC|GAGNNNCTC'         ),
                               ('TstI'       , 'CACN|NNNNNTCC'               ),
                               ('TsuI'       , 'GCGAC|GCGAC'                 ),
                               ('Tth111I'    , 'GACN|NNGTC'                  ),
                               ('Tth111II'   , 'CAARCA|'                     ),
                               ('UbaF11I'    , 'TCGTA|TCGTA'                 ),
                               ('UbaF12I'    , 'CTACNNNGTC|CTACNNNGTC'       ),
                               ('UbaF13I'    , 'GAGNNNNNNCTGG|GAGNNNNNNCTGG' ),
                               ('UbaF14I'    , 'CCANNNNNTCG|CCANNNNNTCG'     ),
                               ('UbaF9I'     , 'TACNNNNNRTGT|TACNNNNNRTGT'   ),
                               ('UbaPI'      , 'CGAACG|CGAACG'               ),
                               ('UcoMSI'     , '|GAGCTC'                     ),
                               ('UnbI'       , '|GGNCC'                      ),
                               ('Van91I'     , 'CCANNNN|NTGG'                ),
                               ('Vha464I'    , 'C|TTAAG'                     ),
                               ('VneI'       , 'G|TGCAC'                     ),
                               ('VpaK11AI'   , '|GGWCC'                      ),
                               ('VpaK11BI'   , 'G|GWCC'                      ),
                               ('VspI'       , 'AT|TAAT'                     ),
                               ('WviI'       , 'CACRAG|'                     ),
                               ('XagI'       , 'CCTNN|NNNAGG'                ),
                               ('XapI'       , 'R|AATTY'                     ),
                               ('XbaI'       , 'T|CTAGA'                     ),
                               ('XceI'       , 'RCATG|Y'                     ),
                               ('XcmI'       , 'CCANNNNN|NNNNTGG'            ),
                               ('XhoI'       , 'C|TCGAG'                     ),
                               ('XhoII'      , 'R|GATCY'                     ),
                               ('XmaI'       , 'C|CCGGG'                     ),
                               ('XmaIII'     , 'C|GGCCG'                     ),
                               ('XmaJI'      , 'C|CTAGG'                     ),
                               ('XmiI'       , 'GT|MKAC'                     ),
                               ('XmnI'       , 'GAANN|NNTTC'                 ),
                               ('XspI'       , 'C|TAG'                       ),
                               ('YkrI'       , 'C|'                          ),
                               ('ZraI'       , 'GAC|GTC'                     ),
                               ('ZrmI'       , 'AGT|ACT'                     ),
                               ('Zsp2I'      , 'ATGCA|T'                     )])


def religateds(r_enzs):
    """
    returns the resulting list of all possible sequences after religation of two
    digested and repaired ends.
    """
    ligations = OrderedDict()
    for r_enz1 in r_enzs:
        for r_enz2 in r_enzs:
            site1 = RESTRICTION_ENZYMES[r_enz1]
            site2 = RESTRICTION_ENZYMES[r_enz2]
            beg1, end1 = site1.split('|')
            _, end2 = site2.split('|')
            site1 = site1.replace('|', '')
            site2 = site2.replace('|', '')
            ligations[(r_enz1, r_enz2)] = beg1 + end1[:len(end1)-len(beg1)] + end2
    return ligations


def repaired(r_enz):
    """
    returns the resulting sequence after reparation of two digested and repaired
    ends, marking dangling ends.
    """
    site = RESTRICTION_ENZYMES[r_enz]
    beg, end = site.split('|')
    site = site.replace('|', '')
    return complementary(beg + site[min(len(beg), len(end)) :
                                    max(len(beg), len(end))])


def complementary(seq):
    trs = dict([(nt1, nt2) for nt1, nt2 in zip('ATGCN', 'TACGN')])
    return ''.join([trs[s] for s in seq[::-1]])


def iupac2regex(restring):
    """
    Convert target sites with IUPAC nomenclature to regex pattern
    """
    restring = restring.replace('R', '[AG]')
    restring = restring.replace('Y', '[CT]')
    restring = restring.replace('M', '[AC]')
    restring = restring.replace('K', '[GT]')
    restring = restring.replace('S', '[CG]')
    restring = restring.replace('W', '[AT]')
    restring = restring.replace('H', '[ACT]')
    restring = restring.replace('B', '[CGT]')
    restring = restring.replace('V', '[ACG]')
    restring = restring.replace('D', '[AGT]')
    restring = restring.replace('N', '[ATGC]')
    return restring


def which(program):
    """
    stackoverflow: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def magic_open(filename, verbose=False, cpus=None):
    """
    To read uncompressed zip gzip bzip2 or tar.xx files

    :param filename: either a path to a file, or a file handler

    :returns: opened file ready to be iterated
    """
    textchars = bytearray({7,8,9,10,12,13,27} | set(range(0x20, 0x100)) - {0x7f})
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))

    if isinstance(filename, basestring) or isinstance(filename, basestring):
        fhandler = open(filename, 'rb')
        inputpath = True
        if tarfile.is_tarfile(filename):
            print('tar')
            thandler = tarfile.open(filename)
            if len(thandler.members) != 1:
                raise NotImplementedError(
                    'Not exactly one file in this tar archieve.')
            return magic_open(thandler.extractfile(thandler.getnames()[0]))
    else:
        fhandler = filename
        filename = fhandler.name
        inputpath = False
        start_of_file = ''
    if filename.endswith('.dsrc'):
        dsrc_binary = which('dsrc')
        if not dsrc_binary:
            raise Exception('\n\nERROR: DSRC binary not found, install it from:'
                            '\nhttps://github.com/lrog/dsrc/releases')
        proc = Popen([dsrc_binary, 'd', '-t%d' % (cpus or cpu_count()),
                      '-s', filename], stdout=PIPE, universal_newlines=True)
        return proc.stdout
    if inputpath:
        start_of_file = fhandler.read(1024)
        fhandler.seek(0)
        if is_binary_string(start_of_file):
            if start_of_file.startswith(b'\x50\x4b\x03\x04'):
                if verbose:
                    print('zip')
                zhandler = TextIOWrapper(zipfile.ZipFile(fhandler))
                if len(zhandler.NameToInfo) != 1:
                    raise NotImplementedError(
                        'Not exactly one file in this zip archieve.')
                return TextIOWrapper(BufferedReader(zhandler.open(list(zhandler.NameToInfo.keys())[0])))
            if is_binary_string(start_of_file) and start_of_file.startswith(b'\x42\x5a\x68'):
                if verbose:
                    print('bz2')
                fhandler.close()
                return TextIOWrapper(BufferedReader(bz2.BZ2File(filename)))
            if is_binary_string(start_of_file) and start_of_file.startswith(b'\x1f\x8b\x08'):
                if verbose:
                    print('gz')
                return TextIOWrapper(BufferedReader(gzip.GzipFile(fileobj=fhandler)))
        else:
            if verbose:
                print('text')
            fhandler.close()
            fhandler = open(filename, 'r')
    return fhandler


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def quality_plot(fnam, r_enz=None, nreads=float('inf'), axe=None,
                 output_prefix=None, extension_list=["png", "svg"], paired=False):
    """
    Plots the sequencing quality of a given FASTQ file. If a restrinction enzyme
    (RE) name is provided, can also represent the distribution of digested and
    undigested RE sites and estimate an expected proportion of dangling-ends.

    Proportion of dangling-ends is inferred by counting the number of times a
    dangling-end site, is found at the beginning of any of the reads (divided by
    the number of reads).

    :param fnam: path to FASTQ file
    :param None nreads: max number of reads to read, not necesary to read all
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    :param False paired: is input FASTQ contains both ends

    :returns: the percentage of dangling-ends (sensu stricto) and the percentage of
       reads with at least a ligation site.
    """
    phred = dict([(c, i) for i, c in enumerate(
        '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~')])
    if isinstance(r_enz, list):
        r_enzs = r_enz
    elif isinstance(r_enz, basestring):
        r_enzs = [r_enz]
    for k in list(RESTRICTION_ENZYMES.keys()):
        for i in range(len(r_enzs)):
            if k.lower() == str(r_enz[i]).lower():
                r_enz[i] = k
    # else let it as None

    quals = []
    henes = []
    sites = {}
    fixes = {}
    liges = OrderedDict()
    ligep = {}
    tkw = dict(size=4, width=1.5)
    fhandler = magic_open(fnam)
    if len(r_enzs) == 1 and r_enzs[0] is None:
        if nreads:
            while True:
                try:
                    next(fhandler)
                except EOFError:
                    break
                seq = next(fhandler)
                if 'N' in seq:
                    henes.extend([i for i, s in enumerate(seq) if s == 'N'])
                next(fhandler)
                line = next(fhandler)
                quals.append([phred[i] for i in line.strip()])
                if len(quals) > nreads:
                    break
        else: # do this because it's faster
            while True:
                try:
                    next(fhandler)
                except EOFError:
                    break
                seq = next(fhandler)
                if 'N' in seq:
                    henes.extend([i for i, s in enumerate(seq) if s == 'N'])
                next(fhandler)
                line = next(fhandler)
                quals.append([phred[i] for i in line.strip()])
    else:
        r_sites = {}
        d_sites = {}
        for r_enz in r_enzs:
            r_sites[r_enz] = RESTRICTION_ENZYMES[r_enz].replace('|', '')
            d_sites[r_enz] = repaired(r_enz)
            sites[r_enz] = []  # initialize dico to store undigested sites
            fixes[r_enz] = []  # initialize dico to store digested sites
        l_sites = religateds(r_enzs)
        l_sites = OrderedDict((k,iupac2regex(l_sites[k])) for k in l_sites)
        site = {}
        fixe = {}
        for r_enz in r_enzs:
            site[r_enz] = re.compile(iupac2regex(r_sites[r_enz]))
            fixe[r_enz] = re.compile(iupac2regex(d_sites[r_enz]))
        # ligation sites should appear in lower case in the sequence
        lige = {}
        for k in l_sites:
            liges[k] = []  # initialize dico to store sites
            ligep[k] = 0   # initialize dico to store sites
            l_sites[k] = l_sites[k].lower()
            lige[k] = re.compile(l_sites[k])
        callback = lambda pat: pat.group(0).lower()
        while len(quals) <= nreads:
            try:
                next(fhandler)
            except StopIteration:
                break
            seq = next(fhandler)
            # ligation sites replaced by lower case to ease the search
            for lig in list(l_sites.values()):
                seq = re.sub(lig.upper(), callback, seq)
            for r_enz in r_enzs:
                sites[r_enz].extend([m.start() for m in site[r_enz].finditer(seq)])
                # TODO: you cannot have a repaired/fixed site in the middle of
                # the sequence, this could be only checked at the beginning
                fixes[r_enz].extend([m.start() for m in fixe[r_enz].finditer(seq)])
            for k in lige:  # for each paired of cut-site
                liges[k].extend([m.start() for m in lige[k].finditer(seq)])
                if lige[k].search(seq):
                    ligep[k] += 1
            # store the number of Ns found in the sequences
            if 'N' in seq:
                henes.extend([i for i, s in enumerate(seq) if s == 'N'])
            next(fhandler)
            line = next(fhandler)
            quals.append([phred[i] for i in line.strip()])
    fhandler.close()
    if not nreads:
        nreads = len(quals)
    quals = zip_longest(*quals, fillvalue=float('nan'))
    meanquals, errorquals = list(zip(*[(nanmean(q), nanstd(q)) for q in quals]))
    max_seq_len = len(meanquals)

    if axe:
        ax = axe
        fig = axe.get_figure()
        ax2 = fig.add_subplot(212)
    else:  # configure plot
        if len(r_enzs) == 1 and r_enzs[0] is None:  # do both plots
            _, ax = plt.subplots(1, 1, figsize=(15, 6))
        else:  # only do the quality_plot plot
            _, (ax, ax2) = plt.subplots(2, 1, figsize=(15, 12))
        ax.patch.set_facecolor('lightgrey')
        ax.patch.set_alpha(0.4)
        ax.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
        ax.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
        ax.set_axisbelow(True)
        # remove tick marks
        ax.tick_params(axis='both', direction='out', top=False, right=False,
                       left=False, bottom=False)
        ax.tick_params(axis='both', direction='out', top=False, right=False,
                       left=False, bottom=False, which='minor')

    ax.errorbar(list(range(max_seq_len)), meanquals,
                linewidth=1, elinewidth=1, color='darkblue',
                yerr=errorquals, ecolor='orange')

    ax.set_xlim((0, max_seq_len))
    ax.set_xlabel('Nucleotidic position')
    ax.set_ylabel('PHRED score')
    ax.set_title('Sequencing Quality (%d reads)' % (nreads))
    ax.yaxis.label.set_color('darkblue')
    ax.tick_params(axis='y', colors='darkblue', **tkw)
    axb = ax.twinx()
    # quality_plot plot
    axb.plot([henes.count(i) for i in range(max_seq_len)], linewidth=1,
             color='black', linestyle='--')
    axb.yaxis.label.set_color('black')
    axb.tick_params(axis='y', colors='black', **tkw)
    axb.set_ylabel('Number of "N" per position')
    try: # no Ns found (yes... it happens)
        axb.set_yscale('log')
        with catch_warnings():
            simplefilter("ignore")
            axb.set_ylim((0, axb.get_ylim()[1] * 1000))
    except ValueError:
        axb.set_yscale('linear')
    ax.set_ylim((0, ax.get_ylim()[1]))
    ax.set_xlim((0, max_seq_len))

    # Hi-C plot
    percent_digested_sites = {}
    if not (len(r_enzs) == 1 and r_enzs[0] is None):
        ax.set_title('Sequencing Quality and deconvolution (%s %d reads)' % (
            ', '.join(map(str, r_enzs)), nreads))
        ax.set_xlabel('')
        plt.setp(ax.get_xticklabels(), visible=False)
        ax2.patch.set_facecolor('lightgrey')
        ax2.patch.set_alpha(0.4)
        ax2.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
        ax2.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
        ax2.set_axisbelow(True)
        ax2.set_xlabel('Nucleotidic position')

        # seq_len is the length of the line to plot. we don't want to plot
        # if there is no room for the cut-site, or ligation site.
        site_len = max((max([len(r_sites[k]) for k in r_sites]),
                        max([len(l_sites[k]) for k in l_sites]),
                        max([len(d_sites[k]) for k in d_sites])))
        seq_len = max_seq_len - site_len

        # transform dictionaries of positions into dictionaries of counts
        for r_enz in sites:
            sites[r_enz] = [sites[r_enz].count(k) for k in range(seq_len)] # Undigested
            fixes[r_enz] = [fixes[r_enz].count(k) for k in range(seq_len)] # DE
        for r1, r2 in liges:
            liges[(r1, r2)] = [liges[(r1, r2)].count(k) for k in range(seq_len)] # OK

        # in case the pattern of the repaired cut-site contains the target
        # cut-site pattern. These sites were counted twice, once in the
        # undigested, and once in the repaired. We remove them from the
        # repaired:
        for r_enz in r_enzs:
            if d_sites[r_enz] in r_sites[r_enz]:
                pos = r_sites[r_enz].find(d_sites[r_enz])

                fixes[r_enz] = (fixes[r_enz][:pos] +
                                [fixes[r_enz][k] - sites[r_enz][k-pos]
                                 for k in range(pos, seq_len)])
        # same for ligated sites
        for r_enz1 in r_enzs:
            for r_enz2 in r_enzs:
                if d_sites[r_enz1] not in l_sites[(r_enz1, r_enz2)]:
                    continue
                pos = l_sites[(r_enz1, r_enz2)].find(d_sites[r_enz1])
                fixes[r_enz1] = (fixes[r_enz1][:pos] +
                                 [fixes[r_enz1][k] - liges[(r_enz1, r_enz2)][k - pos]
                                  for k in range(pos, seq_len)])

        # remove anything that could be in between the two read ends
        if paired:
            for k in sites:
                sites[k][max_seq_len // 2 - site_len:
                         max_seq_len // 2] = [float('nan')] * site_len
                fixes[k][max_seq_len // 2 - site_len:
                         max_seq_len // 2] = [float('nan')] * site_len
            for k in liges:
                liges[k][max_seq_len // 2 - site_len:
                         max_seq_len // 2] = [float('nan')] * site_len

        # plot undigested cut-sites
        color = iter(plt.cm.Reds(linspace(0.3, 0.95, len(r_enzs))))
        for r_enz in sites:
            # print 'undigested', r_enz
            # print sites[r_enz][:20]
            ax2.plot(sites[r_enz], linewidth=2, color=next(color),
                     alpha=0.9,
                     label='Undigested RE site (%s: %s)' % (r_enz, r_sites[r_enz])
                     if any([f > 0 for f in fixes[r_enz]])
                     else 'Undigested & Dangling-Ends (%s: %s)' % (r_enz, r_sites[r_enz]))
        ax2.set_ylabel('Undigested')
        ax2.yaxis.label.set_color('darkred')
        ax2.tick_params(axis='y', colors='darkred', **tkw)

        lines, labels = ax2.get_legend_handles_labels()

        ax3 = ax2.twinx()
        color = iter(plt.cm.Blues(linspace(0.3, 0.95, len(liges))))
        for r1, r2 in liges:
            # print 'ligated', r1, r2
            # print liges[(r1, r2)][:20]
            ax3.plot(liges[(r1, r2)], linewidth=2, color=next(color),
                     alpha=0.9,
                     label = 'Ligated (%s-%s: %s)' % (r1, r2, l_sites[(r1, r2)].upper()))
        ax3.yaxis.label.set_color('darkblue')
        ax3.tick_params(axis='y', colors='darkblue', **tkw)
        ax3.set_ylabel('Ligated')

        tmp_lines, tmp_labels = ax3.get_legend_handles_labels()
        lines.extend(tmp_lines)
        labels.extend(tmp_labels)

        color = iter(plt.cm.Greens(linspace(0.3, 0.95, len(r_enzs))))
        for i, r_enz in enumerate(r_enzs):
            if any([f > 0 for f in fixes[r_enz]]):
                ax4 = ax2.twinx()
                ax4.spines["right"].set_position(("axes", 1.07))
                make_patch_spines_invisible(ax4)
                ax4.spines["right"].set_visible(True)
                # print 'repaired', r_enz
                # print fixes[r_enz][:20]
                ax4.plot(fixes[r_enz], linewidth=2, color=next(color),
                         alpha=0.9,
                         label='Dangling-ends (%s: %s)' % (r_enz, d_sites[r_enz]))
                ax4.yaxis.label.set_color('darkgreen')
                ax4.tick_params(axis='y', colors='darkgreen', **tkw)
                ax4.set_ylabel('Dangling-ends')
                tmp_lines, tmp_labels = ax4.get_legend_handles_labels()
                lines.extend(tmp_lines)
                labels.extend(tmp_labels)
            else:
                ax2.set_ylabel('Undigested & Dangling-ends')
        ax2.set_xlim((0, max_seq_len))

        # Count ligation sites
        lig_cnt = {}
        for k in liges:
            lig_cnt[k] = (nansum(liges[k]) - liges[k][0] -
                              liges[k][max_seq_len // 2])

        # Count undigested sites
        sit_cnt = {}
        for r_enz in r_enzs:
            sit_cnt[r_enz] = (nansum(sites[r_enz]) - sites[r_enz][0] -
                              sites[r_enz][max_seq_len // 2])

        # Count Dangling-Ends
        des = {}
        for r_enz in r_enzs:
            if any([f > 0 for f in fixes[r_enz]]):
                des[r_enz] = ((100. * (fixes[r_enz][0] + (fixes[r_enz][(max_seq_len // 2)]
                                                          if paired else 0))) / nreads)
            else:
                des[r_enz] = (100. * (sites[r_enz][0] + (sites[r_enz][(max_seq_len // 2)]
                                                         if paired else 0))) / nreads

        # Decorate plot
        title = ''

        for r_enz in r_enzs:
            lcnt = float(sum([lig_cnt[(r_enz1, r_enz2)] * (2 if r_enz1 == r_enz2 else 1)
                              for r_enz1 in r_enzs for r_enz2 in r_enzs
                              if r_enz1 == r_enz or r_enz2 == r_enz]))
            title += ('Percentage of digested sites (not considering Dangling-Ends) '
                      '%s: %.1f%%\n' % (r_enz,
                                        100. * float(lcnt) / (lcnt + sit_cnt[r_enz])))
            percent_digested_sites[r_enz] = 100. * float(lcnt) / (lcnt + sit_cnt[r_enz])
        for r_enz in r_enzs:
            title += 'Percentage of dangling-ends %s: %.1f%%\n' % (r_enz, des[r_enz])

        for r_enz1 in r_enzs:
            for r_enz2 in r_enzs:
                title += ('Percentage of reads with ligation site (%s-%s): %.1f%% \n' %
                          (r_enz1, r_enz2, (ligep[(r_enz1, r_enz2)] * 100.) / nreads))
        plt.title(title.strip(), size=10, ha='left', x=0)
        plt.subplots_adjust(right=0.85)
        ax2.legend(lines, labels, bbox_to_anchor=(0.75, 1.0),
                   loc=3, borderaxespad=0., frameon=False, fontsize=9)
    plt.tight_layout()
    if output_prefix:
        for ext in extension_list:
            plt.savefig(output_prefix + "." + ext, )
        plt.close('all')

    elif not axe:
        plt.show()
    for k in ligep:
        ligep[k] = (ligep[k] * 100.) / nreads
    if len(r_enzs) == 1 and r_enzs[0] is None:
        return {}, {}
    return des, ligep, percent_digested_sites


parser = argparse.ArgumentParser()


parser.add_argument("-f", "--input_fastq", action="store", dest="input_fastq", required=True,
                    help="Input fastq file, can be compressed. Required")
parser.add_argument("-e", "--restriction_enzyme_list", action="store", dest="restriction_enzyme_list", default=None,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of the restriction enzymes.")
parser.add_argument("-n", "--read_number", action="store", dest="read_number", default=100000, type=int,
                    help="Number of reads to analyze. Default: 100000")
parser.add_argument("-s", "--percentage_header_prefix", action="store", dest="percentage_header_prefix", default="",
                    help="Prefix of percentage column in the header. "
                         "Useful if you with to analyze multiple files and combine stats in a single table later. "
                         "Default: '', i.e not set")
parser.add_argument("-i", "--interleaved", action="store_true", dest="interleaved", default=False,
                    help="Input fastq is interleaved, i.e contains both forward and reverse reads. Default: not set")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files.")
parser.add_argument("-x", "--figure_format_list", action="store", dest="figure_format_list", default=["png", "svg"],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of extensions for output figures. Default: png,svg")
args = parser.parse_args()

des, ligep, percent_digested_sites = quality_plot(args.input_fastq, r_enz=args.restriction_enzyme_list,
                                                  nreads=args.read_number,
                                                  axe=None, output_prefix=args.output_prefix,
                                                  extension_list=args.figure_format_list, paired=args.interleaved)

with open(args.output_prefix + ".stats", "w") as out_fd:
    out_fd.write("#type\trestrictase(s)\t{0}percentage\n".format(args.percentage_header_prefix))
    for restrictase in percent_digested_sites:
        out_fd.write("{0}\t{1}\t{2:3.2f}\n".format("digested sites", restrictase, percent_digested_sites[restrictase]))
    for restrictase in des:
        out_fd.write("{0}\t{1}\t{2:3.2f}\n".format("dangling ends", restrictase, des[restrictase]))
    for restrictase_couple in ligep:
        out_fd.write("{0}\t{1}-{2}\t{3:3.2f}\n".format("ligation_sites", restrictase_couple[0],
                                                  restrictase_couple[1], ligep[restrictase_couple]))
