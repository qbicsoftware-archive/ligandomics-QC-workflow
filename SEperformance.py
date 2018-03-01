#!/usr/bin/env python
__author__ = 'walzer'
import sys
import argparse
import pyopenms as oms
import pandas as pd
from os import path

VERSION = "0.2"


def __main__():
    parser = argparse.ArgumentParser(version=VERSION)
    parser.add_argument('-in', dest="inf", help='<Required> full path to the input idXML', required=True)
    parser.add_argument('-out', dest="out", help="<Required> full path to the csv to be written", required=True)
    parser.add_argument('-doi', '--dept_of_immuno-filenames', dest='doi', action='store_true', help="If filename conforms with the filename scheme of the dept. of immunology... (in doubt, it probably doesn'T)")
    parser.set_defaults(doi=False)
    parser.add_argument('-qcML', dest="qcml", action='store_true', help="If instead of writing a csv embedding in a qcml file - prerequisite: -out must be a existing qcml file")
    parser.set_defaults(qcml=False)

    options = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    if not (options.inf or options.out):
        parser.print_help()
        sys.exit(1)

    pros = list()
    peps = list()
    f = oms.IdXMLFile()
    f.load(options.inf, pros, peps)

    if options.qcml:
        try:
            qf = oms.QcMLFile()
            qf.load(oms.String(options.out))
        except:
            print "no usable qcml file given"
            sys.exit(1)

    SE = "NA"
    SE = pros[0].getSearchEngine()
    sesmv = "NA"
    if SE == 'MS-GF+':
        sesmv = 'MS:1002053'
    elif SE == 'XTandem':
        sesmv = 'E-Value'
    elif SE == 'Comet':
        sesmv = 'MS:1002257'
    elif SE == 'Mascot':
        sesmv = 'EValue'
    if not sesmv:
        print "no adequate search engine score found"
        #sys.exit(1)

    organ, repli, patient = None, None, None
    if options.doi:
        organ = path.basename(options.inf).split('_')[3]
        repli = path.basename(options.inf).split('_')[6].split('#')[-1]
        patient = path.basename(options.inf).split('_')[2]
    if not organ or not repli or not patient:
        options.doi = False

    sepe = list()
    for pep in peps:
        pep.sort()
        hits = pep.getHits()
        st = pep.getScoreType()
        if st != 'Posterior Error Probability':
            print "Warning, no PEP as score type found!"
        if pep.metaValueExists('spectrum_reference'):
            sn_ref = pep.getMetaValue('spectrum_reference')
            sn = sn_ref.split('=')[-1]
        else:
            sn = '?'
            sn_ref = '?'
        for i, h in enumerate(hits):
            row = dict()
            row['file'] = path.basename(options.inf)
            row['sequence'] = h.getSequence().toUnmodifiedString()
            row['modified'] = h.getSequence().isModified()
            row['rank'] = i + 1
            row['spectrum'] = sn
            row['spectrum_reference'] = sn_ref # match with spectrum_native_id von featureXML and spectrum_reference from mapped PeptideIdentifications
            if pep.getScoreType() == 'Posterior Error Probability':
                row['PEP'] = h.getScore()
            else:
                row['PEP'] = "NA"
            row['SE'] = SE
            if h.metaValueExists(sesmv):
                row['SEscore'] = h.getMetaValue(sesmv)
            if st == 'q-value':
                row['qvalue'] = h.getScore()
            else:
                if h.metaValueExists('q-value_score'):
                    row['qvalue'] = h.getMetaValue('q-value_score')
            if h.metaValueExists('binder'):
                if h.metaValueExists('weak'):
                    row['binder'] = 'weak'
                elif h.metaValueExists('strong'):
                    row['binder'] = 'strong'
                else:
                    row['binder'] = 'yes'
                row['allele'] = h.getMetaValue("binder")
            else:
                row['binder'] = 'no'
            if "decoy" not in h.getMetaValue("target_decoy"):
                row['target'] = 'target'
            else:
                row['target'] = 'decoy'
            if options.doi:
                row['organ'] = organ
                row['replicate'] = repli
                row['patient'] = patient
            sepe.append(row)

    outframe = pd.DataFrame(sepe)

    if options.qcml:
        idxa = oms.Attachment()
        idxa.name = "generic table" #"extended id tab"
        idxa.cvRef = "qcML"
        idxa.cvAcc = "QC:0000049"
        idxa.qualityRef = "QC:0000025"

        idxa.colTypes = list(outframe.columns)
        rows = [[str(row[f]) for f in idxa.colTypes] for ind, row in outframe.iterrows()] #writes nan if cell not filled
        idxa.tableRows = rows

        ls = list()
        qf.getRunNames(ls)
        qf.addRunAttachment(ls[0], idxa)
        #qf.store(oms.String(options.out.replace('%', 'perc')))
        qf.store(oms.String(options.out))

    else:
        outframe.to_csv(options.out, sep='\t', index=False)

if __name__ == '__main__':
    __main__()