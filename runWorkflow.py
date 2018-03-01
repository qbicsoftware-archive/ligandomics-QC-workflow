from CTDopts.CTDopts import _InFile, CTDModel, args_from_file
import sys
import os
import subprocess
import re

pattern = re.compile('Q\w{4}[0-9]{3}[a-zA-Z]\w')

wf_dir = sys.argv[1]
ctd_params = args_from_file(wf_dir + '/WORKFLOW-CTD')
ctd_files = args_from_file(wf_dir + '/IN-FILESTOSTAGE')
r_scripts = '/lustre_cfc/software/qbic/openms/openMS2.0-immuno-src/share/OpenMS/SCRIPTS/{s}'
fasta_decoy_path = '/lustre_cfc/qbic/reference_genomes/Homo_sapiens/Proteome/swissprotHUMANwoi_130927.revCat.fasta'

data_path = '%s/data/' % wf_dir
result_path = '%s/result/' % wf_dir
log_path = '%s/logs/' % wf_dir

mzmlFiles = []

for filePath in ctd_files['Mass Spectrometry Data']:
	fileName = filePath.split('/')[-1]
	mzmlFiles.append('%s%s' % (data_path, fileName))

allele_fileName = ctd_files['MHC Alleles'].split('/')[-1]
allelePath = '%s%s' % (data_path, allele_fileName)

logfilename = 'ligandomicsQC_1_0_workflow.logs'
logfile = open(logfilename, 'w')

param_file = '%s/src/comet_MHCI_2015024.params' % wf_dir

# QC
# mzml processing
for mzml in mzmlFiles:
	if mzml.endswith('.gz'):
		logfile.write("Extracting gzipped content... \n")
		cmd = "gzip -d {f}".format(f=mzml)
		os.system(cmd)
		mzml = mzml.replace('.gz', '')

	idPath = mzml.replace('mzML', 'idXML')

	identifier = mzml.split('/')[-1].split('.')[0]

	if ctd_params['centroided'] == 'false':
		pickpeakcommand = 'PeakPickerHiRes -in {i} -out {o} -algorithm:ms_levels 1'
		subprocess.call(pickpeakcommand.format(i=mzml, o=mzml).split(),stderr=logfile, stdout=logfile)

	ffcom = 'FileFilter -in {i} -out {o}  -peak_options:remove_chromatograms -peak_options:indexed_file true'
	subprocess.call(ffcom.format(i=mzml, o=mzml).split(),stderr=logfile, stdout=logfile)

	commandComet = 'comet.2015024.linux.exe -P{p} {s}'.format(p=param_file, s=mzml) 
	subprocess.call(commandComet.split(),stderr=logfile, stdout=logfile)
	 
	commandFileConverter = 'IDFileConverter -in {f} -out {o} {a}'

	annotation = '-mz_file {f} -add_ionmatch_annotation'
	subprocess.call(commandFileConverter.format(f=mzml.replace('mzML','pep.xml'), o=idPath, a='').split(),stderr=logfile, stdout=logfile)
	subprocess.call(commandFileConverter.format(f=idPath, o=idPath, a=annotation.format(f=mzml)).split(),stderr=logfile, stdout=logfile)

	peptideIndexer = 'PeptideIndexer -in {f} -out {o} -fasta {d} -decoy_string XXX -prefix -enzyme:specificity none'.format(f=idPath, o=idPath, d=fasta_decoy_path)
	subprocess.call(peptideIndexer.split(),stderr=logfile, stdout=logfile)

	### predict hits of fitting length and calc FDR and PEP
	#FDR calc

	falseDiscovery = 'FalseDiscoveryRate -in {f} -out {o} -algorithm:add_decoy_peptides'.format(f=idPath,o=idPath)
	subprocess.call(falseDiscovery.split(),stderr=logfile, stdout=logfile)

	postErrorProb = 'IDPosteriorErrorProbability -in {f} -out {o}'.format(f=idPath, o=idPath)
	subprocess.call(postErrorProb.split(),stderr=logfile, stdout=logfile)

	#predict
	peptideHit = 'python peptidehit_predict_netMHC.py -c {c} -in {f} -out {o} -allele {a}'.format(f=idPath,o=idPath,a=allelePath,c=ctd_params['MHC class'])
	subprocess.call(peptideHit.split(),stderr=logfile, stdout=logfile)

	#QC creation
	imageCreator = 'ImageCreator -in {f} -out {o} -width 640 -height 640'.format(f=mzml,o=mzml.replace('mzML', 'png'))
	subprocess.call(imageCreator.split(),stderr=logfile, stdout=logfile)

	#Feature Finding
	featureFinder = 'FeatureFinderCentroided -in {f} -out {o} -ini {i}'.format(f=mzml,o=mzml.replace('mzML', 'featureXML'),i='FeatureFinderCentroided.ini')
	subprocess.call(featureFinder.split(),stderr=logfile, stdout=logfile)

	#fork FDR filtered file
	filterfile = 'IDFilter -in {f} -out {o} -score:pep 0.05 -remove_decoys'
	subprocess.call(filterfile.format(f=idPath,o=idPath.replace('idXML', 'filtered.idXML')).split(),stderr=logfile, stdout=logfile)

	#QCCalculator
	qcresult = os.path.join(result_path, mzml.replace('mzML', 'qcML').split('/')[-1])
	qccalc = 'QCCalculator -in {f} -feature {ff} -id {i} -out {o}'.format(f=mzml, ff=mzml.replace('mzML', 'featureXML'), i=idPath.replace('idXML', 'filtered.idXML'), o=qcresult)
	subprocess.call(qccalc.split(),stderr=logfile, stdout=logfile)

	qcembedd = 'QCEmbedder -in {f} -out {o} -plot {p} -qp_att_acc {qp} -cv_acc {cv}'

	#QCEmbedder
	subprocess.call(qcembedd.format(f=qcresult, o=qcresult, p=mzml.replace('mzML', 'png'), qp='QC:0000004', cv='QC:0000055').split(),stderr=logfile, stdout=logfile)

	qextract = 'QCExtractor -qp {qp} -out_csv {o} -in {i}'
	subprocess.call(qextract.format(qp='QC:0000022' , o=data_path + identifier + '_qc22.csv', i=qcresult).split(), stderr=logfile, stdout=logfile)

	r_command = 'Rscript --vanilla {s} {i} {o}'
	subprocess.call(r_command.format(s=r_scripts.format(s='ProduceQCFigures_tic.R'), i=data_path + identifier + '_qc22.csv', o=data_path + identifier + 'tic.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(qcembedd.format(qp='QC:0000023', cv='MS:1000235', f=qcresult, o=qcresult, p=data_path + identifier + '_tic.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(qextract.format(qp='QC:0000038', o=data_path + identifier + '_qc38.csv', i=qcresult).split(),stderr=logfile, stdout=logfile)

	subprocess.call(r_command.format(s=r_scripts.format(s='ProduceQCFigures_acc.R'), i=data_path + identifier + '_qc38.csv', o=data_path +  identifier +'_acc.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(qcembedd.format(qp='QC:0000041', cv='QC:0000053', f=qcresult, o=qcresult, p=data_path + identifier + '_acc.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(r_command.format(s=r_scripts.format(s='ProduceQCFigures_rt_acc.R'), i=data_path + identifier +'_qc38.csv', o=data_path + identifier + '_acc_rt.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(qcembedd.format(qp='QC:0000041', cv='QC:0000054', f=qcresult, o=qcresult, p=data_path + identifier + '_acc_rt.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(qextract.format(qp='QC:0000018', o=data_path + identifier + '_qc18.csv', i=qcresult).split(), stderr=logfile, stdout=logfile)

	subprocess.call(r_command.format(s=r_scripts.format(s='ProduceQCFigures_inj.R'), i=data_path + identifier + '_qc18.csv', o=data_path + identifier + '_inj.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(qcembedd.format(qp='QC:0000023', cv='QC:0000051', f=qcresult, o=qcresult, p=data_path + identifier + '_inj.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(qextract.format(qp='QC:0000044', o=data_path + identifier + '_qc44.csv', i=qcresult).split(), stderr=logfile, stdout=logfile)

	subprocess.call(r_command.format(s=r_scripts.format(s='ProduceQCFigures_idmap.R'), i=data_path + identifier + '_qc44.csv' + ' ' + data_path + identifier + '_qc38.csv', o=data_path + identifier + '_idmap.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(qcembedd.format(qp='QC:0000035', cv='QC:0000052', f=qcresult, o=qcresult, p=data_path + identifier + '_idmap.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call('python SEperformance.py -in {i} -out {o}'.format(i=idPath, o=data_path + identifier + '_performance.tsv').split(), stderr=logfile, stdout=logfile)

	subprocess.call(r_command.format(s='ProduceQCFigures_allelecounts.R', i=data_path + identifier + '_performance.tsv', o=data_path + identifier + '_allelecounts.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(qcembedd.format(qp='QC:0000023', cv='QC:0000051', f=qcresult, o=qcresult, p=data_path + identifier + '_allelecounts.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(r_command.format(s='ProduceQCFigures_bindercounts.R', i=data_path + identifier + '_performance.tsv', o=data_path + identifier + '_bindercounts.png').split(), stderr=logfile, stdout=logfile)

	subprocess.call(qcembedd.format(qp='QC:0000023', cv='QC:0000051', f=qcresult, o=qcresult, p=data_path + identifier + '_bindercounts.png').split(), stderr=logfile, stdout=logfile)

	r_command = 'Rscript --vanilla {s} {i} {o} {l}'
	subprocess.call(r_command.format(s='FDRoptimizer.R', i=data_path + identifier + '_performance.tsv', o=data_path + identifier + '_fdroptimizer.png', l=data_path + identifier + '_temp.csv').split(), stderr=logfile, stdout=logfile)

	subprocess.call(qcembedd.format(qp='QC:0000023', cv='QC:0000051', f=qcresult, o=qcresult, p=data_path + identifier + '_fdroptimizer.png').split(), stderr=logfile, stdout=logfile)

logfile.close()
subprocess.call(["mv", logfilename, log_path])
