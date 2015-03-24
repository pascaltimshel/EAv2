#!/usr/bin/env python2.7

import os
import sys
import glob

import pdb

###################################### FILE SNIPPETS ######################################
################## file_annot ##################

### FORMAT ###
# This format is tricky! Be carefull. This file does not seem to behave well!

### EXAMPLE ###
# ^Annotation
# !Annotation_date = Nov 23 2007
# ..... MORE ....
# #GO:Function ID = Gene Ontology Function identifier
# #GO:Process ID = Gene Ontology Process identifier
# #GO:Component ID = Gene Ontology Component identifier
# !platform_table_begin
# ID	Gene title	Gene symbol	Gene ID	UniGene title	UniGene symbol	UniGene ID	Nucleotide Title	GI	GenBank Accession	Platform_CLONEID	Platform_ORF	Platform_SPOTID	Chromosome location	Chromosome annotation	GO:Function	GO:Process	GO:Component	GO:Function ID	GO:Process ID	GO:Component ID
# HEEBO-001-CTRL1A1	actin, beta	ACTB	60				Homo sapiens actin, beta (ACTB), mRNA	5016088	NM_001101				7p15-p12	Chromosome 7, NC_000007.12 (5533312..5536747, complement)	ATP binding///nucleotide binding///protein binding///protein binding///structural constituent of cytoskeleton///structural molecule activity	cell motility///sensory perception of sound	NuA4 histone acetyltransferase complex///cytoplasm///cytoskeleton///cytoskeleton///cytosol///soluble fraction	GO:0005524///GO:0000166///GO:0005515///GO:0005515///GO:0005200///GO:0005198	GO:0006928///GO:0007605	GO:0035267///GO:0005737///GO:0005856///GO:0005856///GO:0005829///GO:0005625
# HEEBO-001-CTRL1A2	actin, beta	ACTB	60				Homo sapiens actin, beta (ACTB), mRNA	5016088	NM_001101				7p15-p12	Chromosome 7, NC_000007.12 (5533312..5536747, complement)	ATP binding///nucleotide binding///protein binding///protein binding///structural constituent of cytoskeleton///structural molecule activity	cell motility///sensory perception of sound	NuA4 histone acetyltransferase complex///cytoplasm///cytoskeleton///cytoskeleton///cytosol///soluble fraction	GO:0005524///GO:0000166///GO:0005515///GO:0005515///GO:0005200///GO:0005198	GO:0006928///GO:0007605	GO:0035267///GO:0005737///GO:0005856///GO:0005856///GO:0005829///GO:0005625
# HEEBO-001-CTRL1A3	actin, beta	ACTB	60				Homo sapiens actin, beta (ACTB), mRNA	5016088	NM_001101				7p15-p12	Chromosome 7, NC_000007.12 (5533312..5536747, complement)	ATP binding///nucleotide binding///protein binding///protein binding///structural constituent of cytoskeleton///structural molecule activity	cell motility///sensory perception of sound	NuA4 histone acetyltransferase complex///cytoplasm///cytoskeleton///cytoskeleton///cytosol///soluble fraction	GO:0005524///GO:0000166///GO:0005515///GO:0005515///GO:0005200///GO:0005198	GO:0006928///GO:0007605	GO:0035267///GO:0005737///GO:0005856///GO:0005856///GO:0005829///GO:0005625
# HEEBO-001-CTRL1A4	actin, beta	ACTB	60				Homo sapiens actin, beta (ACTB), mRNA	5016088	NM_001101				7p15-p12	Chromosome 7, NC_000007.12 (5533312..5536747, complement)	ATP binding///nucleotide binding///protein binding///protein binding///structural constituent of cytoskeleton///structural molecule activity	cell motility///sensory perception of sound	NuA4 histone acetyltransferase complex///cytoplasm///cytoskeleton///cytoskeleton///cytosol///soluble fraction	GO:0005524///GO:0000166///GO:0005515///GO:0005515///GO:0005200///GO:0005198	GO:0006928///GO:0007605	GO:0035267///GO:0005737///GO:0005856///GO:0005856///GO:0005829///GO:0005625
# HEEBO-001-CTRL1A5	actin, beta	ACTB	60				Homo sapiens actin, beta (ACTB), mRNA	5016088	NM_001101				7p15-p12	Chromosome 7, NC_000007.12 (5533312..5536747, complement)	ATP binding///nucleotide binding///protein binding///protein binding///structural constituent of cytoskeleton///structural molecule activity	cell motility///sensory perception of sound	NuA4 histone acetyltransferase complex///cytoplasm///cytoskeleton///cytoskeleton///cytosol///soluble fraction	GO:0005524///GO:0000166///GO:0005515///GO:0005515///GO:0005200///GO:0005198	GO:0006928///GO:0007605	GO:0035267///GO:0005737///GO:0005856///GO:0005856///GO:0005829///GO:0005625
# HEEBO-001-CTRL1A6	actin, beta	ACTB	60				Homo sapiens actin, beta (ACTB), mRNA	5016088	NM_001101				7p15-p12	Chromosome 7, NC_000007.12 (5533312..5536747, complement)	ATP binding///nucleotide binding///protein binding///protein binding///structural constituent of cytoskeleton///structural molecule activity	cell motility///sensory perception of sound	NuA4 histone acetyltransferase complex///cytoplasm///cytoskeleton///cytoskeleton///cytosol///soluble fraction	GO:0005524///GO:0000166///GO:0005515///GO:0005515///GO:0005200///GO:0005198	GO:0006928///GO:0007605	GO:0035267///GO:0005737///GO:0005856///GO:0005856///GO:0005829///GO:0005625
# ..... MORE ....
# HEEBO-004-CTRL4O11	retinoblastoma 1 (including osteosarcoma)	RB1	5925				Homo sapiens retinoblastoma 1 (including osteosarcoma) (RB1), mRNA	108773786	NM_000321				13q14.2	Chromosome 13, NC_000013.9 (47775912..47954023)	androgen receptor binding///kinase binding///molecular_function///transcription coactivator activity///transcription factor activity///transcription factor binding///transcription repressor activity	G1 phase///G1/S transition of mitotic cell cycle///M phase///androgen receptor signaling pathway///cell cycle///cell cycle arrest///cell cycle checkpoint///cell division///chromatin modification///enucleate erythrocyte differentiation///negative regulation of cell growth///negative regulation of cell proliferation///negative regulation of protein kinase activity///negative regulation of transcription from RNA polymerase II promoter///positive regulation of macrophage differentiation///positive regulation of transcription from RNA polymerase II promoter///regulation of lipid kinase activity///regulation of transcription, DNA-dependent///striated muscle cell differentiation///transcription	PML body///chromatin///nucleus///nucleus///spindle///transcription factor complex	GO:0050681///GO:0019900///GO:0003674///GO:0003713///GO:0003700///GO:0008134///GO:0016564	GO:0051318///GO:0000082///GO:0000279///GO:0030521///GO:0007049///GO:0007050///GO:0000075///GO:0051301///GO:0016568///GO:0043353///GO:0030308///GO:0008285///GO:0006469///GO:0000122///GO:0045651///GO:0045944///GO:0043550///GO:0006355///GO:0051146///GO:0006350	GO:0016605///GO:0000785///GO:0005634///GO:0005634///GO:0005819///GO:0005667
# HEEBO-004-CTRL4O12	retinoblastoma 1 (including osteosarcoma)	RB1	5925				Homo sapiens retinoblastoma 1 (including osteosarcoma) (RB1), mRNA	108773786	NM_000321				13q14.2	Chromosome 13, NC_000013.9 (47775912..47954023)	androgen receptor binding///kinase binding///molecular_function///transcription coactivator activity///transcription factor activity///transcription factor binding///transcription repressor activity	G1 phase///G1/S transition of mitotic cell cycle///M phase///androgen receptor signaling pathway///cell cycle///cell cycle arrest///cell cycle checkpoint///cell division///chromatin modification///enucleate erythrocyte differentiation///negative regulation of cell growth///negative regulation of cell proliferation///negative regulation of protein kinase activity///negative regulation of transcription from RNA polymerase II promoter///positive regulation of macrophage differentiation///positive regulation of transcription from RNA polymerase II promoter///regulation of lipid kinase activity///regulation of transcription, DNA-dependent///striated muscle cell differentiation///transcription	PML body///chromatin///nucleus///nucleus///spindle///transcription factor complex	GO:0050681///GO:0019900///GO:0003674///GO:0003713///GO:0003700///GO:0008134///GO:0016564	GO:0051318///GO:0000082///GO:0000279///GO:0030521///GO:0007049///GO:0007050///GO:0000075///GO:0051301///GO:0016568///GO:0043353///GO:0030308///GO:0008285///GO:0006469///GO:0000122///GO:0045651///GO:0045944///GO:0043550///GO:0006355///GO:0051146///GO:0006350	GO:0016605///GO:0000785///GO:0005634///GO:0005634///GO:0005819///GO:0005667
# HEEBO-004-CTRL4O13												secretory protein LOC284013								
# HEEBO-004-CTRL4O14												LOC440859								
# HEEBO-004-CTRL4O15												LOC441859								
# HEEBO-004-CTRL4O24												forkhead box D3								
# HEEBO-004-CTRL4P1	protein tyrosine phosphatase, receptor type, C	PTPRC	5788				Homo sapiens protein tyrosine phosphatase, receptor type, C (PTPRC), transcript variant 1, mRNA	115385975	NM_002838				1q31-q32	Chromosome 1, NC_000001.9 (196874848..196992953)	hydrolase activity///phosphoric monoester hydrolase activity///protein binding///protein kinase binding///protein kinase binding///protein tyrosine phosphatase activity///protein tyrosine phosphatase activity///receptor activity///transmembrane receptor protein tyrosine phosphatase activity	B cell proliferation///B cell receptor signaling pathway///T cell differentiation///T cell receptor signaling pathway///T cell receptor signaling pathway///cell surface receptor linked signal transduction///defense response to virus///immunoglobulin biosynthetic process///negative regulation of T cell mediated cytotoxicity///negative regulation of cytokine and chemokine mediated signaling pathway///negative regulation of protein kinase activity///negative regulation of protein kinase activity///positive regulation of B cell proliferation///positive regulation of T cell proliferation///positive regulation of antigen receptor-mediated signaling pathway///positive regulation of protein kinase activity///protein amino acid dephosphorylation///protein amino acid dephosphorylation///regulation of cell cycle///regulation of progression through S phase///release of sequestered calcium ion into cytosol	focal adhesion///integral to plasma membrane///membrane	GO:0016787///GO:0016791///GO:0005515///GO:0019901///GO:0019901///GO:0004725///GO:0004725///GO:0004872///GO:0005001	GO:0042100///GO:0050853///GO:0030217///GO:0050852///GO:0050852///GO:0007166///GO:0051607///GO:0002378///GO:0001915///GO:0001960///GO:0006469///GO:0006469///GO:0030890///GO:0042102///GO:0050857///GO:0045860///GO:0006470///GO:0006470///GO:0051726///GO:0033261///GO:0051209	GO:0005925///GO:0005887///GO:0016020
# HEEBO-004-CTRL4P2	protein tyrosine phosphatase, receptor type, C	PTPRC	5788				Homo sapiens protein tyrosine phosphatase, receptor type, C (PTPRC), transcript variant 1, mRNA	115385975	NM_002838				1q31-q32	Chromosome 1, NC_000001.9 (196874848..196992953)	hydrolase activity///phosphoric monoester hydrolase activity///protein binding///protein kinase binding///protein kinase binding///protein tyrosine phosphatase activity///protein tyrosine phosphatase activity///receptor activity///transmembrane receptor protein tyrosine phosphatase activity	B cell proliferation///B cell receptor signaling pathway///T cell differentiation///T cell receptor signaling pathway///T cell receptor signaling pathway///cell surface receptor linked signal transduction///defense response to virus///immunoglobulin biosynthetic process///negative regulation of T cell mediated cytotoxicity///negative regulation of cytokine and chemokine mediated signaling pathway///negative regulation of protein kinase activity///negative regulation of protein kinase activity///positive regulation of B cell proliferation///positive regulation of T cell proliferation///positive regulation of antigen receptor-mediated signaling pathway///positive regulation of protein kinase activity///protein amino acid dephosphorylation///protein amino acid dephosphorylation///regulation of cell cycle///regulation of progression through S phase///release of sequestered calcium ion into cytosol	focal adhesion///integral to plasma membrane///membrane	GO:0016787///GO:0016791///GO:0005515///GO:0019901///GO:0019901///GO:0004725///GO:0004725///GO:0004872///GO:0005001	GO:0042100///GO:0050853///GO:0030217///GO:0050852///GO:0050852///GO:0007166///GO:0051607///GO:0002378///GO:0001915///GO:0001960///GO:0006469///GO:0006469///GO:0030890///GO:0042102///GO:0050857///GO:0045860///GO:0006470///GO:0006470///GO:0051726///GO:0033261///GO:0051209	GO:0005925///GO:0005887///GO:0016020
# ..... MORE ....
# HEEBO-128-CTRL128P22																				
# HEEBO-128-CTRL128P23																				
# HEEBO-128-CTRL128P24																				
# !platform_table_end

################## file_biomart ##################

### FORMAT ###
# csv file. 
# *NOTE*: some "HGNC symbol" fields may be EMPTY - not all ENSEMBL IDs have mapping to HGNC

### EXAMPLE ###
# Ensembl Gene ID,Chromosome Name,Gene Start (bp),Gene End (bp),Strand,Gene type,HGNC symbol,Source (gene),Status (gene)
# ENSG00000223116,13,23551994,23552136,-1,miRNA,,ensembl,NOVEL
# ENSG00000233440,13,23708313,23708703,1,pseudogene,HMGA1P6,havana,KNOWN
# ENSG00000207157,13,23726725,23726825,-1,misc_RNA,RNY3P4,ensembl,KNOWN
# ENSG00000229483,13,23743974,23744736,-1,lincRNA,LINC00362,havana,NOVEL
# ENSG00000252952,13,23791571,23791673,-1,snRNA,RNU6-58P,ensembl,KNOWN
# ENSG00000235205,13,23817659,23821323,1,pseudogene,TATDN2P3,havana,KNOWN
# ENSG00000232849,13,93708910,93710179,1,lincRNA,LINC00363,havana,NOVEL
# ENSG00000236803,13,23841645,23843289,1,pseudogene,SDAD1P4,havana,KNOWN

###################################### DESIGN ######################################
#1) Read

###################################### START SCRIPT ######################################

file_annot = "/Users/pascaltimshel/Dropbox/0_projects/p_EAv2/data_BrainCloud/mapping/GPL4611.annot"

file_biomart = "/Users/pascaltimshel/Dropbox/0_Projects/p_snpsnap/data/misc/biomart_download-2015-02-26-snpsnap_production_v2-ensembl-release_GRCh37.p13_processed-GENCODE.chrosome_clean.unique-ENSG_ID.csv"


###################################### Make mapping dict ######################################
dict_



with open("depict_FDR-significant_gene_sets_gwas_catalog.txt", 'w') as fh_out:
	fh_out.write( "gwas_trait\tgene_set\tp_value_nominal\n" )
	for file_in in gene_set_raw_files:
		gwas_trait = os.path.basename(file_in).split("_affymetrix")[0]
		print "gwas_trait = {}".format(gwas_trait)
		with open(file_in, 'r') as fh_in:
			header = next(fh_in) # skipping header (or we make a small validation)
			if not header.split('\t')[0] == "Original gene set ID":
				print "ERROR: got wrong header. Something might be wrong..."
			for line in fh_in:
				fields = line.strip("\n").split("\t")
				if not len(fields) == 4:
					print "ERROR: got wrong number of fields"
				gene_set_id = fields[0]
				p_value_nominal = fields[2]
				FDR_flag = fields[3]
				if FDR_flag == "Yes":
					fh_out.write( "{}\t{}\t{}\n".format(gwas_trait, gene_set_id, p_value_nominal) )

print "DONE"

