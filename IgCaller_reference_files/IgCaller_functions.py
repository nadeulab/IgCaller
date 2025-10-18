# IgCaller_functions [v2.0]

# modules
import subprocess
import sys
import os
import platform
import regex as re
import numpy as np
import itertools
import operator
from collections import Counter
from scipy import stats
from statistics import mean
from statistics import median
import gzip
import pickle
from difflib import SequenceMatcher
from Bio.Align import PairwiseAligner
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns

# dicts
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N', 'X': 'N', '[': ']', ']': '[', '(': ')', ')': '('}

changeAmbiguousBases = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'N', 'Y': 'N', 'S': 'N', 'W': 'N', 'K': 'N', 'M': 'N', 'B': 'N', 'V': 'N', 'D': 'N', 'H': 'N', 'N': 'N', 'X': 'N', '[': '[', ']': ']', '(': '(', ')': ')'}

tripletsToAA = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
				'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
				'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
				'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
				'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
				'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
				'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
				'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
				'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
				'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
				'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
				'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
				'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
				'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
				'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
				'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
				}

# functions
def smithwaterman(x, y, match_score=5, mismatch_cost=4, gap_cost=8):
    # scoring matrix
	M = np.zeros((len(x) + 1, len(y) + 1), int) # +1 because of the zero column and zero row
	for i, j in itertools.product(range(1, M.shape[0]), range(1, M.shape[1])): 
		match = M[i - 1, j - 1] + (match_score if x[i - 1] == y[j - 1] else - mismatch_cost)
		delete = M[i - 1, j] - gap_cost
		insert = M[i, j - 1] - gap_cost
		M[i, j] = max(match, delete, insert, 0)
	# Get maximum score
	M_flip = np.flip(np.flip(M, 0), 1)
	i_, j_ = np.unravel_index(M_flip.argmax(), M_flip.shape)
	i, j = np.subtract(M.shape, (i_ + 1, j_ + 1))
	return M[i,j]

def flagToCustomBinary(flag):
	binary = format(int(flag), "b")
	binary = ("0"*(12-len(binary)))+binary
	return(binary[::-1])

def CustomSequenceMatcher(seqA, seqB, seqType):
	matches = 0
	total = 0
	if len(seqA) == len(seqB) and seqA != "NA" and seqB != "NA":
		for pos in range(len(seqA)):
			if seqType == "nt" and seqA[pos] == "N" or seqB[pos] == "N": 
				continue
			if seqA[pos] == seqB[pos]: 
				matches += 1
			total += 1
	return(0 if total == 0 else matches/total)

def getGeneralInfo(GENE, chrom, genomeVersion, inputsFolder, chrAnnot, annotateSeqDB):
	Dseqs = "NA"
	germline_db_D = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_Dfake.fasta.txt"
	if GENE == "IGH":
		chromGene = chrom+"14"	
		germline_db_J = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_IGHJ.fasta.txt"
		germline_db_D = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_IGHD.fasta.txt"
		germline_db_V = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_IGHV.fasta.txt"
		if genomeVersion == "hg19":
			coordsToSubsetLocus = chromGene+":106329000-107283600"
			bedFile = inputsFolder+"/hg19/"+chrAnnot+"/wgEncodeGencodeBasicV19_hg19_IGH_genes_VJ.bed"
			Dseqs = inputsFolder+"/hg19/"+chrAnnot+"/DB_D_genes_seq_wgEncodeGencodeBasicV19_hg19.txt"
		else:
			coordsToSubsetLocus = chromGene+":105853198-106889844"
			bedFile = inputsFolder+"/hg38/"+chrAnnot+"/GencodeV40_hg38_IGH_genesVJ.bed"
			Dseqs = inputsFolder+"/hg38/"+chrAnnot+"/DB_D_genes_seq_GencodeV40_hg38.txt"
	
	elif GENE == "IGL":
		chromGene = chrom+"22"
		germline_db_J = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_IGLJ.fasta.txt"
		germline_db_V = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_IGLV.fasta.txt"
		if genomeVersion == "hg19": 
			coordsToSubsetLocus = chromGene+":22384800-23263900"
			bedFile = inputsFolder+"/hg19/"+chrAnnot+"/wgEncodeGencodeBasicV19_hg19_IGL_genes_VJ.bed"
		else: 
			coordsToSubsetLocus = chromGene+":22016076-22932913"
			bedFile = inputsFolder+"/hg38/"+chrAnnot+"/GencodeV40_hg38_IGL_genesVJ.bed"

	elif GENE == "IGK":
		chromGene = chrom+"2"
		germline_db_J = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_IGKJ.fasta.txt"
		germline_db_V = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_IGKV.fasta.txt"
		if genomeVersion == "hg19":
			coordsToSubsetLocus = chromGene+":89131589-90274600"
			bedFile = inputsFolder+"/hg19/"+chrAnnot+"/wgEncodeGencodeBasicV19_hg19_IGK_genes_VJ.bed"
		else:
			coordsToSubsetLocus = chromGene+":88822278-90245370"
			bedFile = inputsFolder+"/hg38/"+chrAnnot+"/GencodeV40_hg38_IGK_genesVJ.bed"
	
	elif GENE == "CSR":
		chromGene = chrom+"14"
		germline_db_J = "NA"
		germline_db_D = "NA"
		germline_db_V = "NA"
		if genomeVersion == "hg19":
			coordsToSubsetLocus = chromGene+":106055000-106329000"
			bedFile = inputsFolder+"/hg19/"+chrAnnot+"/hg19_Huebschmann_et_al_switch_regions.bed"
		else:
			coordsToSubsetLocus = chromGene+":105576937-105862797"
			bedFile = inputsFolder+"/hg38/"+chrAnnot+"/hg38_liftOver_Huebschmann_et_al_switch_regions.bed"
	
	elif GENE == "TRA":
		chromGene = chrom+"14"
		germline_db_J = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_TRAJ.fasta.txt"
		germline_db_V = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_TRAV.fasta.txt"
		if genomeVersion == "hg19":
			coordsToSubsetLocus = chromGene+":22089991-23014042"
			bedFile = inputsFolder+"/hg19/"+chrAnnot+"/gencode.v19.annotation_hg19_TRA_genesVJ.bed"
		else:
			coordsToSubsetLocus = chromGene+":21611904-22562132"
			bedFile = inputsFolder+"/hg38/"+chrAnnot+"/GencodeV40_hg38_TRA_genesVJ.bed"

	elif GENE == "TRB":
		chromGene = chrom+"7"
		germline_db_J = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_TRBJ.fasta.txt"
		germline_db_D = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_TRBD.fasta.txt"
		germline_db_V = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_TRBV.fasta.txt"
		if genomeVersion == "hg19":
			coordsToSubsetLocus = chromGene+":141999017-142511084"
			bedFile = inputsFolder+"/hg19/"+chrAnnot+"/gencode.v19.annotation_hg19_TRB_genesVJ.bed"
			Dseqs = inputsFolder+"/hg19/"+chrAnnot+"/DB_TRBD_genes_seq__gencode.v19_hg19.txt"
		else:
			coordsToSubsetLocus = chromGene+":142289011-142823287"
			bedFile = inputsFolder+"/hg38/"+chrAnnot+"/GencodeV40_hg38_TRB_genesVJ.bed"
			Dseqs = inputsFolder+"/hg38/"+chrAnnot+"/DB_TRBD_genes_seq_GencodeV40_hg38.txt"
	
	elif GENE == "TRG":
		chromGene = chrom+"7"
		germline_db_J = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_TRGJ.fasta.txt"
		germline_db_V = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_TRGV.fasta.txt"
		if genomeVersion == "hg19":
			coordsToSubsetLocus = chromGene+":38292981-38407770"
			bedFile = inputsFolder+"/hg19/"+chrAnnot+"/gencode.v19.annotation_hg19_TRG_genesVJ.bed"
		else:
			coordsToSubsetLocus = chromGene+":38230024-38378055"
			bedFile = inputsFolder+"/hg38/"+chrAnnot+"/GencodeV40_hg38_TRG_genesVJ.bed"

	elif GENE == "TRD":
		chromGene = chrom+"14"
		germline_db_J = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_TRDJ.fasta.txt"
		germline_db_D = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_TRDD.fasta.txt"
		germline_db_V = inputsFolder+"/igblast/db/"+annotateSeqDB.upper()+"_HomoSapiens_TRDV.fasta.txt"
		if genomeVersion == "hg19":
			coordsToSubsetLocus = chromGene+":22382239-22948690"
			bedFile = inputsFolder+"/hg19/"+chrAnnot+"/gencode.v19.annotation_hg19_TRD_genesVJ.bed"
			Dseqs = inputsFolder+"/hg19/"+chrAnnot+"/DB_TRDD_genes_seq_gencode.v19_hg19.txt"
		else:
			coordsToSubsetLocus = chromGene+":21914138-22479614"
			bedFile = inputsFolder+"/hg38/"+chrAnnot+"/GencodeV40_hg38_TRD_genesVJ.bed"
			Dseqs = inputsFolder+"/hg38/"+chrAnnot+"/DB_TRDD_genes_seq_GencodeV40_hg38.txt"

	if genomeVersion == "hg19":
		snps_file = inputsFolder+"/hg19/PoN/SNPsInJandVgenes_hg19.txt"
	else:
		snps_file = inputsFolder+"/hg38/PoN/SNPsInJandVgenes_hg38.txt"
	
	return(chromGene, coordsToSubsetLocus, bedFile, Dseqs, germline_db_J, germline_db_D, germline_db_V, snps_file)
	
def convertSamToAnnotatedTable(miniSamT, chromGene, GENE, minimumNumberOfNucleotidesSoft):
	
	samfile = open(miniSamT, "r")
	
	store = {}
	
	for i in samfile:
		
		w = i.rstrip("\n").split("\t")

		if w[5] != "*" and w[6] == "=":
			sa = []
			# analyse values after qualities, starting with SA...
			for x in (w[11:]):
				if x.startswith("SA:Z"):
					sa.append(x)
			# takes values from position 0 to 9th and appends SA... in the 10th
			if len(sa) != 0:
				sa = sa[0]
				saSplit = sa.split(",")
				# chromosome of SA must be equal to chromGene:
				if saSplit[0].split(":")[2] == chromGene:
					# force SA alignment close to insertSize if split-insertSize
					if abs(int(w[8])) > 10000:
						if abs(int(saSplit[1]) - int(w[7])) <= 1000: # 1000 bp considering standard insert size range
							w = w[:10]+[sa]
						else:
							w = w[:10]+["NA"]
					else:
						w = w[:10]+[sa]
				else:
					w = w[:10]+["NA"]
			else:
				w = w[:10]+["NA"]
				
			w.append("NA") # temporal NA to add later "split" or "insertSize" on 11th position

			split1 = re.findall(r'[A-Za-z]|[0-9]+', w[5])
			two1 = [split1[x:x+2] for x in range(0, len(split1),2)]
			split2 = w[10].split(',') # split SA...
			
			if GENE != "CSR" and sum([int(i[0]) for i in two1 if "S" in i]) >= minimumNumberOfNucleotidesSoft: # split if S >= minimumNumberOfNucleotidesSoft in cigar ## removed from code '(w[10].startswith("SA:Z") or )' and from comment 'SA or ' due to incorrect split annotation when "alt" contigs in ref genome [github issue #5]
				# check if split, split-insertSize (based on type of SV) or just an insertSize:
				if abs(int(w[8])) <= 10000: 
					w[11] = "split"
				else:
					# get type of alteration
					if int(w[1]) in [97, 161, 99, 163, 2145, 2209, 2147, 2211] and int(w[8]) > 0: svType = "Deletion" # reads: -> <-
					elif int(w[1]) in [145, 81, 147, 83, 2193, 2193, 2195, 2131] and int(w[8]) < 0: svType = "Deletion" # reads: -> <-
					elif int(w[1]) in [65, 129, 2113, 2177]: svType = "Inversion2" # reads -> ->
					elif int(w[1]) in [113, 177, 2161, 2225]: svType = "Inversion1" # reads <- <-
					else: svType = "other"
					# get if split-insertSize or just insertSize based on cigar and type of alteration
					if min([count for count, item in enumerate(two1) if "M" in item]) < min([count for count, item in enumerate(two1) if "S" in item]): cig = "MS"
					else: cig = "SM"
					if int(w[8]) > 0 and cig == "MS" and svType == "Deletion":
						w[11] = "split-insertSize"
					elif int(w[8]) < 0 and cig == "SM" and svType == "Deletion":
						w[11] = "split-insertSize"
					elif int(w[8]) < 0 and cig == "MS" and svType == "Inversion2":
						w[11] = "split-insertSize"
					elif int(w[8]) > 0 and cig == "SM" and svType == "Inversion1":
						w[11] = "split-insertSize"
					else: w[11] = "insertSize"
				
			elif abs(int(w[8])) > 10000:  # insertSize if insert size > 10000 [ J-V = 70000bp aprox ]
				w[11] = "insertSize"
			
			if w[11] == "split" or w[11] == "split-insertSize": # add 2 columns
				
				# 1st_pos
				two1 = [split1[x:x+2] for x in range(0, len(split1),2)]
				if len([count for count, item in enumerate(two1) if "S" in item]) == 0:
					firstpos = sum([int(i[0]) for i in two1 if "M" in i or "D" in i]) - 1
				elif min([count for count, item in enumerate(two1) if "M" in item]) < min([count for count, item in enumerate(two1) if "S" in item]):
					firstpos = sum([int(i[0]) for i in two1 if "M" in i or "D" in i]) - 1
				else:
					firstpos = 0
				
				# first PosSplit
				a = int(w[3]) + firstpos
				w.append(a)
				
				# second PosSplit
				if w[10].startswith("SA:Z"):
					
					split3 = re.findall(r'[A-Za-z]|[0-9]+', split2[3])
					two2 = [split3[x:x+2] for x in range(0, len(split3),2)]
					
					# get if overlapping same base/s at each break point... we substract this overlapping bases on the first break below
					diffMappingBases = sum([int(i[0]) for i in two2 if "M" in i or "D" in i]) - sum([int(i[0]) for i in two1 if "S" in i]) 
					
					if len([count for count, item in enumerate(two2) if "S" in item]) == 0: 
						secpos = sum([int(i[0]) for i in two2 if "M" in i or "D" in i]) - 1 
					elif min([count for count, item in enumerate(two2) if "M" in item]) < min([count for count, item in enumerate(two2) if "S" in item]):
						secpos = sum([int(i[0]) for i in two2 if "M" in i or "D" in i]) - 1 
					else:
						secpos = 0
						
					b = int(split2[1]) + secpos  # second PosSplit
					w.append(b)
					
					# sort (smaller first)
					w[-2:] = min(w[-2:]), max(w[-2:])
					
					if diffMappingBases > 0: 
						
						if int(flagToCustomBinary(w[1])[4]) == 0: strand = "+"
						else: strand = "-"
						strandSA = w[10].split(",")[2]
						
						# reads <- <-  (Inversion1) => sum diffMappingBases at first read soft clipped position
						if ( int(flagToCustomBinary(w[1])[4]) == 0 and strand == "-" and strandSA == "+" ) or ( int(flagToCustomBinary(w[1])[4]) == 1 and strand == "+" and strandSA == "-" ):
							w[-2] = w[-2] + diffMappingBases
						
						else:  # reads -> <- or -> -> => substract diffMappingBases at first read soft clipped position
							w[-2] = w[-2] - diffMappingBases
				
				# for non split-reads add NA
				else:
					w.append("NA")
				
				if w[11] == "split": w.extend(["NA"]* 2) # no insertSize info
				
			if w[11] == "insertSize" or w[11] == "split-insertSize": 
				
				if w[11] == "insertSize": w.extend(["NA"]* 2) # no split info
				
				if int(flagToCustomBinary(w[1])[4]) == 0: # positive strand 
					two = [split1[x:x+2] for x in range(0, len(split1),2)]
					converted = sum([int(i[0]) for i in two if  "I" not in i and "S" not in i]) - 1 
					d = int(w[3]) + converted
				else: # negative strand  
					d = int(w[3]) 
				
				if 	int(flagToCustomBinary(w[1])[6]) == 1 and int(flagToCustomBinary(w[1])[4]) == 0 and int(w[8]) > 0: # first in pair and positive strand and insertSize > 0 (97 -> del) (65 -> inv)
					w.append(d)
					w.append("NA")
				elif int(flagToCustomBinary(w[1])[6]) == 0 and int(flagToCustomBinary(w[1])[4]) == 1 and int(w[8]) < 0: # second in pair and negative strand and insertSize < 0 (145 <- del) (177 <- inv)
					w.append("NA")
					w.append(d)
				elif int(flagToCustomBinary(w[1])[6]) == 1 and int(flagToCustomBinary(w[1])[4]) == 1 and int(w[8]) < 0: # first in pair and negative strand and insertSize < 0 (81 <- del)
					w.append("NA")
					w.append(d)
				elif int(flagToCustomBinary(w[1])[6]) == 0 and int(flagToCustomBinary(w[1])[4]) == 0 and int(w[8]) > 0: # second in pair and positive strand and insertSize > 0 (161 -> del)
					w.append(d)
					w.append("NA")
				elif int(flagToCustomBinary(w[1])[6]) == 0 and int(flagToCustomBinary(w[1])[4]) == 0 and int(w[8]) < 0: # second in pair and positive strand and insertSize < 0 (129 -> inv)
					w.append("NA")
					w.append(d)
				elif int(flagToCustomBinary(w[1])[6]) == 1 and int(flagToCustomBinary(w[1])[4]) == 1 and int(w[8]) > 0: # first in pair and negative strand and insertSize > 0 (113 <- inv)
					w.append(d)
					w.append("NA")
				else:
					w.append("NA")
					w.append("NA")
					if w[11] == "insertSize": w[11] = "NA" # if no deletion or inversions considered above, remove info if no split				
			
			if w[11] == "NA":
				w.extend(["NA"]* 4) # if it is neither insertSize nor split we add the five empty columns for it
			
			# dict store:
			if w[0] not in store:
				if w[11] == "split" or w[11] == "insertSize" or w[11] == "split-insertSize":
					store[w[0]] = w
			
			else:
				if w[11] == "split" or w[11] == "split-insertSize":
					if store[w[0]][11] == "insertSize":
						if store[w[0]][14] != "NA":
							w[14] = store[w[0]][14] 
						else:
							w[15] = store[w[0]][15]
						
						store[w[0]] = w
				
				if (w[11] == "insertSize" and store[w[0]][11] == "insertSize") or (w[11] == "insertSize" and store[w[0]][11] == "split-insertSize"):
					if w[14] != "NA":
						store[w[0]][14] = w[14]
					else:
						store[w[0]][15] = w[15]
				if (w[11] == "split-insertSize" and store[w[0]][11] == "split-insertSize"):
					if w[14] != "NA":
						store[w[0]][14] = w[14]
						store[w[0]][12] = w[14]
					else:
						store[w[0]][15] = w[15]
						store[w[0]][13] = w[15]
	
	return(store)
	
def findJandVgenes(annot_table, bedFile, GENE, genomeVersion):
	insertsplit = open(annot_table, "r")
	
	JV_list = []

	for i in insertsplit:
		w = i.rstrip("\n").split("\t")

		# recover information from insertSize if second position has not been found
		if w[13] == "NA" and w[15] == "NA" and abs(int(w[8])) > 10000:
			w[15] = w[7]
			if w[14] == "NA": w[14] = w[12]
			if int(w[14]) > int(w[15]):
				w14 = w[14]
				w[14] = w[15]
				w[15] = w14
		
		w.extend(["NA"]* 4) # we add 4 new columns where we will add new information
		for j in range(12,16):
			if w[j] != "NA": 
				VDJ = open(bedFile, "r") # bed file for V/J genes
				for k in VDJ:
					v = k.rstrip("\n").split("\t")
					if GENE == "IGH": # orientation J - V (-)
						if v[3].startswith("IGHJ"):
							window1 = 0 if w[11] != "insertSize" else 150 if v[3] == "IGHJ1" else 300
							window2 = 10
						else:
							window1 = 10
							window2 = 0
					elif GENE == "IGK": # orientation J - V (-)
						if v[3].startswith("IGKJ"):
							window1 = 250 if w[11] == "insertSize" else 0
							window2 = 10
						else:
							window1 = 10
							window2 = 5 # 5bp added for inversions just after V
					elif GENE == "IGL": # orientation V - J (+)
						if v[3].startswith("IGLJ"):
							window1 = 10
							window2 = 300 if w[11] == "insertSize" else 0
						else:
							window1 = 0
							window2 = 10
					elif GENE == "CSR": # class switch
						window1 = 0
						window2 = 0
					
					elif GENE == "TRA": # orientation V - J (+)
						if v[3].startswith("TRAJ"):
							window1 = 10
							window2 = 150 if w[11] == "insertSize" else 0
						else:
							window1 = 0
							window2 = 10
					
					elif GENE == "TRB": # orientation V - J (+) [some V genes in - strand in hg19]
						if v[3].startswith("TRBJ"):
							window1 = 10
							window2 = 50 if w[11] == "insertSize" else 0
						else:
							window1 = 0 if genomeVersion == "hg38" else 5
							window2 = 10
					
					elif GENE == "TRG": # orientation J - V (-)
						if v[3].startswith("TRGJ"):
							window1 = 300 if w[11] == "insertSize" else 0
							window2 = 10
						else:
							window1 = 10
							window2 = 0

					elif GENE == "TRD": # orientation V - J (+)
						if v[3].startswith("TRDJ"):
							window1 = 10
							window2 = 300 if w[11] == "insertSize" else 0
						else:
							window1 = 0
							window2 = 10

					if int(w[j]) >= int(v[1])-window1 and int(w[j]) <= int(v[2])+window2: # we check in which V or J is the position included, considering the defined windows
						w[j+4] = v[3] # we append the corresponding V or J, 4 positions to the right to our table
						break
				VDJ.close()

		# We classify reads according to their orientation (flags):
		if w[11] == "insertSize" or ( w[11] == "split-insertSize" and w[10] == "NA" ):
			
			# reads: -> <-
			if int(w[1]) in [97, 161, 99, 163, 2145, 2209, 2147, 2211] and int(w[8]) > 0: w.append("Deletion")
			elif int(w[1]) in [145, 81, 147, 83, 2193, 2193, 2195, 2131] and int(w[8]) < 0: w.append("Deletion")
			
			# reads -> ->
			elif int(w[1]) in [65, 129, 2113, 2177]: w.append("Inversion2")
		
			# reads <- <-
			elif int(w[1]) in [113, 177, 2161, 2225]: w.append("Inversion1")
			
			# other potential SV not considered:
			else: w.append("NA")
		
		else:
			# no information of soft clipped map:
			if w[10] == "NA": 
				w.append("NotComplete") 
			else: # split			
				# get strands
				if int(flagToCustomBinary(w[1])[4]) == 0: strand = "+"
				else: strand = "-"
				strandSA = w[10].split(",")[2]
				
				# reads: -> <-
				if strand == "+" and strandSA == "+" and int(w[8]) > 0: w.append("Deletion")
				elif strand == "-" and strandSA == "-" and int(w[8]) < 0: w.append("Deletion")
				
				elif (strand == "+" and strandSA == "-") or (strand == "-" and strandSA == "+"):
					# reads -> ->
					if int(w[1]) in [65, 129, 99, 97, 147, 145, 2113, 2177, 2147, 2145]: w.append("Inversion2")

					# reads <- <-
					elif int(w[1]) in [113, 177, 2161, 2225]: w.append("Inversion1")

					else: w.append("NA")
				
				# other potential SV not considered:
				else: w.append("NA")
		
		if len(set(w[-5:-1])) > 1 and w[-1] != "NA": # not NA only and DELETION/INVERSION specification
			JV_list.append("%s\n" %"\t".join([str(x) for x in w]))

	insertsplit.close() # we have a table with ID (1 column), insertsize or split positions (4 columns), corresponding VDJ genes (4 columns)
	
	return(JV_list)

def findCombinationsJandV(annot_table_JV, GENE):
	ANNOT_TABLE_JV = open(annot_table_JV, "r") # we use previous output file as input file
	l = list()

	for i in ANNOT_TABLE_JV:
		w = i.rstrip("\n").split("\t")
		r = []
		s = []
		for j in w[-5:-1]: # we take the four columns corresponding to V and J types for each insertSize and split positions
			if j != "NA": # if there exists information
				if j[3] not in r:
					r.append(j[3]) # we save letters corresponding to V and J ex: IGHV1-3 -> V, IGHJ1 -> J
					s.append([j, w[-5:-1].index(j)]) # we create sublists with V and J variants and the column number they belong
				else: # if list already contains the letter we are analysing (ex: s = [[V, 1],[J, 2]] and we are analysing V), we eliminate all information, as V can only match J and viceversa
					r = []
					s = []
			
			if w[-5:-1].index(j) == 1 or w[-5:-1].index(j) == 3: # second split or second insertSize
				if len(s) == 2: # if we find any type VJ, JV
					if GENE != "CSR":
						if [s[0][0]+" - "+s[1][0], w[-1]] not in l: # if not in list of pairs, add..
							l.append([s[0][0]+" - "+s[1][0], w[-1]])
					else:
						if "M" in r and [s[0][0]+" - "+s[1][0], w[-1]] not in l: # if it is class switch it must contain M
							l.append([s[0][0]+" - "+s[1][0], w[-1]])
				s = []
				r = []
	
	ANNOT_TABLE_JV.close()
	
	# clean l
	if GENE == "IGH": geneCombinations = ["JD", "DV", "JV"]
	elif GENE == "IGK": geneCombinations = ["JV"]
	elif GENE == "IGL": geneCombinations = ["VJ"]
	elif GENE == "CSR": geneCombinations = ["AM", "EM", "GM"]
	elif GENE == "TRA": geneCombinations = ["VJ"]
	elif GENE == "TRB": geneCombinations = ["VJ", "VD", "DJ"]
	elif GENE == "TRG": geneCombinations = ["JV"]
	elif GENE == "TRD": geneCombinations = ["VJ", "VD", "DJ"]
	ll = []
	for i in l:
		XX = i[0].split(" - ")[0][3]+i[0].split(" - ")[1][3]
		if "Kde" in i[0] or "RSS" in i[0] or XX in geneCombinations:
			ll.append(i)
	return(ll)

def assignPositionsToJandV(l, annot_table_JV, seq, GENE):
	
	VJ_positions = {} # we store pairs J-V positions and if they come from split/insertsize or both in some cases
	data = {} # we store count of pairs and individuals J/V by positions (from split) and by gene names (by insertSize) 
	pos = {}
	
	for k12 in l: # iterates over every sublist of pairs in list
		
		ANNOT_TABLE_JV = open(annot_table_JV, "r")
		for j in ANNOT_TABLE_JV:
			w = j.rstrip("\n").split("\t")
			
			for idx in [12, 14]:
				a = idx # first pos split/insert
				b = idx+1 # second pos split/insert
				
				if w[a] == "NA" or w[b] == "NA": continue
				
				srt = w[a+4]+" - "+w[b+4]
				srtPos = w[a]+" - "+w[b]
				readtype = w[20]
				
				if srt == k12[0] and readtype == k12[1]: # k12[0] = each pair, k12[1] = Deletion, Inversion1, Inversion2
					if srt not in VJ_positions:
						VJ_positions[k12[0]] = [[srtPos, "split" if idx == 12 else "insertSize", readtype]] 
					else:
						VJ_positions[k12[0]].append([srtPos, "split" if idx == 12 else "insertSize", readtype])
					
					if idx == 12:
						if srtPos+" - "+k12[1] not in data:
							data[srtPos+" - "+readtype] = 1
						else:
							data[srtPos+" - "+readtype] += 1
							
					else:
						if srt+" - "+readtype not in data:
							data[srt+" - "+readtype] = 1
						else:
							data[srt+" - "+readtype] += 1
		
		ANNOT_TABLE_JV.close()
		
		for key in VJ_positions: # dictionary of IGH from V and J position of start and end (specific for each rearrangement)
			spl = [] # save split position
			ins = [] # save insertSize position
			
			for m in VJ_positions[key]:
				if m[1] == "split":
					spl.append([m[0], m[2]]) # if info comes from split we save it on one list
				else:
					ins.append([m[0], m[2]]) # if info comes from insertSize we save it on another list
			
			# info from paired-split
			if len(spl) != 0:
				UNIQUEspl = []
				for s in spl:
					if s not in UNIQUEspl:
						UNIQUEspl.append(s) 
				pos[key] = UNIQUEspl
			
			# info from single-split (make all possible combinations)
			if len(spl) == 0 or all(data[el[0]+" - "+el[1]] <= 2 for el in UNIQUEspl):
				ANNOT_TABLE_JV = open(annot_table_JV, "r")
				Jpos = []
				Vpos = []
				
				svClassInsert = [] # to get sv class from inserSize pairs
				for j in ANNOT_TABLE_JV:
					w = j.rstrip("\n").split("\t")
					if w[11].startswith("split") and key.split(" - ")[0] == w[16]:
						Jpos.append([w[12], w[20]]) # save all possible J positions comming individual splits
					elif w[11].startswith("split") and key.split(" - ")[1] == w[16]:
						Vpos.append([w[12], w[20]]) # save all possible V positions comming from individual splits
					
					if w[18] == key.split(" - ")[0] and w[19] == key.split(" - ")[1]: # key by insertSize reads to get SV class
						svClassInsert.append(w[20])
				
				ANNOT_TABLE_JV.close()
				
				if svClassInsert != []:
					svClassInsert = Counter(svClassInsert).most_common(1)[0][0] # Simplify to most common sv class

					UNIQUEjpos = []
					for x in Jpos:
						if x[1] == "NotComplete" or x[1] == svClassInsert:
							if x[0] not in UNIQUEjpos:
								UNIQUEjpos.append(x[0])

					UNIQUEvpos = []
					for y in Vpos:
						if y[1] == "NotComplete" or y[1] == svClassInsert:
							if y[0] not in UNIQUEvpos:
								UNIQUEvpos.append(y[0])

					JV = [[x+" - "+y, svClassInsert] for x in UNIQUEjpos for y in UNIQUEvpos] # we create all possible combinations if they have equal read orientation

					if key in pos: pos[key].extend(JV)
					else: pos[key] = JV
				
				# info still no info, get info from paired-insertSize, unpaired insertSize and unpaired split
				if pos[key] == [] or (seq == "amplicon" and GENE == "IGH"): # exception for amplicon and IGH to get additional pairs by insertSize (only IGH due to N-D-N plus SHM)
					ANNOT_TABLE_JV = open(annot_table_JV, "r")
					Jpos = []
					Vpos = []
				
					svClassInsert = [] # to get sv class from inserSize pairs
					for j in ANNOT_TABLE_JV:
						w = j.rstrip("\n").split("\t")
						if key.split(" - ")[0] == w[16]:
							Jpos.append([w[12], w[20]]) 
						elif key.split(" - ")[1] == w[16]:
							Vpos.append([w[12], w[20]]) 
						if key.split(" - ")[0] == w[18]:
							Jpos.append([w[14], w[20]]) 
						elif key.split(" - ")[1] == w[18]:
							Vpos.append([w[14], w[20]]) 
						if key.split(" - ")[1] == w[19]:
							Vpos.append([w[15], w[20]]) 
					
						if w[18] == key.split(" - ")[0] and w[19] == key.split(" - ")[1]: # key by insertSize reads to get SV class
							svClassInsert.append(w[20])
						
					ANNOT_TABLE_JV.close()
					
					if svClassInsert != []:
						svClassInsert = Counter(svClassInsert).most_common(1)[0][0] # Simplify to most common sv class

						UNIQUEjpos = []
						for x in Jpos:
							if x[1] == "NotComplete" or x[1] == svClassInsert:
								if x[0] not in UNIQUEjpos:
									UNIQUEjpos.append(x[0])

						UNIQUEvpos = []
						for y in Vpos:
							if y[1] == "NotComplete" or y[1] == svClassInsert:
								if y[0] not in UNIQUEvpos:
									UNIQUEvpos.append(y[0])
						
						JV = [[x+" - "+y, svClassInsert] for x in UNIQUEjpos for y in UNIQUEvpos] # we create all possible combinations if they have equal read orientation
						
						if pos[key] == []: pos[key] = JV
						else: pos[key].extend([jvPair for jvPair in JV if jvPair not in pos[key]])
	
	return(VJ_positions, data, pos)

def addPositionsAndOccurrences(GENE, pos, bedFile, shortV, data):
	information = [] # list of sublists with pairs and information about them
	
	for key in pos:
		
		for i in range(len(pos[key])):
			
			keyJ = key.split(" - ")[0]
			keyV = key.split(" - ")[1]
			readtype = pos[key][i][1] # read type: INVERSION1,2 or DELETION
			keyPos = pos[key][i][0] # J-V positions split
			keyPosJ = int(keyPos.split(" - ")[0]) # last position for corresponding J
			keyPosV = int(keyPos.split(" - ")[1]) # first position for corresponding V
			partialJDRearrangement = "yes" if ( (GENE == "IGH" and keyV[3] == "D") or (GENE in ["TRB", "TRD"] and keyJ[3] == "D") ) else "no"

			row = [key, readtype] # we start adding V-J combination
			
			a = int(data.get(keyPos+" - "+readtype, 0)) # we find occurrences for corresponding V-J from split
			row.append(a)
			
			b = int(data.get(key+" - "+readtype, 0)) # we find occurrences for corresponding V-J from insertSize
			row.append(b)
			
			# find position in bed file (J)   
			c = "NA" 
			VDJ = open(bedFile, "r")
			for k in VDJ:
				v = k.rstrip("\n").split("\t")
				if keyJ == v[3]:
					if pos[key][i][1] == "Deletion":
						c = int(v[1]) if not keyJ.startswith(("IGL", "TRA", "TRB", "TRD")) or partialJDRearrangement == "yes" or shortV == "no" or int(v[2])-320 > keyPosJ else int(v[2])-320 # if IGL/TRA/TRB/TRD and shortV == "yes" => short V (note that keyJ is indeed V for IGL/TRA/TRB/TRD locus here)
						break
					elif pos[key][i][1] == "Inversion1":
						c = int(v[2])
						break
					elif pos[key][i][1] == "Inversion2":
						c = int(v[1])
						break
			row.extend(sorted([c, keyPosJ]))
			VDJ.close()	
			if row[-1] - row[-2] < 5 or ( keyJ.startswith(("IGL", "TRA", "TRB", "TRD")) and partialJDRearrangement == "no" and row[-1] - row[-2] < 50 ): continue # (note that keyJ is indeed V for IGL/TRA/TRB/TRD locus here)
			
			# add 0 for split reads in J used afterwards
			row.append(0)
			
			# find position in bed file (V) 
			h = "NA" 
			VDJ = open(bedFile, "r")
			for k in VDJ:
				v = k.rstrip("\n").split("\t")
				if keyV == v[3]:
					if pos[key][i][1] == "Deletion":
						h = int(v[2]) if shortV == "no" or keyV.startswith(("IGL", "TRA", "TRB", "TRD")) or partialJDRearrangement == "yes" or keyJ == "IGKKde" or keyJ == "IGKRSS" or int(v[1])+320 < keyPosV else int(v[1])+320 # if not IGL/TRA/TRB/TRD, not IGKKde/IGKRSS, and shortV == "yes" => short V
						break
					elif pos[key][i][1] == "Inversion1":
						h = int(v[2]) if shortV == "no" or keyV.startswith(("IGL", "TRA", "TRB", "TRD")) or partialJDRearrangement == "yes" or keyJ == "IGKKde" or keyJ == "IGKRSS" or int(v[1])+320 < keyPosV else int(v[1])+320 # if not IGL/TRA/TRB/TRD, not IGKKde/IGKRSS, and shortV == "yes" => short V
						break
					elif pos[key][i][1] == "Inversion2":
						h = int(v[1]) if shortV == "no" or keyV.startswith(("IGL", "TRA", "TRB", "TRD")) or partialJDRearrangement == "yes" or keyJ == "IGKKde" or keyJ == "IGKRSS" or int(v[1])+320 > keyPosV else int(v[2])-320 # if not IGL/TRA/TRB/TRD, not IGKKde/IGKRSS, and shortV == "yes" => short V
						break
			VDJ.close()		
			row.extend(sorted([h, keyPosV]))
			if row[-1] - row[-2] < 5 or ( not keyV.startswith(("IGL", "TRA", "TRB", "TRD")) and partialJDRearrangement == "no" and  row[-1] - row[-2] < 50 ): continue
			
			# add 0 for split reads in V used afterwards
			row.append(0)
			
			information.append(row)
	
	# keep only potential pairs with V length > 200 bp
	VendIndex = 8 if not GENE in ["IGL", "TRA", "TRB", "TRD"] else 5
	information = [i for i in information if ( ((i[VendIndex]-i[VendIndex-1])+1) > 200 or (i[0].startswith("IGKKde") or i[0].startswith("IGKRSS")) or "IGHD" in i[0] or "TRBD" in i[0] or "TRDD" in i[0] )]
	
	# return information
	return(information)
	
def cleanPositionsAndOccurrences(GENE, bedFile, information, highSensitivity, seq):

	informationClean = []

	# remove pairs initially supported mostly by insertSize reads if highSensitivity == "no"
	if highSensitivity == "no":
		for i in information:
			if i[2] > 1 or (i[2] == 1 and i[3] > 2):
				informationClean.append(i)

	# check if potential break is found in the expected region if no split-read support
	else:
		for i in information:
			
			geneJ = i[0].split(" - ")[1] if GENE in ["IGL", "TRA", "TRB", "TRD"] else i[0].split(" - ")[0]
			geneV = i[0].split(" - ")[0] if GENE in ["IGL", "TRA", "TRB", "TRD"] else i[0].split(" - ")[1]
			mechanism = i[1]
			partialRearrangement = "yes" if i[0].split(" - ")[0][3] == "D" or i[0].split(" - ")[1][3] == "D" else "no"

			# keep if split reads, partial rearrangement, or Kde/RSS
			if i[2] > 0 or partialRearrangement == "yes" or "IGKKde" in i[0] or "IGKRSS" in i[0]:
				informationClean.append(i)

			# if no split-read support, check if the potential breakpoints are close to the expected regions of the gene
			else:
				if seq != "amplicon" or GENE != "IGH": # no doing this if amplicon and IGH because initially found breakpoints may be far away from the start of gene depending on primer design (only IGH due to N-D-N plus SHM)
					
					# check position of break J
					breakJ = "NA"
					VDJ = open(bedFile, "r")
					for k in VDJ:
						v = k.rstrip("\n").split("\t")
						if geneJ == v[3]:
							breakJ = int(v[1]) if GENE in ["IGL", "TRA", "TRB", "TRD"] else int(v[2])
							leftWinJ = breakJ-10 if GENE in ["IGL", "TRA", "TRB", "TRD"] else breakJ-100
							rightWinJ = breakJ+100 if GENE in ["IGL", "TRA", "TRB", "TRD"] else breakJ+10
							potentialBreakJ = i[7] if GENE in ["IGL", "TRA", "TRB", "TRD"] else i[5]
							break
					VDJ.close()	

					if breakJ == "NA": continue
					if potentialBreakJ < leftWinJ or potentialBreakJ > rightWinJ: continue
					
					# check position of break V
					breakV = "NA"
					VDJ = open(bedFile, "r")
					for k in VDJ:
						v = k.rstrip("\n").split("\t")
						if geneV == v[3]:
							if mechanism == "Deletion":
								breakV = int(v[2]) if GENE in ["IGL", "TRA", "TRB", "TRD"] else int(v[1])
								leftWinV = breakV-100 if GENE in ["IGL", "TRA", "TRB", "TRD"] else breakV-10
								rightWinV = breakV+10 if GENE in ["IGL", "TRA", "TRB", "TRD"] else breakV+100
								potentialBreakV = i[5] if GENE in ["IGL", "TRA", "TRB", "TRD"] else i[7]
							elif mechanism == "Inversion1" and GENE == "TRB":
								breakV = int(v[1])
								leftWinV = breakV-10
								rightWinV = breakV+100
								potentialBreakV = i[4]				
							elif mechanism == "Inversion2" and GENE == "IGK":
								breakV = int(v[2])
								leftWinV = breakV-100
								rightWinV = breakV+10
								potentialBreakV = i[8]
							break
					VDJ.close()
					
					if breakV == "NA": continue
					if potentialBreakV < leftWinV or potentialBreakV > rightWinV: continue

				# if break close to the position where it should be found, append to informationClean
				informationClean.append(i)

		# hard cutoff to avoid excessive running time if highSensitivity == "yes" and too many potential rearrangements
		informationClean.sort(key=lambda p: round(p[2]*2 + p[3], 1), reverse=True)
		if len(informationClean) > 200:
			informationClean = informationClean[:200]
		
	# return informationClean
	return(informationClean)

def addReadNames(GENE, information, annot_table_JV):

	for i in information:
		
		geneJ = i[0].split(" - ")[0]
		geneV = i[0].split(" - ")[1]
		mech = i[1]
		posJ = i[4] if GENE == "TRB" and mech == "Inversion1" else i[5]
		posV = i[8] if GENE == "IGK" and mech == "Inversion2" else i[7]
		countSplit = 0
		countInsert = 0
		readNames = []

		ANNOT_TABLE_JV = open(annot_table_JV, "r")
		for j in ANNOT_TABLE_JV:
			w = j.rstrip("\n").split("\t")
			if w[11] == "insertSize" and geneJ == w[18] and geneV == w[19] and mech == w[20]:
				countInsert += 1
				readNames.append(w[0])
			elif w[12] != "NA" and w[13] != "NA":
				if posJ == int(w[12]) and posV == int(w[13]) and mech == w[20]:
					countSplit += 1
					readNames.append(w[0])
		ANNOT_TABLE_JV.close()
		i[2] = countSplit # will be the same number, but just for uniformity
		i[3] = countInsert # will be a different number - this is recalculated based in pure insertSize reads
		i.append(",".join(readNames))
	
	# sort and return
	information.sort(key=lambda p: round(p[2]*2 + p[3], 1), reverse=True)
	return(information)

def getJandVsequences(round, seq, phaseReadsBasedOnMutations, information, annot_table_JV, GENE, refGenome, snps_file, baseq, chromGene, bamN, pairedMode, miniBamT, miniBamN, depth, altDepth, tumorPurity, vafCutoff, vafCutoffNormal, pathToSamtools, threadsForSamtools, errLogMpileup):
	
	if round == "first" and GENE in ["IGL", "TRA", "TRB", "TRD"]: # if IGL/TRA/TRB/TRD, switch V <-> J info
		for i in information:
			i[0] = i[0].split(" - ")[1]+" - "+i[0].split(" - ")[0]
			i4 = i[4]
			i5 = i[5]
			i[4] = i[7]
			i[5] = i[8]
			i[7] = i4
			i[8] = i5

	for i in information:
		temporary = []
		phaseInfo = "NA"

		if round == "second" and phaseReadsBasedOnMutations == "no":
			i[10:10] = [0]
			i.append("")
			continue
		
		if round == "second" and phaseReadsBasedOnMutations == "yes":
			readsAlreadyUsed = i[16].split(",")+i[17].split(",")+i[18].split(",")
			readsPhasedUsingMutations = [] # to keep reads correctly phased using mutations
			readPairsRemovedFromPhased = [] # to keep reads removed from phasing because they lack mutations or the start of the read is found before the breakpoint

		for z in [4,7]: # to iterate over positions for J and V
			if "Kde" in i[0] or "RSS" in i[0]: # no sequence to retrieve
				temporary.extend(["NA"]*2) 
				
			elif i[z] != "NA" and i[z+1] != "NA":
				
				if refGenome is None:
					fr = " -r "
				else:
					fr = " -f "+refGenome+" -r "
					
				# mpileup tumor all reads
				subprocess.call(pathToSamtools+"samtools mpileup -d 0 -a -A -B -Q "+baseq+fr+chromGene+":"+str(i[z])+"-"+str(i[z+1])+" "+miniBamT+ " > "+miniBamT.replace(".bam", "_output_mpileup.tsv")+" 2>> "+errLogMpileup, shell=True) # allow -A (anomalous read pairs) in tumor sample only
				
				if round == "first":
					# mpileup tumor only reads spanning V-J
					J = i[0].split(" - ")[1] if GENE in ["IGL", "TRA", "TRB", "TRD"] else i[0].split(" - ")[0]
					V = i[0].split(" - ")[0] if GENE in ["IGL", "TRA", "TRB", "TRD"] else i[0].split(" - ")[1]
					readsSpanningJV = []
					ANNOT_TABLE_JV = open(annot_table_JV, "r")
					for readLine in ANNOT_TABLE_JV:
						readList = readLine.rstrip("\n").split("\t")
						if (J == readList[16] and V == readList[17]) or (J == readList[18] and V == readList[19]):
							readsSpanningJV.append(readList[0])
					ANNOT_TABLE_JV.close()

					if len(readsSpanningJV) > 0:
						readNamesFileTxt = miniBamT.replace(".bam", "_readNameReadsSpanningJV.txt")
						readNamesFile = open(readNamesFileTxt, "w")
						readNamesFile.write("\n".join(readsSpanningJV))
						readNamesFile.close()
						subprocess.call(pathToSamtools+"samtools view -@ "+threadsForSamtools+" -h -b -N "+readNamesFileTxt+" -o "+miniBamT.replace("miniBam.bam", "miniBam_readsSpanningJV.bam")+" "+miniBamT, shell=True)
						subprocess.call(pathToSamtools+"samtools index "+miniBamT.replace("miniBam.bam", "miniBam_readsSpanningJV.bam"), shell=True) 
						subprocess.call(pathToSamtools+"samtools mpileup -d 0 -a -A -B -Q "+baseq+fr+chromGene+":"+str(i[z])+"-"+str(i[z+1])+" "+miniBamT.replace("miniBam.bam", "miniBam_readsSpanningJV.bam")+ " > "+miniBamT.replace(".bam", "_output_mpileup_readsSpanningJV.tsv")+" 2>> "+errLogMpileup, shell=True)
						if os.stat(miniBamT.replace(".bam", "_output_mpileup_readsSpanningJV.tsv")).st_size == 0:
							O = open(miniBamT.replace(".bam", "_output_mpileup_readsSpanningJV.tsv"), "w")
							for missingPos in range(i[z], i[z+1]+1):
								O.write("%s\t%s\tNA\t0\tNA\tNA\n" %(chromGene, str(missingPos)))
							O.close()
					else:
						O = open(miniBamT.replace(".bam", "_output_mpileup_readsSpanningJV.tsv"), "w")
						for missingPos in range(i[z], i[z+1]+1):
							O.write("%s\t%s\tNA\t0\tNA\tNA\n" %(chromGene, str(missingPos)))
						O.close()
					
					# mpileup tumor only reads spanning exact breakpoints
					J = str(i[7]) if GENE == "TRB" and i[1] == "Inversion1" else str(i[8]) if GENE in ["IGL", "TRA", "TRB", "TRD"] else str(i[5])
					V = str(i[4]) if GENE in ["IGL", "TRA", "TRB", "TRD"] else str(i[8]) if GENE == "IGK" and i[1] == "Inversion2" else str(i[7])
					readsSpanningBreak = []
					ANNOT_TABLE_JV = open(annot_table_JV, "r")
					for readLine in ANNOT_TABLE_JV:
						readList = readLine.rstrip("\n").split("\t")
						if (J == readList[12] and V == readList[13]) or (J == readList[14] and V == readList[15]):
							readsSpanningBreak.append(readList[0])
					ANNOT_TABLE_JV.close()

					if len(readsSpanningBreak) > 0:
						readNamesFileTxt = miniBamT.replace(".bam", "_readNameReadsSpanningBreak.txt")
						readNamesFile = open(readNamesFileTxt, "w")
						readNamesFile.write("\n".join(readsSpanningBreak))
						readNamesFile.close()
						subprocess.call(pathToSamtools+"samtools view -@ "+threadsForSamtools+" -h -b -N "+readNamesFileTxt+" -o "+miniBamT.replace("miniBam.bam", "miniBam_readsSpanningBreak.bam")+" "+miniBamT, shell=True)
						subprocess.call(pathToSamtools+"samtools index "+miniBamT.replace("miniBam.bam", "miniBam_readsSpanningBreak.bam"), shell=True) 
						subprocess.call(pathToSamtools+"samtools mpileup -d 0 -a -A -B -Q "+baseq+fr+chromGene+":"+str(i[z])+"-"+str(i[z+1])+" "+miniBamT.replace("miniBam.bam", "miniBam_readsSpanningBreak.bam")+ " > "+miniBamT.replace(".bam", "_output_mpileup_readsSpanningBreak.tsv")+" 2>> "+errLogMpileup, shell=True)
						if os.stat(miniBamT.replace(".bam", "_output_mpileup_readsSpanningBreak.tsv")).st_size == 0:
							O = open(miniBamT.replace(".bam", "_output_mpileup_readsSpanningBreak.tsv"), "w")
							for missingPos in range(i[z], i[z+1]+1):
								O.write("%s\t%s\tNA\t0\tNA\tNA\n" %(chromGene, str(missingPos)))
							O.close()
					else:
						O = open(miniBamT.replace(".bam", "_output_mpileup_readsSpanningBreak.tsv"), "w")
						for missingPos in range(i[z], i[z+1]+1):
							O.write("%s\t%s\tNA\t0\tNA\tNA\n" %(chromGene, str(missingPos)))
						O.close()
				
				else: # "second" or "third"
					# mpileup tumor only reads spanning the rearrangement (based on readnames assigned to it)
					if round == "second":
						readsSpanningRearrangement = i[-3].split(",")+i[-2].split(",")+i[-1].split(",")
					else: # "third"
						readsSpanningRearrangement = i[-4].split(",")+i[-3].split(",")+i[-2].split(",")+i[-1].split(",")

					readNamesFileTxt = miniBamT.replace("miniBam.bam", "_readNameReadsSpanningRearrangement.txt")
					readNamesFile = open(readNamesFileTxt, "w")
					readNamesFile.write("\n".join(readsSpanningRearrangement))
					readNamesFile.close()
					subprocess.call(pathToSamtools+"samtools view -@ "+threadsForSamtools+" -h -b -N "+readNamesFileTxt+" -o "+miniBamT.replace("miniBam.bam", "miniBam_readsSpanningRearrangement.bam")+" "+miniBamT, shell=True)					
					subprocess.call(pathToSamtools+"samtools index "+miniBamT.replace("miniBam.bam", "miniBam_readsSpanningRearrangement.bam"), shell=True) 
					subprocess.call(pathToSamtools+"samtools mpileup -d 0 -a -A -B -Q "+baseq+fr+chromGene+":"+str(i[z])+"-"+str(i[z+1])+" "+miniBamT.replace("miniBam.bam", "miniBam_readsSpanningRearrangement.bam")+ " > "+miniBamT.replace(".bam", "_output_mpileup_readsSpanningRearrangement.tsv")+" 2>> "+errLogMpileup, shell=True)
					
					if os.stat(miniBamT.replace(".bam", "_output_mpileup_readsSpanningRearrangement.tsv")).st_size == 0:
						O = open(miniBamT.replace(".bam", "_output_mpileup_readsSpanningRearrangement.tsv"), "w")
						for missingPos in range(i[z], i[z+1]+1):
							O.write("%s\t%s\tNA\t0\tNA\tNA\n" %(chromGene, str(missingPos)))
						O.close()
					
					# empty mpileup to iterate (to match structure in round 1)
					O = open(miniBamT.replace(".bam", "_output_mpileup_fake.tsv"), "w")
					for missingPos in range(i[z], i[z+1]+1):
						O.write("%s\t%s\tNA\t0\tNA\tNA\n" %(chromGene, str(missingPos)))
					O.close()

					# if seq is amplicon, store input V sequence to keep previously recovered nucleotides
					if seq == "amplicon" and z == 7:
						if round == "second":
							vSeqIn = i[12]
							vSeqGermIn = i[13]
						else: # third
							vSeqIn = i[13]
							vSeqGermIn = i[14]
						if i[1] == "Inversion1" or i[1] == "Inversion2":
							vSeqIn = ''.join(complement[base] for base in reversed(vSeqIn))
							vSeqGermIn = ''.join(complement[base] for base in reversed(vSeqGermIn))
						tumorSeqBK = {}
						germSeqBK = {}
						nuc = int(i[z+1])-int(i[z]) + 1
						for c in range(nuc):
							tumorSeqBK[int(i[z])+c] = [vSeqIn[c]]
							germSeqBK[int(i[z])+c] = [vSeqGermIn[c]]
				
				# check if mpileup result with all reads
				if os.stat(miniBamT.replace(".bam", "_output_mpileup.tsv")).st_size != 0:
					
					# Normal seq:
					if bamN is not None and pairedMode == "paired":
						subprocess.call(pathToSamtools+"samtools mpileup -d 0 -B -Q "+baseq+fr+chromGene+":"+str(i[z])+"-"+str(i[z+1])+" "+miniBamN+ " > "+miniBamN.replace(".bam", "_output_mpileup.tsv")+" 2>> "+errLogMpileup, shell=True)
						
						normal = open(miniBamN.replace(".bam", "_output_mpileup.tsv"), "r")
						wild = {} # normal patient sequence

						sq = int(i[z]) # starts at 1st position interval
						
						passar = 0

						for h in normal:
							w = h.rstrip("\n").split("\t")
							if passar == 0:
								g = 0 # saves possible snps in patient
								dna = {} # save possible snps 
								w[4] = w[4].replace(",", w[2]).replace(".", w[2]) # change reference nucleotide
								w[4] = w[4].upper() # convert all nucleotides to uppercase 
								
								while g < len(w[4]): # look for ACGT in each position
									if w[4][g].isdigit():
										if w[4][g+1].isdigit(): s = 2 # insertion/deletion of > 9 bases
										else: s = 1 # insertion/deletion of < 10 bases
											
										prev = w[4][g-1] # previous shows + for insertions or - for deletions
										now = int(w[4][g:g+s]) # current number of nucleotides being added or removed
										indels = (w[4][g+s:g+s+now]) # nucleotides being added or removed
										
										if prev == "+": # insertion
											ins = w[4][g-2]+"["+indels+"]"
											if ins not in dna: # appends inserted region to dictionary
												dna[ins] = 1  
												dna[w[4][g-2]] -= 1
											else: # adds an occurrence in dictionary
												dna[ins] += 1
												dna[w[4][g-2]] -= 1
											g += now + s # jump as many positions as number shows (number shows nucleotides inserted) plus "s"
											
										elif prev == "-": # deletion
											dele = w[4][g-2]+"("+indels+")"
											if dele not in dna: # appends deleted region to dictionary
												dna[dele] = 1
												dna[w[4][g-2]] -= 1
											else: # adds an occurrence in dictionary
												dna[dele] += 1
												dna[w[4][g-2]] -= 1
											g += now + s # jump as many positions as number shows (number shows nucleotides deleted) plus "s"
									
									elif w[4][g] in ("A", "C", "G", "T"):
										if w[4][g] not in dna: # appends mutation to dictionary
											dna[w[4][g]] = 1
										else: # adds an occurrence in dictionary
											dna[w[4][g]] += 1
										g += 1
									
									else: # ^, $, N, etc.
										if w[4][g] == "^": # to include first base (^6A) (6=ASCI quality; A base of interest)
											g += 2
										else: # last base in read is encoded A$ => A is kept in the previous round, here skip $
											g += 1
								
								dna = sorted(dna.items(), key=lambda dna: dna[1], reverse = True) # from high to low, the possible mutations
								
								tp = []
								if len(dna) > 0:
									for c in dna:
										if int(w[3]) >= depth and c[1] >= altDepth and (c[1]/int(w[3])) >= vafCutoffNormal: # consider base if position depth >= depth, mut count >= altDepth and vaf mutation >= vafCutoffNormal
											tp.append(c[0])
								
								if len(tp) == 0:
									tp.append(w[2]) 

								
								if sq == int(w[1]): # positions in interval with info
									wild[sq] = tp # append mutation to sequence instead of reference nucleotide
									
								else: # non existing positions
									while sq < int(w[1]):
										wild[sq] = ["N"]
										sq += 1
									wild[sq] = tp
									
								if len([True for x in wild[sq] if "(" in x]) > 0:
									wild[sq] = [x for x in wild[sq] if "(" in x]
									passar = len(wild[sq][0]) - 3  # ex: A(CA) -> 3 == A()
								
							else:
								passar -= 1
								
							sq += 1				
						
						# adjust length normal seq if nucleotides are missing
						while sq <= int(i[z+1]):
							wild[sq] = ["N"]
							sq += 1
							
						normal.close()
					
					else:
						wild = {}
						nuc = int(i[z+1])-int(i[z]) + 1
						for c in range(nuc):
							wild[int(i[z])+c] = ["N"]
					
					# Tumor seq:				
					current = open(miniBamT.replace(".bam", "_output_mpileup.tsv"), "r")
					if round == "first":
						currentJV = open(miniBamT.replace(".bam", "_output_mpileup_readsSpanningJV.tsv"), "r")
						currentBreak = open(miniBamT.replace(".bam", "_output_mpileup_readsSpanningBreak.tsv"), "r")
					else: # "second" or "third":
						currentJV = open(miniBamT.replace(".bam", "_output_mpileup_fake.tsv"), "r") # fake mpileup just for iteration
						currentBreak = open(miniBamT.replace(".bam", "_output_mpileup_readsSpanningRearrangement.tsv"), "r") # readnames spanning rearrangement
					
					tumSeq = []
					normSeq = []
					countMuts = 0
					countMutsPhased = 0
					countSNPsNotPhased = 0
					mutPhased = "" # to keep last mutation "chr:pos-pos_nuc" to phase on the fly based on last mutation
					mutPhasedDone = ""
					mutOneRead = 0 # to count number of mutations supported by only one read

					passar = 0 # used to jump sequences if there is a deletion in tumor sequence
					sq = int(i[z])
					for j, j2, j3 in zip(current, currentJV, currentBreak):
						phaseDepth = 0
						inPhase = ""
						
						# phased with exact breakpoints
						if int(j3.rstrip("\n").split("\t")[3]) >= depth: 
							v = j3.rstrip("\n").split("\t")
							phaseDepth = int(v[3])
							vafCutoffToUse = float(vafCutoff[0]) 
							inPhase = "yes"

						# phased with J-V reads
						elif int(j2.rstrip("\n").split("\t")[3]) >= depth: 
							v = j2.rstrip("\n").split("\t")
							phaseDepth = int(v[3])
							vafCutoffToUse = float(vafCutoff[0]) 
							inPhase = "yes"

						# phase reads in V gene if second round using last mutation phased
						if z == 7 and round == "second":
							depthMutPhase = 0
							if mutPhased != "": # try to phase reads by last mutation (if any)
								if mutPhased != mutPhasedDone: # if mut not phased yet... do it
									mutPhasedDone = mutPhased
									readNamePhaseMutTemp = []
									readNamePhaseMut = []
									subprocess.call(pathToSamtools+"samtools mpileup -d 0 -B --output-QNAME -Q "+baseq+fr+mutPhased.split("_")[0]+" "+miniBamT+ " > "+miniBamT.replace(".bam", "_output_mpileup_MutPhased.tsv")+" 2>> "+errLogMpileup, shell=True) # get reads spaining last mutation
									MUTPHASE = open(miniBamT.replace(".bam", "_output_mpileup_MutPhased.tsv"), "r")
									for mutPhaseLine in MUTPHASE:
										vp = mutPhaseLine.rstrip("\n").split("\t")
										vp[4] = vp[4].replace(",", vp[2]).replace(".", vp[2]) # change reference nucleotide
										vp[4] = vp[4].upper() # convert all nucleotides to uppercase
										k = 0 # shows the number of sequences which have to be jumped
										readIndex = 0 # keep the index (ie order/position) of the read
										while k < len(vp[4]): # look for ACGT in each position
											if vp[4][k].isdigit():
												if vp[4][k+1].isdigit(): s = 2 # insertion/deletion of > 9 bases
												else: s = 1 # insertion/deletion of < 10 bases

												prev = vp[4][k-1] # previous shows + for insertions or - for deletions
												now = int(vp[4][k:k+s]) # current number of nucleotides being added or removed
												indels = (vp[4][k+s:k+s+now]) # nucleotides being added or removed
												k += now + s # jump as many positions as number shows (number shows nucleotides inserted/deleted) plus "s"
											
											elif vp[4][k] in ("A", "C", "G", "T", "N", "*", "#", ">", "<"):
												if vp[4][k] == mutPhased.split("_")[1]:
													if vp[6].split(",")[readIndex] not in readPairsRemovedFromPhased:
														readNamePhaseMutTemp.append(vp[6].split(",")[readIndex]) # append readname to readNamePhaseMutTemp if it has seen the mutation
												else:
													readPairsRemovedFromPhased.append(vp[6].split(",")[readIndex]) # append readname to readPairsRemovedFromPhased if it has not seen the mutation
												k += 1
												readIndex += 1 # sum 1 read position
											
											else: # ^ or $
												if vp[4][k] == "^": # to include first base (^6A) (6=ASCI quality; A base of interest)
													k += 2
												else: # last base in read is encoded A$ => A is kept in the previous round, here skip $
													k += 1
									MUTPHASE.close()

									# make sure phased reads do not extend beyond the breakpoint of the V gene
									if len(readNamePhaseMutTemp) > 0:
										breakInV = i[8] if (GENE in ["IGL", "TRA", "TRD"] or (GENE == "TRB" and i[1] == "Deletion") or (GENE == "IGK" and i[1] == "Inversion2")) else i[7]

										readNamesFileTxt = miniBamT.replace("miniBam.bam", "_readNamePhaseMut.txt")
										readNamesFile = open(readNamesFileTxt, "w")
										readNamesFile.write("\n".join(readNamePhaseMutTemp))
										readNamesFile.close()
										subprocess.call(pathToSamtools+"samtools view -@ "+threadsForSamtools+" -N "+readNamesFileTxt+" "+miniBamT+" "+mutPhased.split("_")[0]+" > "+miniBamT.replace("miniBam.bam", "miniSam_readNameMutPhased.sam"), shell=True)
										
										readNamePhaseMut_sam = open(miniBamT.replace("miniBam.bam", "miniSam_readNameMutPhased.sam"), "r")
										for rr in readNamePhaseMut_sam:
											rrList = rr.rstrip("\n").split("\t")
											if rrList[0] in readPairsRemovedFromPhased: 
												continue
											elif GENE in ["IGL", "TRA", "TRD"] or (GENE == "TRB" and i[1] == "Deletion") or (GENE == "IGK" and i[1] == "Inversion2"): 
												rrSplit = re.findall(r'[A-Za-z]|[0-9]+', rrList[5])
												rrSplitTwo = [rrSplit[x:x+2] for x in range(0, len(rrSplit),2)] 
												lastPositionRead = int(rrList[3]) + sum([int(rrSub[0]) for rrSub in rrSplitTwo if "M" in rrSub or "D" in rrSub]) - 1
												if lastPositionRead > breakInV:
													readPairsRemovedFromPhased.append(rrList[0])
											else:
												if int(rrList[3]) < breakInV:
													readPairsRemovedFromPhased.append(rrList[0])
											
										readNamePhaseMut_sam.close()

										readNamePhaseMut = [rNameTmp for rNameTmp in readNamePhaseMutTemp if rNameTmp not in readPairsRemovedFromPhased]
										readsPhasedUsingMutations = readsPhasedUsingMutations + readNamePhaseMut # append to readsPhasedUsingMutations
									
									# make bam with reads phased with last mutation
									readNamesFileTxt = miniBamT.replace("miniBam.bam", "_readNamePhaseMut.txt")
									readNamesFile = open(readNamesFileTxt, "w")
									readNamesFile.write("\n".join(readNamePhaseMut))
									readNamesFile.close()
									subprocess.call(pathToSamtools+"samtools view -@ "+threadsForSamtools+" -h -b -N "+readNamesFileTxt+" -o "+miniBamT.replace("miniBam.bam", "miniSam_readNameMutPhased.bam")+" "+miniBamT, shell=True)
									subprocess.call(pathToSamtools+"samtools index "+miniBamT.replace("miniBam.bam", "miniSam_readNameMutPhased.bam"), shell=True) 
								
								# mpileup for the ongoing position only using bam with reads phased with last mutation
								subprocess.call(pathToSamtools+"samtools mpileup -d 0 -a -A -B -Q "+baseq+fr+chromGene+":"+j.rstrip("\n").split("\t")[1]+"-"+j.rstrip("\n").split("\t")[1]+" "+miniBamT.replace("miniBam.bam", "miniSam_readNameMutPhased.bam")+ " > "+miniBamT.replace(".bam", "_output_mpileup_readNameMutPhased.tsv")+" 2>> "+errLogMpileup, shell=True)
								currentJV_MUTPHASE = open(miniBamT.replace(".bam", "_output_mpileup_readNameMutPhased.tsv"), "r")
								for mutPhaseLine in currentJV_MUTPHASE:
									depthMutPhase = int(mutPhaseLine.rstrip("\n").split("\t")[3]) # keep at depthMutPhase
									j1 = mutPhaseLine # keep at j1
								currentJV_MUTPHASE.close()
								
								# check if enough phased-on-the-fly reads
								if depthMutPhase >= depth and depthMutPhase > phaseDepth:
									v = j1.rstrip("\n").split("\t")
									vafCutoffToUse = float(vafCutoff[0])
									inPhase = "yes"
						
						# if not phased...
						if inPhase == "":
							v = j.rstrip("\n").split("\t")
							vafCutoffToUse = float(vafCutoff[1])
							inPhase = "no"
						
						v[4] = v[4].replace(",", v[2]).replace(".", v[2]) # we change reference nucleotide
						v[4] = v[4].upper() # convert all nucleotides to uppercase
						
						# check if missing positions
						if int(v[3]) == 0:
							if passar == 0:
								# if indel in normal, we consider the indel independently of the tumor seq
								if len([True for x in wild[sq] if "(" in x or "[" in x]) > 0:
									if len([True for x in wild[sq] if "(" in x]) > 0:
										snps = [x for x in wild[sq] if "(" in x]
										passar = len(snps[0]) - 3  # ex: A(CA) -> 3 == A()
									else:
										snps = [x for x in wild[sq] if "]" in x]
									tumSeq.append(snps[0])
									normSeq.append(snps[0])
								else:
									if round != "first" and seq == "amplicon" and z == 7:
										tumSeq.append(tumorSeqBK[sq][0])
										normSeq.append(germSeqBK[sq][0])
									else:
										tumSeq.append("N")
										normSeq.append(wild[sq][0])
							else:
								passar -= 1
								
						elif passar == 0:
							
							newmut = {} # diccionary with nucleotide:n.ocurrences
							k = 0 # shows the number of sequences which have to be jumped
							
							while k < len(v[4]): # look for ACGT in each position
								if v[4][k].isdigit():
									if v[4][k+1].isdigit(): s = 2 # insertion/deletion of > 9 bases
									else: s = 1 # insertion/deletion of < 10 bases
									
									prev = v[4][k-1] # previous shows + for insertions or - for deletions
									now = int(v[4][k:k+s]) # current number of nucleotides being added or removed
									indels = (v[4][k+s:k+s+now]) # nucleotides being added or removed
									
									if prev == "+": # if it is insertion
										ins = v[4][k-2]+"["+indels+"]"
										if ins not in newmut: # appends inserted region to dictionary
											newmut[ins] = 1  
											newmut[v[4][k-2]] -= 1
										else: # adds an occurrence in dictionary
											newmut[ins] += 1
											newmut[v[4][k-2]] -= 1
										k += now + s # jump as many positions as number shows (number shows nucleotides inserted) plus "s"
										
									elif prev == "-": # if it is deletion
										dele = v[4][k-2]+"("+indels+")"
										if dele not in newmut: # appends deleted region to dictionary
											newmut[dele] = 1
											newmut[v[4][k-2]] -= 1
										else: # adds an occurrence in dictionary
											newmut[dele] += 1
											newmut[v[4][k-2]] -= 1
										k += now + s # jump as many positions as number shows (number shows nucleotides deleted) plus "s"
								
								elif v[4][k] in ("A", "C", "G", "T"):
									if v[4][k] not in newmut: # appends mutation to dictionary
										newmut[v[4][k]] = 1
									else: # adds an occurrence in dictionary
										newmut[v[4][k]] += 1
									k += 1
								
								else: # ^, $, N, etc.
									if v[4][k] == "^": # to include first base (^6A) (6=ASCI quality; A base of interest)
										k += 2
									else: # last base in read is encoded A$ => A is kept in the previous round, here skip $
										k += 1
							
							newmut = sorted(newmut.items(), key=lambda newmut: newmut[1], reverse = True) # from high to low, the possible mutations
							
							snps = []
							other = []
							
							# check if nucleotide in normal seq, if not add the one from reference
							if wild[sq] == ["N"]:							
								wild[sq] = [v[2].upper()]
							
							# if nucleotides in tumor newmut:
							if len(newmut) > 0: 
								for c in newmut:
									# min coverage >= depth, mut count >= altDepth 
									if int(v[3]) >= depth and c[1] >= altDepth: 
										# min VAF corrected by tumorPurity (if available and not phased reads) > vafCutoff OR it is SNP in phased reads (to allow ref reads in low coverage phased analyses)
										if ( c[1]/int(v[3])/(1 if inPhase == "yes" else tumorPurity) ) >= vafCutoffToUse or ( inPhase == "yes" and c[0] in wild[sq] ): 
											if c[0] in wild[sq]:
												snps.append(c[0])
											else:
												other.append(c[0])
												# count if the mutation is supported by one read
												if len(other) == 1 and c[1] == 1: 
													mutOneRead += 1
								
								# if indel in normal, we consider the indel independently of the tumor seq
								if len([True for x in wild[sq] if "(" in x or "[" in x]) > 0:
									if len([True for x in wild[sq] if "(" in x]) > 0:
										snps = [x for x in wild[sq] if "(" in x]
										passar = len(snps[0]) - 3  # ex: A(CA) -> 3 == A()
									else:
										snps = [x for x in wild[sq] if "]" in x]
									tumSeq.append(snps[0])
									normSeq.append(snps[0])
								
								# no indel in normal, and mutation in tumor
								elif len(other) > 0:
									tumSeq.append(other[0])
									if "(" in other[0] or "[" in other[0]:
										normSeq.append(other[0])
										if "(" in other[0]:
											passar = len(other[0]) - 3  # ex: A(CA) -> 3 == A()
									else:
										normSeq.append(wild[sq][0])
									
									# update count phased muts
									countMuts += 1
									if inPhase == "yes":
										countMutsPhased += 1
										# add info to mutPhased to be used in subsequent phasing if round is second and phasing needs to be done
										if round == "second":
											snpInPopulation = "no"
											if bamN is None or pairedMode == "unpaired": # check in SNPs list only if no normal BAM or unpaired mode
												SNPs = open(snps_file, "r")
												for iSnp in SNPs:
													iSnpLst = iSnp.rstrip("\n").split("\t")
													if v[0].replace("chr", "") == iSnpLst[0] and v[1] == iSnpLst[1] and other[0] == iSnpLst[3]:
														snpInPopulation = "yes"
														break
												SNPs.close()
											if snpInPopulation == "no":
												mutPhased = v[0]+":"+v[1]+"-"+v[1]+"_"+other[0] 
									
								# no indel in normal neither mutation in tumor, add first SNP in both:
								elif len(snps) > 0:
									tumSeq.append(sorted(snps)[0]) # if more than one nucleotide, sort them and keep the 1st for reproducibility and because coverage is not informative of the rearranged snp if not phased
									normSeq.append(sorted(snps)[0])
									if len(snps) > 1: countSNPsNotPhased += 1 # count number of SNPs not phased
								
								# else... add N:
								else:
									tumSeq.append("N")
									normSeq.append(wild[sq][0])
							
							# in case no newmut... 
							else:
								if len([True for x in wild[sq] if "(" in x or "[" in x]) > 0:
									if len([True for x in wild[sq] if "(" in x]) > 0:
										snps = [x for x in wild[sq] if "(" in x]
										passar = len(snps[0]) - 3  # ex: A(CA) -> 3 == A()
									else:
										snps = [x for x in wild[sq] if "]" in x]
									tumSeq.append(sorted(snps)[0]) # if more than one nucleotide, sort them and keep the 1st for reproducibility and because coverage is not informative of the rearranged snp if not phased
									normSeq.append(sorted(snps)[0])
									if len(snps) > 1: countSNPsNotPhased += 1 # count number of SNPs not phased
								else:
									tumSeq.append("N")
									normSeq.append(wild[sq][0])
						
						# passar != 0
						else:
							passar -= 1
						
						sq += 1
					
					current.close()
					currentJV.close()
					currentBreak.close()
					
					# Adjust length tumor and normal if nucleotides missing:
					if sq-1 < int(i[z+1]):
						tumSeq.extend("N" * (int(i[z+1])-sq-1))
						
					if len(normSeq) < len(tumSeq):
						normSeq.extend("N" * (len(tumSeq)-len(normSeq)))
					
					if z == 7: # V
						temporary.append("".join(tumSeq)) # V seq tumor
						temporary.append("".join(normSeq)) # V seq normal
						phaseInfo = str(countMutsPhased)+"/"+str(countMuts)+" - "+str(mutOneRead)+" - "+str(countSNPsNotPhased)
					
					else: # J -> add sequence of J and NA for D
						temporary.append(''.join(tumSeq)) # J seq
						temporary.append("") # empty for D seq
						
				else: # if no mpileup file
					temporary.extend(["NA"]*2) 
				
			else: # if no start-end
				temporary.extend(["NA"]*2) 
			
		# Add phaseInfo to temporary and append/update temporary to i
		temporary.append(phaseInfo)
		if round == "first": 
			i[-1:-1] = temporary
		elif round == "second":
			i[10] = temporary[0]
			i[12] = temporary[2]
			i[13] = temporary[3]
			i[14] = temporary[4]
			readsPhasedUsingMutationsString = ",".join(rName for rName in set(readsPhasedUsingMutations) if rName not in readsAlreadyUsed and rName not in readPairsRemovedFromPhased)
			countReadsPhasedMutations = 0 if readsPhasedUsingMutationsString == "" else len(readsPhasedUsingMutationsString.split(","))
			i[10:10] = [countReadsPhasedMutations]
			i.append(readsPhasedUsingMutationsString)
		else: # "third"
			i[11] = temporary[0]
			i[13] = temporary[2]
			i[14] = temporary[3]
			i[15] = temporary[4]

		if "Kde" not in i[0] and "RSS" not in i[0] and i[4] != "NA":
			# do reverse complement if needed
			if i[1] == "Inversion1" or i[1] == "Inversion2": # TRB inversion1 V in strand positive or IGK inversion2 V in strand negative
				idxSum = 0 if round == "first" else 1 # 1 if second or third round
				i[12+idxSum] = ''.join(complement[base] for base in reversed(i[12+idxSum]))
				i[13+idxSum] = ''.join(complement[base] for base in reversed(i[13+idxSum]))

			# update complete seq if second or third round
			if round in ["second", "third"]:
				if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
					totseqW = i[11]+i[12]+i[13]
				else:
					totseqW = i[13]+i[12]+i[11]
				totseqW = re.sub("\(.*?\)", "", totseqW.replace("[", "").replace("]", ""))
				i[16] = totseqW
	
	return(information)
	
def createConsensusD(DseqTemp, GENE, i, Dseqs):
	
	geneNames = i[0]
	
	DseqConsensus = list()
	
	for c in range(0, len(DseqTemp[0])):
		DseqConsensus.append(Counter([item[c] for item in DseqTemp]).most_common(1)[0][0])
	DseqConsensus = "".join(DseqConsensus)
	
	# check if insertion occurs at breaks (then the insertion is added at J/V in the mpileup and also present in the "D" sequence...remove it from "D"):
	if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
		j = i[10]
		v = i[12]
	else:
		j = i[12] # it is V
		v = i[10] # it is J
	
	if j[-1] == "]":
		if DseqConsensus.startswith(j[max([a.start()+1 for a in re.finditer("\[", j)]):-1]):
			DseqConsensus = DseqConsensus[len(j)-2-max([a.start()+1 for a in re.finditer("\[", j)])+1:]
			
	if v[0] == "[":
		if DseqConsensus.endswith(v[1:min([a.start()+1 for a in re.finditer("\]", v)])-1]):
			DseqConsensus = DseqConsensus[:len(DseqConsensus)-min([a.start()+1 for a in re.finditer("\]", v)])+2]
	
	# if IGH or TRB or TRD, and not partial rearrangement (J-D), check D gene and update geneNames:
	if GENE in ["IGH", "TRB", "TRD"] and i[0].split(" - ")[0][3] != "D" and i[0].split(" - ")[1][3] != "D":
		sqs = open(Dseqs, "r")
		t = []
		for sqsLine in sqs:
			w = sqsLine.rstrip("\n").split("\t")
			score = smithwaterman(DseqConsensus, w[1])
			t.append([w[0], score])
		
		dGeneName = sorted(t, key=operator.itemgetter(1), reverse = True)[0][0]
		geneNames = " - ".join([i[0].split(" - ")[0], dGeneName, i[0].split(" - ")[1]]) # update J-V to J-D-V

	# add N-D-N / N
	return(geneNames, DseqConsensus)		

def fillInGapBetweenDandV(dSeq, potentialReadNames, i, bedFile, minimumNumberOfNucleotidesSoft, refGenome, miniBamT, pathToSamtools, errLogMpileup):
    
	# define minimum length to consider for alignment
	minimumNumberOfNucleotidesSoftAmplicon = minimumNumberOfNucleotidesSoft + 2
	windowInBreak = 5

	# get info from i
	vGene = i[0].split(" - ")[1]
	vBreak = i[7]
	vEnd = i[8]
	readNamesList = potentialReadNames

	# get start V gene brom bedFile
	VDJ = open(bedFile, "r")
	for k in VDJ:
		v = k.rstrip("\n").split("\t")
		if vGene == v[3]:
			chrom = v[0]
			startV = int(v[1]) - windowInBreak
			break
	VDJ.close()
	
	# subset miniBamT with reads of rearrangement
	fillInReadsTxt = miniBamT.replace(".bam", "_fillInReads.txt")
	miniBamTreadsFillIn = miniBamT.replace(".bam", "_fillInReads.bam")
	O = open(fillInReadsTxt, "w")
	for rn in readNamesList:
		O.write(rn+"\n")
	O.close()
	comm = pathToSamtools+"samtools view -b -h -N "+fillInReadsTxt+" -o "+miniBamTreadsFillIn+" "+miniBamT
	subprocess.call(comm, shell=True)
	subprocess.call(pathToSamtools+"samtools index "+miniBamTreadsFillIn, shell=True)
	
	# mpileup to get reference sequence
	comm = pathToSamtools+"samtools mpileup -d 0 -a -A -f "+refGenome+" -r "+str(chrom)+":"+str(startV)+"-"+str(vEnd)+" "+miniBamTreadsFillIn+" > "+miniBamTreadsFillIn.replace(".bam", "_mpileup.tsv")+" 2>> "+errLogMpileup
	subprocess.call(comm, shell=True)
	
	# get potential gap sequence
	MPILEUP = open(miniBamTreadsFillIn.replace(".bam", "_mpileup.tsv"), "r")
	potentialGap = ""
	for mLine in MPILEUP:
		mList = mLine.rstrip("\n").split("\t")
		if int(mList[1]) < vBreak+minimumNumberOfNucleotidesSoftAmplicon:
			potentialGap = potentialGap+mList[2]
		else:
			break
	MPILEUP.close()
	
	# align dSeq against potentialGap to get missingNucleotides
	vStartNucleotidesGermline = "NotFound"
	vStartNucleotidesTumor = "NotFound"
	gapFound = "NotFound"
	nGapFound = "NotFound"
	if potentialGap != "":
		aligner = PairwiseAligner()
		# force stitching without internal gap (potentialGap [a]: ---xx[nnn], dSeq [b]: nnnxx----)
		aligner.mode = "global"
		aligner.match_score = 2
		aligner.mismatch_score = -1
		aligner.open_gap_score = -1000
		aligner.extend_gap_score = -1000
		aligner.left_open_gap_score = 0
		aligner.right_open_gap_score = 0
		aligner.left_extend_gap_score = 0
		aligner.right_extend_gap_score = 0
		dAlignment = aligner.align(potentialGap, dSeq)
		numAlignments = 0
		try: numAlignments = len(dAlignment)
		except OverflowError: numAlignments = 1000 
		if numAlignments > 0:
			# get alignment that ends at the minimum position (closer to the start position of the V gene) and check for stitching
			stitchingOK = "no"
			if str(dAlignment[0]).startswith("target"):
				bestFirstAlignment = min(dAlignment, key=lambda aln: aln.aligned.shape[1])
				alignA = [el for el in str(bestFirstAlignment).split("\n")[0].split(" ") if el != ""][-2]
				alignS = [el for el in str(bestFirstAlignment).split("\n")[1].split(" ") if el != ""][-2]
				alignB = [el for el in str(bestFirstAlignment).split("\n")[2].split(" ") if el != ""][-2]
				if not alignA.endswith("-") and not alignB.startswith("-"):
					stitchingOK = "yes"
			else:
				bestFirstAlignment = min(dAlignment, key=lambda aln: aln.aligned[1][-1][1])
				alignA = str(bestFirstAlignment).split("\n")[0].replace(" ", "-")
				alignS = str(bestFirstAlignment).split("\n")[1].replace(" ", "-")
				alignB = str(bestFirstAlignment).split("\n")[2].replace(" ", "-")
				if not alignA.endswith("-") and not alignB.startswith("-"):
					stitchingOK = "yes"

			if stitchingOK == "yes":
				startIndx = alignS.find("|")
				endIndx = alignS.rfind("|")+1
				numMatches = sum(1 for a, b in zip(alignA[startIndx:endIndx], alignB[startIndx:endIndx]) if a == b and a != "-" and b != "-")
				lenAlignment = endIndx - startIndx
				ratioMatchesMissmatches = numMatches/lenAlignment
				# force minimum length and matches in alignment based on position
				empySeqBstart = startIndx - (len(re.match(r"^-+", alignA).group()) if re.match(r"^-+", alignA) else 0)
				condiA = empySeqBstart < 5+windowInBreak and lenAlignment >= minimumNumberOfNucleotidesSoft and ratioMatchesMissmatches >= 0.8
				condiB = empySeqBstart >= 5+windowInBreak and lenAlignment >= minimumNumberOfNucleotidesSoftAmplicon and ratioMatchesMissmatches >= 0.85
				if condiA or condiB:
					gapFound = alignA[endIndx:-minimumNumberOfNucleotidesSoftAmplicon]
					nGapFound = "N"*len(gapFound)
					vStartNucleotidesGermline = alignA[startIndx:endIndx]
					vStartNucleotidesTumor = alignB[startIndx:endIndx]
	
	return(vStartNucleotidesGermline, vStartNucleotidesTumor, gapFound, nGapFound)
	
def getDsequence(information, annot_table_JV, GENE, Dseqs, minimumNumberOfNucleotidesSoft, seq, bedFile, refGenome, miniBamT, pathToSamtools, errLogMpileup):
	
	readsAlreadyRecovered = [] # list to append readNames when already recovered
	toAddInInformation = [] # list to append to Information if same D with same length
	ndnBasesToAdd = 5 if GENE in ["IGH", "TRB", "TRD"] else 0 # number to count for N-D-N nucleotides while recovering J-V splits (works together with minimumNumberOfNucleotidesSoft) 

	# Round 1: Recover V reads, sort, and recover J reads
	for geneRound in [1, 2]:
		for i in information:
			
			# Add two elements in i for readNames rescued at J and V
			if geneRound == 1:
				i.append("")
				i.append("")
			# Get soft clipped start/end J-V :  
			if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
				breakJ = int(i[5])
				breakV = int(i[7]) if i[1] == "Deletion" else int(i[8]) # for IGK Inversion2
			else:
				breakV = int(i[4]) # it is J
				breakJ = int(i[8]) if i[1] == "Deletion" else int(i[7]) # it is V  # for TRB Inversion1
			
			DseqTemp = []
			nonMatchedDsInJ = [] # list to append tupple of readName - D seq not matched in J

			if not "IGKKde" in i[0] and not "IGKRSS" in i[0]:
				ANNOT_TABLE_JV = open(annot_table_JV, "r")
				for j in ANNOT_TABLE_JV:
					w = j.rstrip("\n").split("\t")
					if w[0] in readsAlreadyRecovered: continue

					if w[11].startswith("split"):

						# if information last value J and first value V coincide with w split values
						if breakJ == int(w[12].replace("NA", "0")) and breakV == int(w[13].replace("NA", "0")):
							
							# split or split-insertSize with SA info in alignment
							if w[10] != "NA":
								split = re.findall(r'[A-Za-z]|[0-9]+', w[5])
								cigar1 = [split[x:x+2] for x in range(0, len(split), 2)]
								split = re.findall(r'[A-Za-z]|[0-9]+', w[10].split(",")[3])
								cigar2 = [split[x:x+2] for x in range(0, len(split), 2)]
								# MS cigar
								if min([count for count, item in enumerate(cigar1) if "M" in item]) < min([count for count, item in enumerate(cigar1) if "S" in item]):
									mStart = sum([ int(x[0]) if x[1] in ["M", "I"] else 0 for x in cigar1 ]) # we add numbers previous to M and I
									if w[20] == "Deletion": subcigar2 = cigar2[:min([i for i in range(len(cigar2)) if cigar2[i][1] == "M"])] # limit cigar2 to first M
									else: subcigar2 = cigar2
									mEnd = sum([ int(x[0]) if x[1] == "S" else 0 for x in subcigar2 ]) # we add numbers previous to S
								# SM cigar
								else:
									mStart = sum([ int(x[0]) if x[1] in ["M", "I"] else 0 for x in cigar2 ]) # we add numbers previous to M and I
									if w[20] == "Deletion": subcigar1 = cigar1[:min([i for i in range(len(cigar1)) if cigar1[i][1] == "M"])] # limit cigar1 to first M
									else: subcigar1 = cigar1
									mEnd = sum([ int(x[0]) if x[1] == "S" else 0 for x in subcigar1 ]) # we add numbers previous to S
								
								DseqTemp.append(w[9][mStart:mEnd]) # we analyse from M,I+seq until seq-everything but S
								readsAlreadyRecovered.append(w[0])

							# split-insertSize and split-insertSize without SA info in alignment
							else: 
								# get soft clipped in J and V breaks
								jUnmap = ""
								vUnmap = ""
								comm = pathToSamtools+"samtools view "+miniBamT+" | grep "+w[0]+" > "+miniBamT.replace(".bam", "_spliInsert2reads.txt")
								subprocess.call(comm, shell=True)
								
								SPLIT2READS = open(miniBamT.replace(".bam", "_spliInsert2reads.txt"), "r")
								for split2read in SPLIT2READS:
									split2readList = split2read.rstrip("\n").split("\t")
									split = re.findall(r'[A-Za-z]|[0-9]+', split2readList[5])
									cigar1 = [split[x:x+2] for x in range(0, len(split), 2)]
									if not any("M" in cigEl for cigEl in cigar1) or not any("S" in cigEl for cigEl in cigar1): continue
									# MS cigar
									if min([count for count, item in enumerate(cigar1) if "M" in item]) < min([count for count, item in enumerate(cigar1) if "S" in item]):
										mStart = sum([ int(x[0]) if x[1] in ["M", "I"] else 0 for x in cigar1 ]) - 1
									# SM cigar
									else:
										if w[20] == "Deletion": subcigar1 = cigar1[:min([i for i in range(len(cigar1)) if cigar1[i][1] == "M"])] # limit cigar1 to first M
										else: subcigar1 = cigar1
										mEnd = sum([ int(x[0]) if x[1] == "S" else 0 for x in subcigar1 ])
										mStart = 0
									if int(w[12]) == int(split2readList[3])+mStart:
										jUnmap = split2readList[9][mStart+1:]
									elif int(w[13]) == int(split2readList[3])+mStart:
										vUnmap = split2readList[9][:mEnd]
								SPLIT2READS.close()
								
								# align soft clipped to get entire sequence
								if jUnmap != "" and vUnmap != "":
									aligner = PairwiseAligner()
									# force stitching without internal gap
									aligner.mode = "global"
									aligner.match_score = 2
									aligner.mismatch_score = -1
									aligner.open_gap_score = -1000
									aligner.extend_gap_score = -1000
									aligner.left_open_gap_score = 0
									aligner.right_open_gap_score = 0
									aligner.left_extend_gap_score = 0
									aligner.right_extend_gap_score = 0
									unmapAlignment = aligner.align(jUnmap, vUnmap)
									numAlignments = 0
									try: numAlignments = len(unmapAlignment)
									except OverflowError: numAlignments = 1000 
									if numAlignments > 0:
										unmapAlignment = unmapAlignment[0]
										if str(unmapAlignment).startswith("target"):
											alignA = [el for el in str(unmapAlignment).split("\n")[0].split(" ") if el != ""][-2]
											alignS = [el for el in str(unmapAlignment).split("\n")[1].split(" ") if el != ""][-2]
											alignB = [el for el in str(unmapAlignment).split("\n")[2].split(" ") if el != ""][-2]
										else:
											alignA = str(unmapAlignment).split("\n")[0].replace(" ", "-")
											alignS = str(unmapAlignment).split("\n")[1].replace(" ", "-")
											alignB = str(unmapAlignment).split("\n")[2].replace(" ", "-")
										# stitching conditions: condiA (a: xxxxx, b: xxxxx) or condiB (a: nnnxx--, b:---xxnnn) or condiC (a: xx---, b:xxnnn) or condiD (a: nnnxx, b:---xx)
										condiA = not alignA.startswith("-") and not alignA.endswith("-") and not alignB.startswith("-") and not alignB.endswith("-")
										condiB = not alignA.startswith("-") and alignA.endswith("-") and alignB.startswith("-") and not alignB.endswith("-")
										condiC = not alignA.startswith("-") and alignA.endswith("-") and not alignB.startswith("-") and not alignB.endswith("-")
										condiD = not alignA.startswith("-") and not alignA.endswith("-") and alignB.startswith("-") and not alignB.endswith("-")
										if condiA or condiB or condiC or condiD:
											startIndx = alignS.find("|")
											endIndx = alignS.rfind("|")+1
											numMatches = sum(1 for a, b in zip(alignA[startIndx:endIndx], alignB[startIndx:endIndx]) if a == b and a != "-" and b != "-")
											lenAlignment = endIndx - startIndx
											ratioMatchesMissmatches = numMatches/lenAlignment
											# force minimum length and matches in alignment
											if lenAlignment >= minimumNumberOfNucleotidesSoft and ratioMatchesMissmatches >= 0.8:
												unmapNucs = []
												for a, b in zip(alignA, alignB):
													if a != "-" and b == "-":
														unmapNucs.append(a)
													elif a == "-" and b != "-":
														unmapNucs.append(b)
													elif a != "-" and b != "-" and a == b:
														unmapNucs.append(a)
													else:
														unmapNucs.append("N")
												unmapSeq = "".join(unmapNucs)
												# append to DseqTemp
												DseqTemp.append(unmapSeq)
												readsAlreadyRecovered.append(w[0])
						
						# if information last position J:
						elif breakJ == int(w[12].replace("NA", "0")) and w[13] == "NA" and ((geneRound == 1 and GENE in ["IGL", "TRA", "TRB", "TRD"]) or (geneRound == 2 and GENE not in ["IGL", "TRA", "TRB", "TRD"])):
							# force insertSize read close to breakpoint
							if w[11] == "split-insertSize" and abs(breakV - int(w[7])) > 1000: continue 
							# cigar
							split = re.findall(r'[A-Za-z]|[0-9]+', w[5])
							cigar1 = [split[x:x+2] for x in range(0, len(split), 2)]
							# MS cigar
							if min([count for count, item in enumerate(cigar1) if "M" in item]) < min([count for count, item in enumerate(cigar1) if "S" in item]):
								if w[20] == "Deletion": subcigar1 = cigar1[max([i for i in range(len(cigar1)) if cigar1[i][1] == "M"]):] # limit cigar1 to start at last M
								else: subcigar1 = cigar1
								mStart = sum([ int(x[0]) if x[1] == "S" else 0 for x in subcigar1 ]) # we add numbers previous to S
								J = w[9][-mStart:]
							# SM cigar
							else:
								if w[20] == "Deletion": subcigar1 = cigar1[:min([i for i in range(len(cigar1)) if cigar1[i][1] == "M"])] # limit cigar1 to first M
								else: subcigar1 = cigar1
								mEnd = sum([ int(x[0]) if x[1] == "S" else 0 for x in subcigar1 ]) # we add numbers previous to S
								J = w[9][:mEnd]
							
							if len(J) < ndnBasesToAdd + minimumNumberOfNucleotidesSoft: continue

							# Vseq: remove deleted nucleotides, check insertion at first bases, keep insertions not at first base:
							if GENE not in ["IGL", "TRA", "TRB", "TRD"]: 
								vSeq = re.sub("\(.*?\)", "",  i[12])
								if vSeq[0] == "[":
									vSeq = vSeq[min([a.start()+1 for a in re.finditer("\]", vSeq)]):]
							
							else: # it is J in IGL/TRA/TRB/TRD
								vSeq = re.sub("\(.*?\)", "",  i[10])
								if vSeq[-1] == "]":
									vSeq = vSeq[:max([a.start() for a in re.finditer("\[", vSeq)]):]
							
							vSeq = vSeq.replace("[", "").replace("]", "")

							if i[1] == "Inversion1": J = ''.join(complement[base] for base in reversed(J))

							j = ndnBasesToAdd
							while j <= len(J)-minimumNumberOfNucleotidesSoft:
								if vSeq.startswith(J[j:j+minimumNumberOfNucleotidesSoft]):
									DseqTemp.append(J[:j])
									if GENE not in ["IGL", "TRA", "TRB", "TRD"]: 
										i[6] += 1 # count split J
										i[16] = w[0] if i[16] == "" else i[16]+","+w[0]
									else: 
										i[9] += 1 # count split V
										i[17] = w[0] if i[17] == "" else i[17]+","+w[0]
									readsAlreadyRecovered.append(w[0])
									break
								j += 1
							
							# add D-seq to nonMatchedDsInJ if not recovered:
							if w[0] not in readsAlreadyRecovered: 
								nonMatchedDsInJ.append((w[0], J))
						
						# if information first position V:
						elif breakV == int(w[12].replace("NA", "0")) and w[13] == "NA" and ((geneRound == 1 and GENE not in ["IGL", "TRA", "TRB", "TRD"]) or (geneRound == 2 and GENE in ["IGL", "TRA", "TRB", "TRD"])):
							# force insertSize read close to breakpoint
							if w[11] == "split-insertSize" and abs(breakJ - int(w[7])) > 1000: continue 
							# cigar
							split = re.findall(r'[A-Za-z]|[0-9]+', w[5])
							cigar1 = [split[x:x+2] for x in range(0, len(split),2)]
							# MS cigar
							if min([count for count, item in enumerate(cigar1) if "M" in item]) < min([count for count, item in enumerate(cigar1) if "S" in item]):
								if w[20] == "Deletion": subcigar1 = cigar1[max([i for i in range(len(cigar1)) if cigar1[i][1] == "M"]):] # limit cigar1 to start at last M
								else: subcigar1 = cigar1
								mStart = sum([ int(x[0]) if x[1] == "S" else 0 for x in subcigar1 ]) # we add numbers previous to S
								V = w[9][-mStart:]
							# SM cigar
							else:
								if w[20] == "Deletion": subcigar1 = cigar1[:min([i for i in range(len(cigar1)) if cigar1[i][1] == "M"])] # limit cigar1 to first M
								else: subcigar1 = cigar1
								mEnd = sum([ int(x[0]) if x[1] == "S" else 0 for x in subcigar1 ]) # we add numbers previous to S
								V = w[9][:mEnd]
							
							if len(V) < ndnBasesToAdd + minimumNumberOfNucleotidesSoft: continue

							# jSeq: remove deleted nucleotides, check insertion at last bases, keep insertions not at last base:
							if GENE not in ["IGL", "TRA", "TRB", "TRD"]: 
								jSeq = re.sub("\(.*?\)", "",  i[10])
								if jSeq[-1] == "]":
									jSeq = jSeq[:max([a.start() for a in re.finditer("\[", jSeq)]):]
							else:  # it is V in IGL/TRA/TRB/TRD
								jSeq = re.sub("\(.*?\)", "",  i[12])
								if jSeq[0] == "[":
									jSeq = jSeq[min([a.start()+1 for a in re.finditer("\]", jSeq)]):]
									
							jSeq = jSeq.replace("[", "").replace("]", "")

							if i[1] == "Inversion2": V = ''.join(complement[base] for base in reversed(V))
							
							v = len(V) - ndnBasesToAdd
							while v >= minimumNumberOfNucleotidesSoft:
								if jSeq.endswith(V[v-minimumNumberOfNucleotidesSoft:v]):
									DseqTemp.append(V[v:])
									if GENE not in ["IGL", "TRA", "TRB", "TRD"]: 
										i[9] += 1 # count split V
										i[17] = w[0] if i[17] == "" else i[17]+","+w[0]
									else: 
										i[6] += 1 # count split J
										i[16] = w[0] if i[16] == "" else i[16]+","+w[0]
									readsAlreadyRecovered.append(w[0])
									break
								v -= 1
							
							# add D-seq to nonMatchedDsInJ if not recovered:
							if w[0] not in readsAlreadyRecovered: 
								nonMatchedDsInJ.append((w[0], V))

				ANNOT_TABLE_JV.close()
			
			if DseqTemp == []: DseqTemp.append("NotFoundInIter")
			if geneRound == 1:
				i.append(",".join(DseqTemp))
			else:
				i[-1] = i[-1]+","+",".join(DseqTemp)
				i.append(",".join([f"{key};{value}" for key, value in nonMatchedDsInJ]))

		# sort information to account for recovered V reads (geneRound 1) and also J reads (geneRound 2)
		information.sort(key=lambda p: round(p[2]*2 + p[3] + p[6]*2 + p[9]*2, 1), reverse=True)
	
	# Report Ds:
	readsRecoveredNonMatched = []
	for i in information:
		DseqTemp = i[-2].split(",") # get DseqTemp from i
		DseqTemp = [dseq for dseq in DseqTemp if dseq != "NotFoundInIter"] # remove NotFoundInIter element if present
		nonMatchedDsInJ = i[-1].split(",") # get nonMatchedDsInJ from i
		del i[-2:] # delete DseqTemp and nonMatchedDsInJ from i
		AorBdone = "no"

		if len(DseqTemp) > 0:
			
			## A) all possible "D"s have different lengths... keep them all...
			if len(set([len(s) for s in DseqTemp])) == len(DseqTemp): 
				
				countToAdd = 1
				for DseqTempSimple in DseqTemp:
					geneNames, DseqConsensus = createConsensusD([DseqTempSimple], GENE, i, Dseqs)
					if countToAdd < len(DseqTemp):
						iToAddInToAddInInformationlist = i.copy()
						iToAddInToAddInInformationlist[0] = geneNames
						iToAddInToAddInInformationlist[11] = DseqConsensus
						if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
							totseqW = iToAddInToAddInInformationlist[10]+iToAddInToAddInInformationlist[11]+iToAddInToAddInInformationlist[12]
						else:
							totseqW = iToAddInToAddInInformationlist[12]+iToAddInToAddInInformationlist[11]+iToAddInToAddInInformationlist[10]
						totseqW = re.sub("\(.*?\)", "", totseqW.replace("[", "").replace("]", ""))
						iToAddInToAddInInformationlist.insert(-3, totseqW)
						
						toAddInInformation.append(iToAddInToAddInInformationlist)
						
					else:
						i[0] = geneNames
						i[11] = DseqConsensus
						if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
							totseqW = i[10]+i[11]+i[12]
						else:
							totseqW = i[12]+i[11]+i[10]
						totseqW = re.sub("\(.*?\)", "", totseqW.replace("[", "").replace("]", ""))
						i.insert(-3, totseqW)
					
					countToAdd += 1
				
				AorBdone = "yes"
			
			## B) if not, and first and second D lengths have the same number of supporting reads -> keep both!
			elif len(set([len(s) for s in DseqTemp])) >= 2: 
				
				if Counter([len(s) for s in DseqTemp]).most_common(2)[0][1] == Counter([len(s) for s in DseqTemp]).most_common(2)[1][1]:
					
					# get first
					DseqTempSimple = [ss for ss in DseqTemp if len(ss) == Counter([len(s) for s in DseqTemp]).most_common(2)[0][0]] # get Dseqs with the same length	
					geneNames, DseqConsensus = createConsensusD(DseqTempSimple, GENE, i, Dseqs)
					iToAddInToAddInInformationlist = i.copy()
					iToAddInToAddInInformationlist[0] = geneNames
					iToAddInToAddInInformationlist[11] = DseqConsensus
					if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
						totseqW = iToAddInToAddInInformationlist[10]+iToAddInToAddInInformationlist[11]+iToAddInToAddInInformationlist[12]
					else:
						totseqW = iToAddInToAddInInformationlist[12]+iToAddInToAddInInformationlist[11]+iToAddInToAddInInformationlist[10]
					totseqW = re.sub("\(.*?\)", "", totseqW.replace("[", "").replace("]", ""))
					iToAddInToAddInInformationlist.insert(-3, totseqW)
					toAddInInformation.append(iToAddInToAddInInformationlist)
					
					# get second
					DseqTempSimple = [ss for ss in DseqTemp if len(ss) == Counter([len(s) for s in DseqTemp]).most_common(2)[1][0]] # get Dseqs with the same length	
					geneNames, DseqConsensus = createConsensusD(DseqTempSimple, GENE, i, Dseqs)
					i[0] = geneNames
					i[11] = DseqConsensus
					if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
						totseqW = i[10]+i[11]+i[12]
					else:
						totseqW = i[12]+i[11]+i[10]
					totseqW = re.sub("\(.*?\)", "", totseqW.replace("[", "").replace("]", ""))
					i.insert(-3, totseqW)
					
					AorBdone = "yes"
					
			## C) if not A or B, get Dseqs with the same length
			if AorBdone == "no":
				DseqTemp = [ss for ss in DseqTemp if len(ss) == Counter([len(s) for s in DseqTemp]).most_common(1)[0][0]] # get Dseqs with the same length		
				geneNames, DseqConsensus = createConsensusD(DseqTemp, GENE, i, Dseqs)
				i[0] = geneNames
				i[11] = DseqConsensus
				
				if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
					totseqW = i[10]+i[11]+i[12]
				else:
					totseqW = i[12]+i[11]+i[10]
				totseqW = re.sub("\(.*?\)", "", totseqW.replace("[", "").replace("]", ""))
				i.insert(-3, totseqW)
		
		## D) No D... just concatenate sequence
		else:
			if "IGKKde" in i[0] or "IGKRSS" in i[0]:
				totseqW = "NA"
			else:
				# if seq is amplicon, try to recover IGH rearrangement by trying to align 'D' against the germline V sequence before the breakpoint (only IGH due to N-D-N plus SHM)
				if seq == "amplicon" and GENE == "IGH" and nonMatchedDsInJ != ['']:
					# get nonMatched nucleotides for reads not used
					passingNonMatchedDsInJ = []
					listNonMatchedDsInJ = []
					for tupp in nonMatchedDsInJ:
						tuppList = tupp.split(";")
						if tuppList[0] not in readsAlreadyRecovered:
							passingNonMatchedDsInJ.append(tuppList[1]) # sequence
							listNonMatchedDsInJ.append(tupp) # tupp
					if len(passingNonMatchedDsInJ) > 0:
						# get most common nonMatched nucleotides
						DseqNotMatchedMostCommon = Counter(passingNonMatchedDsInJ).most_common(1)[0]
						# get readNames of these nonMatched nucleotides
						readNamesDseqNotMatchedMostCommon = []
						for tupp in listNonMatchedDsInJ:
							tuppList = tupp.split(";")
							if tuppList[1] == DseqNotMatchedMostCommon[0]:
								readNamesDseqNotMatchedMostCommon.append(tuppList[0])
						# try to get missing nucleotides
						vStartNucleotidesGermline, vStartNucleotidesTumor, gapFound, nGapFound = fillInGapBetweenDandV(DseqNotMatchedMostCommon[0], readNamesDseqNotMatchedMostCommon, i, bedFile, minimumNumberOfNucleotidesSoft, refGenome, miniBamT, pathToSamtools, errLogMpileup)
						if gapFound != "NotFound":
							# assign nonMatched (ie D) in 11; count of reads in 6; and readNames in -2 (split J count)
							i[11] = DseqNotMatchedMostCommon[0][:-len(vStartNucleotidesTumor)]
							i[6] = DseqNotMatchedMostCommon[1]
							i[-2] = ",".join(readNamesDseqNotMatchedMostCommon)
							# add missing nucleotides to the begining of V (somatic and germline) and adjust breakpoint of V
							i[12] = vStartNucleotidesTumor+nGapFound+i[12]
							i[13] = vStartNucleotidesGermline+gapFound+i[13]
							i[7] = i[7]-len(vStartNucleotidesTumor)-len(nGapFound)
							# add readNames to readsRecoveredNonMatched
							readsRecoveredNonMatched.extend(readNamesDseqNotMatchedMostCommon)

				if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
					totseqW = i[10]+i[11]+i[12]
				else:
					totseqW = i[12]+i[11]+i[10]
			totseqW = re.sub("\(.*?\)", "", totseqW.replace("[", "").replace("]", ""))
			i.insert(-3, totseqW)
	
	information.extend(toAddInInformation) # extend information with duplicated entries with different D (from previous A and B)
	if seq == "amplicon" and GENE == "IGH":
		for readRecovered in readsRecoveredNonMatched:
			if readRecovered not in readsAlreadyRecovered:
				readsAlreadyRecovered.append(readRecovered)

	# Round 2: recover reads J-D and D-V (allow assignment to multiple potential rearrangements sharing the same junction)
	for i in information:

		# Get soft clipped start/end J-V :  
		if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
			breakJ = int(i[5])
			breakV = int(i[7]) if i[1] == "Deletion" else int(i[8]) # for IGK Inversion2
		else:
			breakV = int(i[4]) # it is J 
			breakJ = int(i[8]) if i[1] == "Deletion" else int(i[7]) # it is V  # for TRB Inversion1
		
		if "IGKKde" in i[0] or "IGKRSS" in i[0]:
			continue
		
		else:
			ANNOT_TABLE_JV = open(annot_table_JV, "r")
			for j in ANNOT_TABLE_JV:
				w = j.rstrip("\n").split("\t")
				if w[0] in readsAlreadyRecovered: continue

				if w[11].startswith("split"):

					# if information last position J:
					if breakJ == int(w[12].replace("NA", "0")) and w[13] == "NA":							
						split = re.findall(r'[A-Za-z]|[0-9]+', w[5])
						cigar1 = [split[x:x+2] for x in range(0, len(split),2)]
						# MS cigar
						if min([count for count, item in enumerate(cigar1) if "M" in item]) < min([count for count, item in enumerate(cigar1) if "S" in item]):
							if w[20] == "Deletion": subcigar1 = cigar1[max([i for i in range(len(cigar1)) if cigar1[i][1] == "M"]):] # limit cigar1 to start at last M
							else: subcigar1 = cigar1
							mStart = sum([ int(x[0]) if x[1] == "S" else 0 for x in subcigar1 ]) # we add numbers previous to S
							J = w[9][-mStart:]
						# SM cigar
						else:
							if w[20] == "Deletion": subcigar1 = cigar1[:min([i for i in range(len(cigar1)) if cigar1[i][1] == "M"])] # limit cigar1 to first M
							else: subcigar1 = cigar1
							mEnd = sum([ int(x[0]) if x[1] == "S" else 0 for x in subcigar1 ]) # we add numbers previous to S
							J = w[9][:mEnd]
						
						# Vseq: remove deleted nucleotides, check insertion at first bases, keep insertions not at first base:
						if GENE not in ["IGL", "TRA", "TRB", "TRD"]: 
							vSeq = re.sub("\(.*?\)", "",  i[12])
							if vSeq[0] == "[":
								vSeq = vSeq[min([a.start()+1 for a in re.finditer("\]", vSeq)]):]
						
						else: # it is J in IGL/TRA/TRB/TRD
							vSeq = re.sub("\(.*?\)", "",  i[10])
							if vSeq[-1] == "]":
								vSeq = vSeq[:max([a.start() for a in re.finditer("\[", vSeq)]):]
						
						vSeq = vSeq.replace("[", "").replace("]", "")
						dvSeq = i[11]+vSeq # add D prior Vseq

						if i[1] == "Inversion1": J = ''.join(complement[base] for base in reversed(J))

						if dvSeq.startswith(J):
							if GENE not in ["IGL", "TRA", "TRB", "TRD"]: 
								i[6] += 1 # count split J
								i[17] = w[0] if i[17] == "" else i[17]+","+w[0]
							else: 
								i[9] += 1 # count split V
								i[18] = w[0] if i[18] == "" else i[18]+","+w[0]
							# readsAlreadyRecovered.append(w[0]) # commented to allow here read to be considered in multiple potential rearrangements
					
					# if information first position V:
					elif breakV == int(w[12].replace("NA", "0")) and w[13] == "NA":
						split = re.findall(r'[A-Za-z]|[0-9]+', w[5])
						cigar1 = [split[x:x+2] for x in range(0, len(split),2)]
						# MS cigar
						if min([count for count, item in enumerate(cigar1) if "M" in item]) < min([count for count, item in enumerate(cigar1) if "S" in item]):
							if w[20] == "Deletion": subcigar1 = cigar1[max([i for i in range(len(cigar1)) if cigar1[i][1] == "M"]):] # limit cigar1 to start at last M
							else: subcigar1 = cigar1
							mStart = sum([ int(x[0]) if x[1] == "S" else 0 for x in subcigar1 ]) # we add numbers previous to S
							V = w[9][-mStart:]
						# SM cigar
						else:
							if w[20] == "Deletion": subcigar1 = cigar1[:min([i for i in range(len(cigar1)) if cigar1[i][1] == "M"])] # limit cigar1 to first M
							else: subcigar1 = cigar1
							mEnd = sum([ int(x[0]) if x[1] == "S" else 0 for x in subcigar1 ]) # we add numbers previous to S
							V = w[9][:mEnd]
						
						# jSeq: remove deleted nucleotides, check insertion at last bases, keep insertions not at last base:
						if GENE not in ["IGL", "TRA", "TRB", "TRD"]: 
							jSeq = re.sub("\(.*?\)", "",  i[10])
							if jSeq[-1] == "]":
								jSeq = jSeq[:max([a.start() for a in re.finditer("\[", jSeq)]):]
						else:  # it is V in IGL/TRA/TRB/TRD
							jSeq = re.sub("\(.*?\)", "",  i[12])
							if jSeq[0] == "[":
								jSeq = jSeq[min([a.start()+1 for a in re.finditer("\]", jSeq)]):]
								
						jSeq = jSeq.replace("[", "").replace("]", "")
						jdSeq = jSeq+i[11] # add D after Jseq

						if i[1] == "Inversion2": V = ''.join(complement[base] for base in reversed(V))
						
						if jdSeq.endswith(V):
							if GENE not in ["IGL", "TRA", "TRB", "TRD"]: 
								i[9] += 1 # count split V
								i[18] = w[0] if i[18] == "" else i[18]+","+w[0]
							else: 
								i[6] += 1 # count split J
								i[17] = w[0] if i[17] == "" else i[17]+","+w[0]							
							# readsAlreadyRecovered.append(w[0]) # commented to allow here read to be considered in multiple potential rearrangements
			
			ANNOT_TABLE_JV.close()

	# sort and return information
	information.sort(key=lambda p: round(p[2]*2 + p[3] + p[6]*2 + p[9]*2, 1), reverse=True)
	return(information)

def removeLowSupportRearrangements(information, tumorPurity, scoreCutoffFilter, keepInsertSizeOnlyRearrangements, annot_table_JV, seqDepth):
	
	informationClean = [] # list to store information that pass filters
	informationCleanTmp = [] # list to store information that pass filters temporarily
	informationCleanTmpInsert = [] # list to store information of insert size only rearrangements
	pairsOfGenes = [] # list to store pairs of genes that pass filters
	readNamesUsed = [] # list to store readNames used
	
	# 1. keep rearrangements with split-reads
	for i in information:
		spl = int(i[2]) + int(i[6]) + int(i[9])
		seqDepthHighCutoffSpl = 1 if keepInsertSizeOnlyRearrangements == "no" else 0
		if (seqDepth != "high" and spl > 0) or (seqDepth == "high" and spl > seqDepthHighCutoffSpl):
			informationCleanTmp.append(i)
			if i[0] not in pairsOfGenes: pairsOfGenes.append(i[0])
			if i[-3] != "": readNamesUsed.extend(i[-3].split(","))
			if i[-2] != "": readNamesUsed.extend(i[-2].split(","))
			if i[-1] != "": readNamesUsed.extend(i[-1].split(","))
			
	# 2. if keepInsertSizeOnlyRearrangements == "yes", recover split-insertSize reads not used during D-seq and keep (split-)insertSize only rearrangements
	if keepInsertSizeOnlyRearrangements == "yes":
		# get rearrangements not supported by split reads
		for i in information:
			if i[0] not in pairsOfGenes:
				informationCleanTmpInsert.append(i)
				pairsOfGenes.append(i[0])
		# recover insertSize and split-insertSize reads not used in split-read-supported rearrangements
		if len(informationCleanTmpInsert) > 0:
			for i in informationCleanTmpInsert:
				jGene = i[0].split(" - ")[0]
				vGene = i[0].split(" - ")[1]
				jvReadsInsert = []
				ANNOT_TABLE_JV = open(annot_table_JV, "r")
				for jvLine in ANNOT_TABLE_JV: 
					jvList = jvLine.rstrip("\n").split("\t")
					if jvList[0] in readNamesUsed: continue
					if (jGene == jvList[18] and vGene == jvList[19]) or (jGene == jvList[19] and vGene == jvList[18]):
						jvReadsInsert.append(jvList[0])	
						readNamesUsed.append(jvList[0])		
				ANNOT_TABLE_JV.close()
				if len(jvReadsInsert) > 0:
					i[3] = len(jvReadsInsert)
					i[-3] = ",".join(jvReadsInsert)
			# extend to informationCleanTmp
			informationCleanTmp.extend(informationCleanTmpInsert)

	# 3. remove rearrangements if score is 0 or lower than scoreCutoffFilter
	for i in informationCleanTmp:
		spl_ins = round( ( int(i[2])*2 + int(i[3]) + int(i[6])*2 + int(i[9])*2 ) / tumorPurity , 1 )
		if spl_ins == 0: continue
		if spl_ins < scoreCutoffFilter: continue
		informationClean.append(i)
	
	# 4. return
	return(informationClean)

def doCollapseSequences(information, collapseSequencesSimilarity, GENE):	

	# initiate dict and list
	updatedInformation = []
	seqDict = dict()
	JVpairsCounted = []

	# iterate information list
	for p in information:
		
		if p[0] == "IGKKde - IGKRSS": 
			updatedInformation.append(p)
			continue
		
		if p[16] != "NA":
			seq = p[16] # p[16] = nucleotide sequence
			seqIs = "sequence"
		else: # if p[16] == "NA" --> IGKKde - IGKVX(D)-X rearrangement
			seq = p[0].replace("D-", "-")
			seqIs = "geneIDs"
		
		foundInDict = "no"
		if seq in seqDict:
			foundInDict = "yes"
			seqInDictToMatch = seq
		elif seqIs == "sequence":
			for seqInDict in seqDict:
				if seqDict[seqInDict][16] == "NA": 
					continue
				# check if same length, similar sequence
				elif len(seq) == len(seqInDict):
					if CustomSequenceMatcher(seq, seqInDict, "nt") > collapseSequencesSimilarity:
						foundInDict = "yes"
						seqInDictToMatch = seqInDict
						break
				# check if partial sequence matches exactly a complete sequence
				elif seqInDict.startswith(seq) or seqInDict.endswith(seq):
					foundInDict = "yes_butPartial"
					seqInDictToMatch = seqInDict
					break
				# check if complete sequence matches exactly a partial sequence
				elif seq.startswith(seqInDict) or seq.endswith(seqInDict):
					foundInDict = "yes_butPartial_needToChangeSeq"
					seqInDictToMatch = seqInDict
					break
				# check if different length, similar sequence
				else:
					aligner = PairwiseAligner()
					aligner.mode = "global"
					aligner.match_score = 1
					aligner.mismatch_score = -1
					aligner.open_gap_score = -3
					aligner.extend_gap_score = -1
					alignments = aligner.align(seq, seqInDict)
					numAlignments = 0
					try: numAlignments = len(alignments)
					except OverflowError: numAlignments = 1000 
					if numAlignments > 0:
						if str(alignments[0]).startswith("target"):
							align_seq = str(alignments[0][0])
							align_seqInDict = str(alignments[0][1])
						else:
							align_seq = str(alignments[0]).split("\n")[0]
							align_seqInDict = str(alignments[0]).split("\n")[2]
						ratio = sum(1 for a, b in zip(align_seq, align_seqInDict) if a == b) / len(align_seq)
						if ratio > collapseSequencesSimilarity:
							# if so, get sequence with less Ns...
							if sum(1 for nuc in seqInDict if nuc == "N") <= sum(1 for nuc in seq if nuc == "N"):
								foundInDict = "yes_butNotIdenticalSeq"
								seqInDictToMatch = seqInDict
							else:
								foundInDict = "yes_butNotIdenticalSeq_needToChangeSeq"
								seqInDictToMatch = seqInDict
							break
		
		if foundInDict == "no":
			seqDict[seq] = p
		else:
			# update mechanism, number of reads, phasing, sequence (if Ns), and readnames
			## combine mechanism if deletion and inversion2 (IGK)
			if p[1] not in seqDict[seqInDictToMatch][1]:
				seqDict[seqInDictToMatch][1] = seqDict[seqInDictToMatch][1]+";"+p[1]
			## sum J-V N_split
			seqDict[seqInDictToMatch][2] = seqDict[seqInDictToMatch][2]+p[2] 
			## sum N_insertSize if different J-V pair and not already counted
			J = p[0].split(" - ")[0]
			V = p[0].split(" - ")[-1]
			JVpair = J+" - "+V
			if ( seqDict[seqInDictToMatch][0].split(" - ")[0] != J or seqDict[seqInDictToMatch][0].split(" - ")[-1] != V ) and not JVpair in JVpairsCounted:
				seqDict[seqInDictToMatch][3] = seqDict[seqInDictToMatch][3]+p[3]
				JVpairsCounted.append(JVpair)
			## merge readNames
			seqDict[seqInDictToMatch][17] = ",".join([readName for readName in set(seqDict[seqInDictToMatch][17].split(",") + p[17].split(",")) if readName != ""])
			seqDict[seqInDictToMatch][18] = ",".join([readName for readName in set(seqDict[seqInDictToMatch][18].split(",") + p[18].split(",")) if readName != ""])
			seqDict[seqInDictToMatch][19] = ",".join([readName for readName in set(seqDict[seqInDictToMatch][19].split(",") + p[19].split(",")) if readName != ""])
			seqDict[seqInDictToMatch][20] = ",".join([readName for readName in set(seqDict[seqInDictToMatch][20].split(",") + p[20].split(",")) if readName != ""])
			## N_split_rescued_J	
			seqDict[seqInDictToMatch][6] = 0 if seqDict[seqInDictToMatch][18] == "" else len(seqDict[seqInDictToMatch][18].split(","))
			## N_split_rescued_V			
			seqDict[seqInDictToMatch][9] = 0 if seqDict[seqInDictToMatch][19] == "" else len(seqDict[seqInDictToMatch][19].split(","))
			## Num of rescued reads during phasing of mutations
			seqDict[seqInDictToMatch][10] = 0 if seqDict[seqInDictToMatch][20] == "" else len(seqDict[seqInDictToMatch][20].split(","))
			## assign max percentage of phasing
			if foundInDict == "yes" and p[15] != "NA":
				phaseDict = 0 if seqDict[seqInDictToMatch][15].startswith("NA") else 100 if seqDict[seqInDictToMatch][15].startswith("0/0 -") else float(seqDict[seqInDictToMatch][15].split(" - ")[0].split("/")[0]) / float(seqDict[seqInDictToMatch][15].split(" - ")[0].split("/")[1])
				phaseNew = 0 if p[15].startswith("NA") else 100 if p[15].startswith("0/0 -") else float(p[15].split(" - ")[0].split("/")[0]) / float(p[15].split(" - ")[0].split("/")[1])
				if phaseDict < phaseNew: seqDict[seqInDictToMatch][15] = p[15]
			## adjust Ns in sequence
			if foundInDict == "yes" and seqDict[seqInDictToMatch][16] != "NA" and "N" in seqDict[seqInDictToMatch][16]:
				seqInDict = seqDict[seqInDictToMatch][16]
				adjustedSeq = ''.join(seq[i] if seqInDict[i] == 'N' else seqInDict[i] for i in range(len(seqInDict)))
				seqDict[seqInDictToMatch][16] = adjustedSeq # update JDV tumor seq
				if GENE in ["IGL", "TRA", "TRB", "TRD"]:
					seqDict[seqInDictToMatch][11] = adjustedSeq[-len(seqDict[seqInDictToMatch][11]):] # update J
					seqDict[seqInDictToMatch][13] = adjustedSeq[:len(seqDict[seqInDictToMatch][13])] # update V tumor
				else:
					seqDict[seqInDictToMatch][11] = adjustedSeq[:len(seqDict[seqInDictToMatch][11])] # update J
					seqDict[seqInDictToMatch][13] = adjustedSeq[-len(seqDict[seqInDictToMatch][13]):] # update V tumor

			## if partial sequence in dict, update with complete sequence information or if very similar sequence coupled with need to change it
			if foundInDict == "yes_butPartial_needToChangeSeq" or foundInDict == "yes_butNotIdenticalSeq_needToChangeSeq":
				seqDict[seqInDictToMatch][0] = p[0]
				seqDict[seqInDictToMatch][4] = p[4]
				seqDict[seqInDictToMatch][5] = p[5]
				seqDict[seqInDictToMatch][7] = p[7]
				seqDict[seqInDictToMatch][8] = p[8]
				seqDict[seqInDictToMatch][11] = p[11]
				seqDict[seqInDictToMatch][12] = p[12]
				seqDict[seqInDictToMatch][13] = p[13]
				seqDict[seqInDictToMatch][14] = p[14]
				seqDict[seqInDictToMatch][15] = p[15]
				seqDict[seqInDictToMatch][16] = p[16]
				seqDict[seq] = seqDict.pop(seqInDictToMatch)
	
	# return updatedInformation	
	updatedInformation.extend(list(seqDict.values()))
	updatedInformation.sort(key=lambda p: round(p[2]*2 + p[3] + p[6]*2 + p[9]*2, 1), reverse=True)
	return(updatedInformation)

def calculateHomology(hom, homN):
	ct = 0 # sum total nucleotides
	ce = 0 # sum equal nucleotides
	j = 0
	
	while j < len(hom):
		if hom[j] != "N" and homN[j] != "N" and homN[j] != "-":	
			ct += 1
			if homN[j] == hom[j] or hom[j] == "-":
				ce += 1
		j += 1
	if ct > len(hom)/2:
		pctHomology = str(round((ce/ct)*100, 2)) # homology
		numHomology = str(ce)+"/"+str(ct) # number nucleotides homology
	else:
		pctHomology = "NA"
		numHomology = "NA"
	
	return(pctHomology, numHomology)

def productivityAndCDR3(vdj, nCDR3, GENE):
	
	# set output
	productiu = "NA"
	aaCDR3 = "NA"
	
	vdjTrip = re.findall('.{3}', vdj)

	# stop codons in seq:
	if "TAA" in vdjTrip or "TAG" in vdjTrip or "TGA" in vdjTrip:
		productiu1 = "Unproductive (stop codons)"
	else:
		productiu1 = "Phe118 not identified (check at IMGT/IgBlast)"

	if 'TTT' in nCDR3 or 'TTC' in nCDR3 or 'TGG' in nCDR3:
		
		for wf in reversed(list(re.finditer('TTT|TTC|TGG', nCDR3, overlapped=True))): 
			
			j118 = "NA"
			
			productiu = productiu1	# re-set productiu
			
			po = wf.start()
			
			tripMotif = re.findall('.{3}', nCDR3[po:po+12])
			
			if len(tripMotif) == 4:
				if "N" in tripMotif[1]: continue
		
				if "N" in tripMotif[3]: 
					if tripletsToAA[tripMotif[1]] == "G":
						j118 = po+3
		
				elif GENE.startswith("IG") and (tripletsToAA[tripMotif[1]] == "G" or tripletsToAA[tripMotif[3]] == "G"):
					j118 = po+3
				
				elif GENE.startswith("TR") and tripletsToAA[tripMotif[1]] == "G" and tripletsToAA[tripMotif[3]] == "G":
					j118 = po+3
			
			elif len(tripMotif) > 1 and "N" not in tripMotif[1]:
				if tripletsToAA[tripMotif[1]] == "G":
					j118 = po+3
			
			if j118 != "NA": 
				nCDR3_a = nCDR3[:j118]
		
				# out-of-frame
				if not (len(nCDR3_a)/3).is_integer(): 				
					if productiu == "Unproductive (stop codons)": productiu = "Unproductive (stop codons, out-of-frame junction)"
					else: productiu = "Unproductive (out-of-frame junction)"
			
					if ((len(nCDR3_a)+1)/3).is_integer(): add = "."
					else: add = ".."
					nCDR3_a = nCDR3_a[:-9]+add+nCDR3_a[-9:]
					tripCDR3 = re.findall('.{3}', nCDR3_a)  
					aaCDR3 = "".join([tripletsToAA[aa] if "." not in aa and "N" not in aa else "?" if "N" in aa else "#" for aa in tripCDR3])

				# in frame
				else: 
					if productiu == "Unproductive (stop codons)":
						productiu = "Unproductive (stop codons, in-frame junction)"
					else: 
						productiu = "Productive (no stop codon and in-frame junction)"
					
					tripCDR3 = re.findall('.{3}', nCDR3_a)
					aaCDR3 = "".join([tripletsToAA[aa] if "N" not in aa else "?" for aa in tripCDR3])
			
					if productiu == "Productive (no stop codon and in-frame junction)": break

			elif j118 == "NA":
				nCDR3_a = nCDR3
				tripCDR3 = re.findall('.{3}', nCDR3_a)
				aaCDR3 = "".join([tripletsToAA[aa] if "N" not in aa else "?" for aa in tripCDR3])
				if productiu == "Unproductive (stop codons)": productiu = "Unproductive (stop codons, Phe118 not identified)"
	
	else:
		productiu = productiu1
		if productiu == "Unproductive (stop codons)": productiu = "Unproductive (stop codons, Phe118 not identified)"
		nCDR3_a = nCDR3
		tripCDR3 = re.findall('.{3}', nCDR3_a)
		aaCDR3 = "".join([tripletsToAA[aa] if "N" not in aa else "?" for aa in tripCDR3])
		
	return(productiu, aaCDR3)

def checkHomologyAndFunctionality(information, GENE):
	
	for p in information:
		lstProductiu = ["NA"]*5
		productiu, aaCDR3, pctHomology, numHomology, hom = lstProductiu
		
		if len(p[0].split(" - ")) == 2 and (p[0].split(" - ")[0][3] == "D" or p[0].split(" - ")[1][3] == "D"): # partialRearrangement (J-D only)
			productiu = "Partial rearrangement"

		elif p[4] == "NA":
			productiu = "Rearrangement without junction coverage"

		elif "Kde" not in p[0] and "RSS" not in p[0]:
			
			productiu = "No junction found"
			
			# remove insertions from V to check productivity (but keep deletions):
			tumSeq = re.sub("\[.*?\]", "", p[13]).replace("(", "").replace(")", "")
			normSeq = re.sub("\[.*?\]", "", p[14]).replace("(", "").replace(")", "")

			# complement if needed
			if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
				t = ''.join(complement[base] for base in reversed(tumSeq))+''.join(complement[base] for base in  reversed(p[12]))  # reverse complement of V + D
				n = ''.join(complement[base] for base in reversed(normSeq))+''.join(complement[base] for base in  reversed(p[12])) # reverse complement of V normal + D tumor
			else:
				t = tumSeq+p[12] # V tumor + D tumor
				n = normSeq+p[12] # V normal + D tumor
			
			# Use N instead of ambiguous bases
			t = ''.join(changeAmbiguousBases[base] for base in t)
			n = ''.join(changeAmbiguousBases[base] for base in n)
			
			# Check: Cys23 (Trp41) Cys104 to define FR1-FR3 for homology		
			for cys in reversed(list(re.finditer('TGT|TGC', t, overlapped=True))): 
				
				countAAfound = 1 # count cys23, trp41, cys104, phe118 found
				
				if productiu.startswith(("Productive", "Phe118 not identified")): break # if already found
			
				posCys104 = "NA"
				
				po = int(cys.start()) # Cys23 
				substractLeft = 21 if (GENE == "IGH" or (GENE.startswith("IG") and t[po-21*3:].startswith("CAG"))) else 22 # aa10 missing in some IGH/K/L; check if missing by considering IGH/K/L often starts with "CAG"
				if po-substractLeft*3 < 0: continue # make sure potential Cys23 is not just at the the begining of the sequence
				
				nextpo = po + 14*3 if GENE == "IGH" else po + 9*3 # Trp41 (9 aa x 3bp)
				
				if nextpo > len(t): continue
				
				# Trp41
				if "TGG" in re.findall('.{3}', t[nextpo:nextpo+30]): 
					
					countAAfound += 1 
					
					nextpo += max([i for i in range(len(re.findall('.{3}', t[nextpo:nextpo+30]))) if re.findall('.{3}', t[nextpo:nextpo+30])[i] == 'TGG']) * 3 
					nextnextpo = nextpo + 58*3 if GENE in ["IGH", "TRG"] else nextpo + 48*3 # for cys 104
					
					if nextnextpo > len(t): continue
					
					# Cys104
					if 'TGT' in re.findall('.{3}', t[nextnextpo:nextnextpo+36]) or 'TGC' in re.findall('.{3}', t[nextnextpo:nextnextpo+36]):
						countAAfound += 1 
						posCys104 = [ nextnextpo + cysUnref * 3 + 3 for cysUnref in [i for i in range(len(re.findall('.{3}', t[nextnextpo:nextnextpo+36]))) if re.findall('.{3}', t[nextnextpo:nextnextpo+36])[i] == 'TGT' or re.findall('.{3}', t[nextnextpo:nextnextpo+36])[i] == 'TGC'] ] #+3 for index 0
				
				# No Trp41, check Cys104 directly	
				else: 
					nextnextpo = nextpo + 58*3 if GENE in ["IGH", "TRG"] else nextpo + 48*3 # for cys 104
					if 'TGT' in re.findall('.{3}', t[nextnextpo:nextnextpo+36]) or 'TGC' in re.findall('.{3}', t[nextnextpo:nextnextpo+36]):
						countAAfound += 1 
						
						posCys104 = [ nextnextpo + cysUnref * 3 + 3 for cysUnref in [i for i in range(len(re.findall('.{3}', t[nextnextpo:nextnextpo+36]))) if re.findall('.{3}', t[nextnextpo:nextnextpo+36])[i] == 'TGT' or re.findall('.{3}', t[nextnextpo:nextnextpo+36])[i] == 'TGC'] ] #+3 for index 0
						
				# Check homolgy and functionality if Cys104 found
				if posCys104 != "NA":
					
					posCys104_list = posCys104    
					
					for posCys104 in posCys104_list:		
						
						vSeq = t[po-substractLeft*3:] # to be used for cdr3: store FR1 - endOfV
						newPosCys104 = posCys104-3-(po-substractLeft*3) # to be used for cdr3: get position Cys104 from FR1, not from beggining of V gene
						hom = t[po-substractLeft*3:posCys104]
						homN = n[po-substractLeft*3:posCys104]							
						
						# calculateHomology:
						pctHomology, numHomology = calculateHomology(hom, homN)
						
						# if complete (D found for IGH/TRB/TRD or not IGH/TRB/TRD and FR1-FR3 found previously), check if productive and get CDR3:
						if ((GENE in ["IGH", "TRB", "TRD"] and p[12] != "") or GENE not in ["IGH", "TRB", "TRD"]) and hom != "NA":
							if GENE not in ["IGL", "TRA", "TRB", "TRD"]:
								j = ''.join(complement[base] for base in reversed(re.sub("\(.*?\)", "", p[11]).replace("[", "").replace("]", ""))) # remove deletions only from J to check productivity
							else:
								j = re.sub("\(.*?\)", "", p[11].replace("[", "").replace("]", "")) # remove deletions only from J to check productivity
							
							# VDJ sequence and CDR3
							vdj = vSeq+j
							nCDR3 = vdj[newPosCys104:].rstrip("N") 
							
							# productivityAndCDR3
							productiu, aaCDR3 = productivityAndCDR3(vdj, nCDR3, GENE)
							
							if lstProductiu == ["NA"]*5 or productiu.startswith("Productive"):
								lstProductiu = [ productiu, aaCDR3, pctHomology, numHomology, hom ]

							if productiu.startswith("Productive"): break
			
			productiu, aaCDR3, pctHomology, numHomology, hom = lstProductiu
			
			
			# if productive rearrangement not found starting with Cys23, start at Trp41 (IG ONLY):
			lstProductiu_2 = ["NA"]*5	
			productiu_2, aaCDR3_2, pctHomology_2, numHomology_2, hom_2 = lstProductiu_2
			
			if not productiu.startswith("Productive") and GENE.startswith("IG"): 
				# Trp41
				for trp in reversed(list(re.finditer('TGG', t, overlapped=True))): 
					
					if productiu_2.startswith(("Productive", "Phe118 not identified")): break # if already found
					
					posCys104 = "NA"
					
					nextpo = int(trp.start()) # Trp41
					if GENE == "IGH":
						substractLeft = 35 # not 40 because aa10 is often missing and 4 aa often missing at CDR1
					elif t[nextpo-33*3:].startswith("CAG"):
						substractLeft = 33 # not 40 because aa10 is often missing and 6 aa often missing at CDR1 in IGK/IGL
					else:
						substractLeft = 34 # 6 aa often missing at CDR1 in IGK/IGL

					if nextpo-substractLeft*3 < 0: continue # make sure potential Trp41 is not at the the begining of the sequence
					
					nextnextpo = nextpo + 58*3 if GENE == "IGH" else nextpo + 48*3 # for cys 104
					
					# Cys104
					if 'TGT' in re.findall('.{3}', t[nextnextpo:nextnextpo+36]) or 'TGC' in re.findall('.{3}', t[nextnextpo:nextnextpo+36]): 
						
						posCys104_list = [ nextnextpo + cysUnref * 3 + 3 for cysUnref in [i for i in range(len(re.findall('.{3}', t[nextnextpo:nextnextpo+36]))) if re.findall('.{3}', t[nextnextpo:nextnextpo+36])[i] == 'TGT' or re.findall('.{3}', t[nextnextpo:nextnextpo+36])[i] == 'TGC'] ] #+3 for index 0
						
						# check homolgy and functionality	
						for posCys104 in posCys104_list:
							
							# check homolgy and functionality								
							vSeq = t[nextpo-substractLeft*3:] # to be used for cdr3: store FR1 - endOfV
							newPosCys104 = posCys104-3-(nextpo-substractLeft*3) # to be used for cdr3: get position Cys104 from FR1, not from beggining of V gene
							hom_2 = t[nextpo-substractLeft*3:posCys104]
							homN = n[nextpo-substractLeft*3:posCys104]
							
							# calculateHomology
							pctHomology_2, numHomology_2 = calculateHomology(hom_2, homN)
							
							# if complete (D found for IGH or not IGH and FR1-FR3 found previously), check if productive and get CDR3 sequence:
							if ((GENE == "IGH" and p[12] != "") or GENE != "IGH") and hom_2 != "NA":
								if GENE != "IGL":
									j = ''.join(complement[base] for base in reversed(re.sub("\(.*?\)", "", p[11]).replace("[", "").replace("]", ""))) # remove deletions only from J to check productivity
								else:
									j = re.sub("\(.*?\)", "", p[11]).replace("[", "").replace("]", "") # remove deletions only from J to check productivity
					
								# VDJ sequence and CDR3
								vdj = vSeq+j
								nCDR3 = vdj[newPosCys104:].rstrip("N")
								
								# productivityAndCDR3
								productiu_2, aaCDR3_2 = productivityAndCDR3(vdj, nCDR3, GENE)
								
								if lstProductiu_2 == ["NA"]*5 or productiu_2.startswith("Productive"):
									lstProductiu_2 = [ productiu_2, aaCDR3_2, pctHomology_2, numHomology_2, hom_2 ]
								
								if productiu_2.startswith(("Productive", "Phe118 not identified")): break	
			
			productiu_2, aaCDR3_2, pctHomology_2, numHomology_2, hom_2 = lstProductiu_2
			
			# merge starting at Cys23 and Trp41:
			if productiu_2.startswith("Productive") or ( productiu.startswith("Unproductive") and productiu_2.startswith("Phe118 not identified") ) or  ( productiu == "Unproductive (stop codons, out-of-frame junction)" and productiu_2 == "Unproductive (stop codons, in-frame junction)" ):
				productiu = productiu_2
				aaCDR3 = aaCDR3_2
				pctHomology = pctHomology_2
				numHomology = numHomology_2
				hom = hom_2
			
			elif productiu.startswith("Unproductive") and productiu_2.startswith("Unproductive"):
				if len(aaCDR3_2) < len(aaCDR3):
					productiu = productiu_2
					aaCDR3 = aaCDR3_2
					pctHomology = pctHomology_2
					numHomology = numHomology_2
					hom = hom_2
			
			
			# if no homolgy calculated before, determine homology entire V gene:			
			if pctHomology == "NA":
				t = tumSeq # V tumor
				n = normSeq # V normal				
				pctHomology, numHomology = calculateHomology(t, n)
			
			if "[" in p[13] or "(" in p[13]:
				productiu = productiu+" [indel(s) in seq]"
				
		p.append(pctHomology)
		p.append(numHomology)
		p.append(productiu)
		p.append(aaCDR3)
		p.append(p[0]) # final gene annotation (for uniformity when running IgBlast)
		
	return(information)

def removePrimerFromSeq(seqToIgblast, GENE, primerFasta):
	
	strand = "pos" if GENE in ["IGL", "TRA", "TRB", "TRD"] else "neg"
	
	# remove N from V gene
	if strand == "pos":
		seqToIgblast_woN = seqToIgblast.lstrip("N")
	else: 
		seqToIgblast_woN = seqToIgblast.rstrip("N")

	# initialize variables	
	nextLine = False
	matchFound = False
	maxScore = "NA"
	lenMaxScore = "NA"

	# iterate primers and align
	FASTA = open(primerFasta, "r")
	for fLine in FASTA:
		if fLine.startswith(">"+GENE):
			nextLine = True
		elif fLine.startswith(">"):
			nextLine = False
		elif nextLine:
			pSeq = fLine.rstrip("\n")
			if strand == "pos":
				seqToMap = seqToIgblast_woN[:len(pSeq)]
			else: 
				pSeq = ''.join(complement[base] for base in reversed(pSeq))
				seqToMap = seqToIgblast_woN[-len(pSeq):]

			# align
			aligner = PairwiseAligner()
			# force stitching without internal gap
			aligner.mode = "global"
			aligner.match_score = 2
			aligner.mismatch_score = -1
			aligner.open_gap_score = -1000
			aligner.extend_gap_score = -1000
			aligner.left_open_gap_score = 0
			aligner.left_extend_gap_score = 0
			aligner.right_open_gap_score = 0
			aligner.right_extend_gap_score = 0
			alignments = aligner.align(seqToMap, pSeq)
			numAlignments = 0
			try: numAlignments = len(alignments)
			except OverflowError: numAlignments = 1000 
			if numAlignments > 0:
				alignment_score = float(alignments[0].score)
				if str(alignments[0]).startswith("target"):
					align_seqToMap = str(alignments[0][0])
					align_pSeq = str(alignments[0][1])
				else:
					align_seqToMap = str(alignments[0]).split("\n")[0]
					align_pSeq = str(alignments[0]).split("\n")[2]
				
				# remove gaps from primer sequence alignment
				if align_pSeq.startswith("-"): start_count = len(align_pSeq) - len(align_pSeq.lstrip('-'))
				else: start_count = 0
				if align_pSeq.endswith("-"): 
					end_count = len(align_pSeq) - len(align_pSeq.rstrip('-'))
					corrected_seqToMap = align_seqToMap[start_count:-end_count]
				else: 
					corrected_seqToMap = align_seqToMap[start_count:]
				
				# remove gaps from aligned sequenced
				corrected_seqToMap = corrected_seqToMap.rstrip("-").lstrip("-").replace("-", "")

				if maxScore == "NA" or alignment_score > maxScore:
					maxScore = alignment_score
					lenMaxScore = len(corrected_seqToMap)
					matchFound = True
			
			# restart nextLine
			nextLine = False
	FASTA.close()
	
	# return clean sequence
	if matchFound:
		if strand == "pos": seqToIgblast_woN = seqToIgblast_woN[lenMaxScore:]
		else: seqToIgblast_woN = seqToIgblast_woN[:-lenMaxScore]
	
	return(seqToIgblast_woN)

def annotateAGS(igblastList):

	locationsList = []
	agsList = []

	# get aa sequence from igblast
	fr1 = igblastList[37]
	cdr1 = igblastList[39]
	fr2 = igblastList[41]
	cdr2 = igblastList[43]
	fr3 = igblastList[45]
	cdr3 = igblastList[49]
	fr4 = igblastList[47]
	aaSeq = fr1+cdr1+fr2+cdr2+fr3+cdr3+fr4
	aaGermSeq = igblastList[17]

	# scan for N-glycosylation motif (N-X-T/S, where X is any amino acid except proline)
	motif = r'N[^P][TS]'
	matches = re.finditer(motif, aaSeq, overlapped=True)

	# iterate, chek if same motif in germline, annotate location/gene
	for match in matches:
		
		# is motif in seq is not the same in germline
		if aaSeq[match.start():match.end()] != aaGermSeq[match.start():match.end()]:
			
			mGroup = match.group()
			mStart = match.start() + 1
			mEnd = match.end()

			if mEnd <= len(fr1):
				mLocation = "FR"
				mGene = "FR1"
			elif mStart <= len(fr1+cdr1):
				mLocation = "CDR"
				mGene = "CDR1"
			elif mEnd <= len(fr1+cdr1+fr2):
				mLocation = "FR"
				mGene = "FR2"
			elif mStart <= len(fr1+cdr1+fr2+cdr2):
				mLocation = "CDR"
				mGene = "CDR2"
			elif mEnd <= len(fr1+cdr1+fr2+cdr2+fr3):
				mLocation = "FR"
				mGene = "FR3"
			elif mStart <= len(fr1+cdr1+fr2+cdr2+fr3+cdr3):
				mLocation = "CDR"
				mGene = "CDR3"
			else:
				mLocation = "FR"
				mGene = "FR4"
			
			locationsList.append(mLocation)
			agsList.append(mGroup+":"+mGene)
	
	# prepare output and return
	if locationsList == []:
		agsToReturn = "[No-AGS]"
	else:
		agsType = "CDR-AGS" if "CDR" in locationsList else "FR-AGS"
		agsToReturn = "["+agsType+" ("+";".join(agsList)+")]"
	return agsToReturn

def igBlastAnnotate(information, GENE, wkDir, inputsFolder, germline_db_J, germline_db_D, germline_db_V, seq, primer, primerStringency, primerFasta, agsAnnotation):

	if len(information) > 0:
		platformSystem = "macos" if platform.system() == "Darwin" else "linux"
		ig_seqtype = "Ig" if GENE in ["IGH", "IGK", "IGL"] else "TCR"
		num_alignments_D = "3" if GENE == "TRD" else "1" if GENE in ["IGH", "TRB"] else "0" # based on IMGT

		# Prepare input for IgBlast
		igblast_in = wkDir+"/tmp/igblast_in.txt"
		IGBLAST = open(igblast_in, "w")
		countSeq = 0
		for i in information:
			IGBLAST.write(">seq_"+str(countSeq)+"\n")
			seqToIgblast = "NA" if len(i[0].split(" - ")) == 2 and (i[0].split(" - ")[0][3] == "D" or i[0].split(" - ")[1][3] == "D") else i[16]
			# if seq is amplicon and sequences of primers are specified, remove primer annealing to V gene from sequence
			if seq == "amplicon" and primer == "Vseq" and seqToIgblast != "NA": 
				seqToIgblast = removePrimerFromSeq(seqToIgblast, GENE, primerFasta) # remove primer
				i[16] = seqToIgblast # update sequence in information
			IGBLAST.write(seqToIgblast+"\n")
			countSeq += 1
		IGBLAST.close()

		# Run IgBlast
		igblastn = inputsFolder+"/igblast/"+platformSystem+"/ncbi-igblast-1.22.0/bin/igblastn"
		aux_data = inputsFolder+"/igblast/"+platformSystem+"/ncbi-igblast-1.22.0/optional_file/human_gl.aux"

		igblast_out = wkDir+"/tmp/igblast_out.txt"

		os.environ['IGDATA'] = inputsFolder+"/igblast/"+platformSystem+"/ncbi-igblast-1.22.0"

		comm = igblastn+" -germline_db_V "+germline_db_V+" -germline_db_J "+germline_db_J+" -germline_db_D "+germline_db_D+" -organism human -ig_seqtype "+ig_seqtype+" -auxiliary_data "+aux_data+" -show_translation -outfmt 19 -extend_align5end -extend_align3end -num_alignments_D "+num_alignments_D+" -num_threads 1 -query "+igblast_in+" > "+igblast_out
		e = subprocess.call(comm, shell=True, env=dict(os.environ))
		if e != 0:
			print("IgCaller: error message (1/4)... IgBlast not working properly on your system. Please, see error message above.")
			print("IgCaller: error message (2/4)... Please, solve this issue and run the following column to test IgBlast:")
			print("IgCaller: error message (3/4)... "+igblastn+" --help")
			print("IgCaller: error message (4/4)... Once this error is solved, please, re-run IgCaller.")
			sys.exit(1)
		
		# If seq is amplicon and fasta sequences of primers are not provided, re-run IgBlast using only 'amplified' sequence starting after the primer location
		if seq == "amplicon" and primer != "Vseq":
			IGBLAST = open(igblast_in, "w")
			IGBLAST_INT = open(igblast_out, "r")
			for igblastLine in IGBLAST_INT:
				if igblastLine.startswith("sequence_id"): continue
				igblastList = igblastLine.rstrip("\n").split("\t")
				igblastSeqName = igblastList[0]
				igblastSeq = igblastList[1]
				if igblastSeq != "NA":
					if primerStringency == "strict":
						if primer == "leader" and igblastList[78] != "": posStart =  78 # fr1 start
						elif primer == "fr1" and igblastList[80] != "": posStart = 80 # cdr1 start
						elif primer == "cdr1" and igblastList[82] != "": posStart = 82 # fr2 start
						elif primer == "fr2" and igblastList[84] != "": posStart = 84 # cdr2 start
						elif primer == "cdr2" and igblastList[86] != "": posStart = 86 # fr3 start
						else: posStart = "NA" # NA when primer is at fr3 or region to start is empty
					else: # permissive
						if primer == "leader" and igblastList[78] != "": posStart =  78 # fr1 start
						elif (primer == "fr1" or (primer in ["leader"] and igblastList[78] == "")) and igblastList[80] != "": posStart = 80 # cdr1 start
						elif (primer == "cdr1" or (primer in ["leader", "fr1"] and igblastList[80] == "")) and igblastList[82] != "": posStart = 82 # fr2 start
						elif (primer == "fr2" or (primer in ["leader", "fr1", "cdr1"] and igblastList[82] == "")) and igblastList[84] != "": posStart = 84 # cdr2 start
						elif (primer == "cdr2" or (primer in ["leader", "fr1", "cdr1", "fr2"] and igblastList[84] == "")) and igblastList[86] != "": posStart = 86 # fr3 start
						else: posStart = "NA" # NA when primer is at fr3 or previous regions are empty
					if posStart != "NA":
						igblastSeq = igblastSeq[(int(igblastList[posStart])-1):]
						information[int(igblastList[0].split("_")[1])][16] = igblastSeq # update sequence in information
				IGBLAST.write(">"+igblastSeqName+"\n")
				IGBLAST.write(igblastSeq+"\n")
			IGBLAST_INT.close()
			IGBLAST.close()
			
			comm = igblastn+" -germline_db_V "+germline_db_V+" -germline_db_J "+germline_db_J+" -germline_db_D "+germline_db_D+" -organism human -ig_seqtype "+ig_seqtype+" -auxiliary_data "+aux_data+" -show_translation -outfmt 19 -extend_align5end -extend_align3end -num_alignments_D "+num_alignments_D+" -num_threads 1 -query "+igblast_in+" > "+igblast_out
			subprocess.call(comm, shell=True, env=dict(os.environ))

		# Merge IgBlast output with "information"
		IGBLAST = open(igblast_out, "r")
		for igblastLine in IGBLAST:
			if igblastLine.startswith("sequence_id"): continue
			
			igblastList = igblastLine.rstrip("\n").split("\t")

			if igblastList[1] != "NA":
				# JDV
				igBlast_JDV = igblastList[13]+" - "+((igblastList[12]+" - ") if igblastList[12] != "" else "")+igblastList[11]

				# CDR3
				aaCDR3 = igblastList[52]
				
				# Productive/unproductive
				if aaCDR3 == "":
					aaCDR3 = "NA"
					if igblastList[7] == "F" and igblastList[4] == "T":
						productiu = "Unproductive (stop codons, no junction found)"
					elif igblastList[7] == "F":
						productiu = "Unproductive (no junction found)"
					elif igblastList[7] == "T" and igblastList[5] == "T" and igblastList[4] == "F" and igblastList[9] == "T": # 7 = productive, 5 = vj_in_frame, 4 = stop_codon, 9 = complete
						productiu = "Complete rearrangement (no junction found)"
					else:
						productiu = "No junction found [confirm using IMGT/V-QUEST]"
				elif igblastList[7] == "T" and igblastList[5] == "T" and igblastList[4] == "F": # 7 = productive, 5 = vj_in_frame, 4 = stop_codon
					productiu = "Productive (no stop codon and in-frame junction)"
				elif igblastList[7] == "F" and igblastList[5] == "F" and igblastList[4] == "F": # 7 = productive, 5 = vj_in_frame, 4 = stop_codon
					productiu = "Unproductive (out-of-frame junction)"
				elif igblastList[7] == "F" and igblastList[5] == "F" and igblastList[4] == "T": # 7 = productive, 5 = vj_in_frame, 4 = stop_codon
					productiu = "Unproductive (stop codons, out-of-frame junction)"
				elif igblastList[7] == "F" and igblastList[5] == "T" and igblastList[4] == "T": # 7 = productive, 5 = vj_in_frame, 4 = stop_codon
					productiu = "Unproductive (stop codons)"
				elif igblastList[7] == "T":
					productiu = "Productive [confirm using IMGT/V-QUEST]"
				elif igblastList[7] == "F":
					productiu = "Unproductive [confirm using IMGT/V-QUEST]"
				else:
					productiu = "Check at IMGT/V-QUEST"
				
				# Homology
				fr1fr3_t = "NA"
				fr1fr3_gl = "NA"
				if igblastList[87] != "" and igblastList[78] != "": # FR1-CDR1-FR2-CDR2-FR3
					fr1fr3_length = int(igblastList[87])-int(igblastList[78])+1
					fr1fr3_t = igblastList[14][:fr1fr3_length]
					fr1fr3_gl = igblastList[15][:fr1fr3_length]
				else:
					fr1fr3_t = igblastList[24] # v_sequence_alignment from igblast
					fr1fr3_gl = igblastList[26] # v_germline_alignment from igblast

				if fr1fr3_t != "NA":
					pctHomology, numHomology = calculateHomology(fr1fr3_t, fr1fr3_gl)
				else:
					pctHomology = "NA"
					numHomology = "NA"
					if seq == "amplicon" and primerStringency == "strict" and primer != "fr3":
						productiu = "NA"
						aaCDR3 = "NA"
						igBlast_JDV = "NA"

				# Add "potentially" if indel in sequence:
				if productiu == "Productive (no stop codon and in-frame junction)" and ("-" in fr1fr3_t or "-" in fr1fr3_gl):
					productiu = "Potentially productive (no stop codon and in-frame junction) [indel detected]"
				
				# Annotate AGS:
				if agsAnnotation == "yes" and GENE == "IGH" and (productiu.startswith("Productive") or productiu.startswith("Potentially productive")):
					ags = annotateAGS(igblastList)
					igBlast_JDV = igBlast_JDV+" "+ags

				# Add to information
				p = information[int(igblastList[0].split("_")[1])]
				if p[2] == 0 and p[6] == 0 and p[9] == 0 and p[12] == "": productiu = "Rearrangement without junction coverage"
				p.append(pctHomology)
				p.append(numHomology)
				p.append(productiu)
				p.append(aaCDR3)
				p.append(igBlast_JDV)
			
			else: # sequence == "NA" (i.e. IGKKde - XX or Partial rearrangement or only insertSize reads)
				p = information[int(igblastList[0].split("_")[1])]
				p.append("NA")
				p.append("NA")
				if len(p[0].split(" - ")) == 2 and (p[0].split(" - ")[0][3] == "D" or p[0].split(" - ")[1][3] == "D"): txtToAppend = "Partial rearrangement"
				else: txtToAppend = "NA"
				p.append(txtToAppend)
				p.append("NA")
				p.append(p[0])
	
	return(information)
	
def addMapQualAndScore(information, miniBamT, tumorPurity, pathToSamtools, threadsForSamtools):
	
	for i in information:

		spl_ins = round( ( int(i[2])*2 + int(i[3]) + int(i[6])*2 + int(i[9])*2 ) / tumorPurity , 1 )
		quals = []

		listOfReadsInRearrangement = [readName for readName in i[17].split(",")+i[18].split(",")+i[19].split(",") if readName != ""]
		
		readNamesFileTxt = miniBamT.replace("miniBam.bam", "_readNamesForMQ.txt")
		readNamesFile = open(readNamesFileTxt, "w")
		readNamesFile.write("\n".join(listOfReadsInRearrangement))
		readNamesFile.close()
		subprocess.call(pathToSamtools+"samtools view -@ "+threadsForSamtools+" -N "+readNamesFileTxt+" "+miniBamT+" > "+miniBamT.replace("miniBam.bam", "miniSam_readNamesForMQ.sam"), shell=True)

		reads = open(miniBamT.replace("miniBam.bam", "miniSam_readNamesForMQ.sam"), "r")
		for read in reads:
			readList = read.rstrip("\n").split("\t")
			quals.append(int(readList[4]))
		reads.close()
		
		if quals == []: MQ = "NA"
		else: MQ = str(round(mean(quals),1))+" ("+str(min(quals))+"-"+str(max(quals))+")"
		
		i.append(spl_ins)
		i.append(MQ)

	return(information)

def adjustReadNamesAndPartialRearrangements(information, GENE):
	
	for i in information:
		
		listOfReadsInRearrangement = [readName for readName in i[17].split(",")+i[18].split(",")+i[19].split(",")+i[20].split(",") if readName != ""]
		del i[17:21]
		i.append(",".join(listOfReadsInRearrangement))
		
		if i[19] == "Partial rearrangement":

			partialJDRearrangement = "yes" if i[0].split(" - ")[1][3] == "D" else "no"
			
			if partialJDRearrangement == "yes": # J-D rearrangement
				i[12] = i[12]+i[13] if GENE == "IGH" else i[13]+i[12] # TRB and TRD are in positive strand
				i[13] = "NA"
				i[14] = "NA"
				i[15] = "NA"
			else: # D-V rearrangement
				i[12] = i[11]+i[12] if GENE == "IGH" else i[12]+i[11] # TRB and TRD are in positive strand
				i[11] = "NA"

	return(information)

def annotateCLLsubsets(information, subsetsAnnotation):
	for i in information:
		if (i[19].startswith("Productive") or i[19].startswith("Potentially productive")) and i[17] != "NA":
			
			vGenes = [ vGene.split("*")[0] for vGene in i[21].split(" - ")[-1].split(",") ]
			identity = float(i[17])
			cdr3aa = i[20][1:-1]
			cdr3len = len(cdr3aa)
			
			subset = ""
			if subsetsAnnotation == "imgt": # based on IMGT defintions excluding IGHJ genes (https://www.imgt.org/IMGT_vquest/user_guide#afunccll)
				if "IGHV3-21" in vGenes and cdr3len == 9 and bool(re.search('^..[DE]...MDV$', cdr3aa)): subset = "[CLL#2]"
				elif "IGHV4-39" in vGenes and cdr3len == 19 and bool(re.search('^A.....SS.W.....WFDP$', cdr3aa)) and identity > 98: subset = "[CLL#8]"
			elif subsetsAnnotation == "major": # based on Agathangelidis et al Blood 2021 definitions excluding IGHJ genes and V-identity (PMID: 32992344)
				if bool(re.search('IGHV[157]', str(vGenes))) and cdr3len == 13 and bool(re.search('^AR.QWL....FDY$', cdr3aa)): subset = "[CLL#1]"
				elif "IGHV3-21" in vGenes and cdr3len == 9 and bool(re.search('^A.[DE]...MDV$', cdr3aa)): subset = "[CLL#2]"
				elif "IGHV4-34" in vGenes and cdr3len == 20 and bool(re.search('^[AV]RG.......[RK]RYYYYGMDV$', cdr3aa)): subset = "[CLL#4]"
				elif "IGHV1-69" in vGenes and cdr3len == 20 and bool(re.search('^AR....GV[IV]...YYYY[GY]MDV$', cdr3aa)): subset = "[CLL#5]"
				elif "IGHV1-69" in vGenes and cdr3len == 21 and bool(re.search('^ARGG.YDY[VI]WGSYR.NDAFDI$', cdr3aa)): subset = "[CLL#6]"
				elif "IGHV4-39" in vGenes and cdr3len == 19 and bool(re.search('^A[RST]...YSSSWY...NWFDP$', cdr3aa)): subset = "[CLL#8]"
				elif "IGHV4-39" in vGenes and cdr3len == 18 and bool(re.search('^A[RST]...YSSSWY...WFDP$', cdr3aa)): subset = "[CLL#8B]"
				elif "IGHV4-39" in vGenes and cdr3len == 22 and bool(re.search('^AR[HD]R.GYCSSTSCYYYYYGMDV$', cdr3aa)): subset = "[CLL#10]"
				elif any(vG in vGenes for vG in ["IGHV1-2", "IGHV1-46"]) and cdr3len == 19 and bool(re.search('^ARD..YYDSSGYY[ST]..FDY$', cdr3aa)): subset = "[CLL#12]"
				elif bool(re.search('IGHV[246]', str(vGenes))) and cdr3len == 10 and bool(re.search('^[AV]RGG.W.FD.$', cdr3aa)): subset = "[CLL#14]"
				elif "IGHV4-34" in vGenes and cdr3len == 24 and bool(re.search('^A.RFYCSG..C....YYYYYG[LM]D[VA]$', cdr3aa)): subset = "[CLL#16]"
				elif "IGHV3-48" in vGenes and cdr3len == 21 and bool(re.search('^AR[DE].DFWSGYY.YYYYY[GY]MDV$', cdr3aa)): subset = "[CLL#31]"
				elif any(vG in vGenes for vG in ["IGHV1-58", "IGHV1-69"]) and cdr3len == 12 and bool(re.search('^A...DFWSGY..$', cdr3aa)): subset = "[CLL#59]"
				elif bool(re.search('IGHV3', str(vGenes))) and cdr3len == 21 and bool(re.search('^A[KR][DE][ST][PL]LVV[PV][AT]AI[FY]YYYYGMDV$', cdr3aa)): subset = "[CLL#64B]"
				elif bool(re.search('IGHV3', str(vGenes))) and cdr3len == 12 and bool(re.search('^A[KR]D....[WY]..DY$', cdr3aa)): subset = "[CLL#73]"
				elif any(vG in vGenes for vG in ["IGHV4-4", "IGHV4-59"]) and cdr3len == 14 and bool(re.search('^[AV]RG[PA][DN].[ST]GW..[FL].Y$', cdr3aa)): subset = "[CLL#77]"
				elif bool(re.search('IGHV[157]', str(vGenes))) and cdr3len == 14 and bool(re.search('^AR.QWL.....FDY$', cdr3aa)): subset = "[CLL#99]"
				elif "IGHV1-18" in vGenes and cdr3len == 10 and bool(re.search('^AR.SGG..[DE].$', cdr3aa)): subset = "[CLL#111]"
				elif "IGHV3-48" in vGenes and cdr3len == 9 and bool(re.search('^AR[DE]......$', cdr3aa)): subset = "[CLL#169]"
				elif "IGHV3-72" in vGenes and cdr3len == 17 and bool(re.search('^[AV]R..YC[ST][SG][TG][TS]CR..[FL]D.$', cdr3aa)): subset = "[CLL#188]"
				elif "IGHV4-34" in vGenes and cdr3len == 17 and bool(re.search('^ARR...W.....D[AG]FD.$', cdr3aa)): subset = "[CLL#201]"
				elif "IGHV3-30" in vGenes and cdr3len == 19 and bool(re.search('^AK[VI]...G.F....YYGMD[VA]$', cdr3aa)): subset = "[CLL#252]"
				elif "IGHV2-5" in vGenes and cdr3len == 17 and bool(re.search('^AHR......W..G.FDY$', cdr3aa)): subset = "[CLL#148B]"
				elif "IGHV1-2" in vGenes and cdr3len == 17 and bool(re.search('^AR.[YL]SGSYYYYYYGMDV$', cdr3aa)): subset = "[CLL#28A]"
				elif "IGHV1-69" in vGenes and cdr3len == 22 and bool(re.search('^A...DIVVVPAA..YYYYGMDV$', cdr3aa)): subset = "[CLL#3C2]"
				elif "IGHV1-69" in vGenes and cdr3len == 22 and bool(re.search('^AR..PDIVVVPAAI.[YR]YYGMDV$', cdr3aa)): subset = "[CLL#3C3]"
				elif "IGHV1-69" in vGenes and cdr3len == 23 and bool(re.search('^A[RST]....DFWSGYYPNYYYYGMDV$', cdr3aa)): subset = "[CLL#7C2]"
				elif "IGHV1-69" in vGenes and cdr3len == 24 and bool(re.search('^A.....[YGD]DFWSGYYPNYYYY[GY]MDV$', cdr3aa)): subset = "[CLL#7D3]"
				elif bool(re.search('IGHV3', str(vGenes))) and cdr3len == 14 and bool(re.search('^ARG..GDY...FD[YIV]$', cdr3aa)): subset = "[CLL#202]"
			if subset != "": i[21] = i[21]+" "+subset
	return(information)

def annotateR110mutation(information):
	for i in information:
		if "IGLV3-21" in i[21] and i[16].endswith("C"):
			i[21] = i[21]+" [R110]"
	return(information)

def predefinedFilter(information, seq, seqDepth, reportOnlyProductive, scoreCutoff, genomeVersion, collapseSequencesSimilarity):

	trip = {} # dict to save passing rearrangements
	kdeCount = 2 # max number of Kde-RSS in filtered file

	for line in information:
		
		pr = 0
		mech = line[1]
		spl_ins = line[22]
		spl_ins_phased = len(line[24].split(",")) # get total number of reads
		original_spl_ins = line[2]+line[3]
		line.append(spl_ins_phased) # append total number of reads to line
		mq = float(line[23].split(" ")[0].replace("NA", "0"))
		phasing_pct = 0 if line[15] == "NA" else float(line[15].split("/")[0])/float(line[15].split(" ")[0].split("/")[1])*100 if line[15].split(" ")[0].split("/")[1] != "0" else 100
		muts_low_confidence = 300 if line[15] == "NA" else float(line[15].split(" - ")[1])
		ident = 100 if line[17] == "NA" else float(line[17])
		nw2 = line[21].split(" - ") # annotated values based on IgBlast
		nw2 = [g.split("*")[0] for g in nw2] # remove alleles and keep first gene 

		# 1) Hard filters
		## Mechanism
		if line[0].startswith(("IGH", "IGL", "TRA", "TRG", "TRD")) and mech != "Deletion": continue # only Deletion considered for IGH/IGL/TRA/TRG/TRD
		if line[0].startswith("IGK") and mech == "Inversion1": continue # only Deletion or Inversion2 considered for IGK
		if line[0].startswith("TRB") and mech != "Deletion" and genomeVersion == "hg38": continue # only Deletion for TRB in hg38
		if line[0].startswith("TRB") and mech == "Inversion2" and genomeVersion == "hg19": continue # only Deletion or Inversion1 for TRB in hg19
		
		## Score
		if spl_ins < scoreCutoff: continue 
	
		## Unknown (N) nucleotides in sequence
		if "Kde" not in line[0] and "RSS" not in line[0] and line[19] != "Partial rearrangement" and line[19] != "Rearrangement without junction coverage":
			if seq == "wgs" and sum(1 for i in line[13] if i == "N")/len(line[13]) > 0.5: continue # remove rearrangement if >50% of the V sequence are "N"s for WGS-derived samples
			if sum(1 for i in line[11] if i == "N")/len(line[11]) > 0.5: continue # remove rearrangement if >50% of the J sequence are "N"s
		
		## Kde-RSS filtering
		if "Kde" in line[0] and "RSS" in line[0]: # remove Kde-RSS if <X split reads and <X insertSize reads
			if seq == "wgs" or seq == "wes":
				if seqDepth == "low" and int(line[2]) < 2 and int(line[3]) < 4: continue
				if seqDepth == "int" and int(line[2]) < 4 and int(line[3]) < 6: continue
				if seqDepth == "high" and int(line[2]) < 8 and int(line[3]) < 8: continue
			else:
				if seqDepth == "low" and int(line[2]) < 2: continue
				if seqDepth == "int" and int(line[2]) < 8: continue
				if seqDepth == "high" and int(line[2]) < 15: continue

		# 2) Soft filters
		if len(trip) == 0:
			trip[line[0]] = line[1:] 
			
		else:
			## a) InsertSize only rearrangements:
			if line[19] == "Rearrangement without junction coverage":
				trip[line[0]] = line[1:]
			
			## b) Partial rearrangements (only consider the one with the highest score)
			elif line[19] == "Partial rearrangement":
				if line[0] not in trip:
					partialJDRearrangement = "yes" if line[0].split(" - ")[1][3] == "D" else "no"
					if partialJDRearrangement == "yes":
						trip[line[0]] = line[1:]
					else:
						Dgene = line[0].split(" - ")[0]
						VgeneFamily = line[0].split(" - ")[1].split("-")[0]
						sameDandVfamily = "no"
						for keyTrip in trip:
							if keyTrip.startswith(Dgene+" - "+VgeneFamily):
								sameDandVfamily = "yes"
								break
						if sameDandVfamily == "no":
							trip[line[0]] = line[1:]					
				elif spl_ins > int(trip[line[0]][21]):
					trip[line[0]] = line[1:]
				else:
					continue

			## c) Kde - RSS:
			elif "Kde" in line[0] and "RSS" in line[0]:
				if line[0] not in trip:
					trip[line[0]] = line[1:]
				elif int(line[2]) == 0 and int(trip[line[0]][1]) == 0 and spl_ins >= int(trip[line[0]][21]): # both 0 split reads, keep new - overwrite previous if >= score
					trip[line[0]] = line[1:]
				elif int(line[2]) > 0 and int(trip[line[0]][1]) == 0: # if new has split reads, keep new, overwrite previous
					trip[line[0]] = line[1:]
				elif int(line[2]) < 2 and int(trip[line[0]][1]) > 1: # if new has <2 split reads, keep old
					continue
				elif int(line[2]) > 1 and int(trip[line[0]][1]) < 2: # if new has 2+ split reads and old <2, keep new, overwrite previous
					trip[line[0]] = line[1:]
				else: # keep both
					if kdeCount == 2:
						line[21] = line[21]+" ("+str(kdeCount)+")"
						trip[line[0]+" ("+str(kdeCount)+")"] = line[1:]
						kdeCount += 1
				
			## d) no Partial and no Kde-RSS:
			else:
				### d1) check if they share >75% of the readNames (ie one with mostly the same insertSize reads), if so, keep only the one with higher score if score is 2x (IGHV1-8 exception)
				filteredByReadNames = False
				keyToUse = ""
				readNamesList = line[24].split(",")
				for keys in trip:
					readNamesListTrip = trip[keys][23].split(",")
					numReadsInCommon = sum([1 if rNameIt in readNamesListTrip else 0 for rNameIt in readNamesList])
					fracReadsShared = numReadsInCommon / min(len(readNamesList), len(readNamesListTrip))
					if fracReadsShared > 0.75:
						filteredByReadNames = True
						keyToUse = keys
						break
				if filteredByReadNames:
					keys = keyToUse
					dict_spl_ins = trip[keys][21]
					dict_phasing_pct = 0 if trip[keys][14] == "NA" else float(trip[keys][14].split("/")[0])/float(trip[keys][14].split(" ")[0].split("/")[1])*100 if trip[keys][14].split(" ")[0].split("/")[1] != "0" else 100
					ts2 = trip[keys][20].split(" - ") # annotated values based on IgBlast
					ts2 = [g.split("*")[0] for g in ts2] # remove alleles and keep first gene 
					# IGHV1-8 gene filtering (not well annotated in hg19/hg38)
					if "IGHV1-8" in nw2 and not "IGHV1-8" in ts2 and phasing_pct >= dict_phasing_pct and spl_ins >= dict_spl_ins*0.25:
						del trip[keys]
					# different scores, keep highest score
					elif spl_ins > dict_spl_ins*2:
						del trip[keys]
					elif spl_ins*2 < dict_spl_ins:
						pr = 1
					else:
						filteredByReadNames = False
				
				if not filteredByReadNames:
					### d2) check if same CDR3 is already annotated...	
					if line[20] != "NA" and line[20] in [trip[keys][19] for keys in trip]:
						for keys in [k for k in trip]: # make list of keys to avoid dictionary changed size during iteration
							if trip[keys][18] == "Partial rearrangement": continue
							dict_spl_ins = trip[keys][21]
							dict_spl_ins_phased = len(trip[keys][23].split(","))
							dict_mq = float(trip[keys][22].split(" ")[0].replace("NA", "0"))
							dict_phasing_pct = 0 if trip[keys][14] == "NA" else float(trip[keys][14].split("/")[0])/float(trip[keys][14].split(" ")[0].split("/")[1])*100 if trip[keys][14].split(" ")[0].split("/")[1] != "0" else 100
							dict_muts_low_confidence = 300 if trip[keys][14] == "NA" else float(trip[keys][14].split(" - ")[1])
							dict_ident = 100 if trip[keys][16] == "NA" else float(trip[keys][16])
							ts2 = trip[keys][20].split(" - ") # annotated values based on IgBlast
							ts2 = [g.split("*")[0] for g in ts2] # remove alleles and keep first gene 
							if line[20] == trip[keys][19]:
								# IGHV1-8 gene filtering (not well annotated in hg19/hg38)
								if "IGHV1-8" in nw2 and not "IGHV1-8" in ts2 and phasing_pct >= dict_phasing_pct and spl_ins >= dict_spl_ins*0.25:
									del trip[keys]
								# other genes
								elif (spl_ins >= dict_spl_ins*0.75 and spl_ins <= dict_spl_ins*1.25) or (ident < 96 and dict_ident < 96 and spl_ins_phased >= dict_spl_ins_phased*0.75 and spl_ins_phased <= dict_spl_ins_phased*1.25): # if similar scores
									if phasing_pct > dict_phasing_pct: # based on phasing
										del trip[keys]
									elif dict_phasing_pct > phasing_pct:
										pr = 1
									elif muts_low_confidence < 3 and dict_muts_low_confidence > 3 and muts_low_confidence < dict_muts_low_confidence*2: # based on number of potential mutations not correctly phased, keep lowest number
										del trip[keys]
									elif muts_low_confidence > 3 and dict_muts_low_confidence < 3 and muts_low_confidence > dict_muts_low_confidence*2:
										pr = 1
									elif mq > 50 and dict_mq < 10: # based on map qual
										del trip[keys]
									elif mq < 10 and dict_mq > 50:
										pr = 1
									elif ident < 96 and dict_ident < 96 and spl_ins_phased > dict_spl_ins_phased:
										del trip[keys]
									elif ident < 96 and dict_ident < 96 and spl_ins_phased < dict_spl_ins_phased:
										pr = 1
									elif spl_ins > dict_spl_ins: # different scores, keep highest score
										del trip[keys]
									else:
										pr = 1
								elif seqDepth == "high" and phasing_pct == 100 and phasing_pct > dict_phasing_pct and ident < 96 and dict_ident < 96 and spl_ins_phased >= dict_spl_ins_phased*0.45: # if high-depth, prioritize phasing over score
									del trip[keys]
								elif seqDepth == "high" and phasing_pct == 100 and phasing_pct > dict_phasing_pct and spl_ins >= dict_spl_ins*0.45: # if high-depth, prioritize phasing over score
									del trip[keys]
								elif ident < 96 and dict_ident < 96 and spl_ins_phased > dict_spl_ins_phased: # different scores, keep highest score
									del trip[keys]
								elif spl_ins > dict_spl_ins: # different scores, keep highest score
									del trip[keys]
								else:
									pr = 1
					
					else:
						### d3) check if exactly the same VDJ is already annotated...
						if line[0] in trip:
							dict_spl_ins = trip[line[0]][21]
							dict_mq = float(trip[line[0]][22].split(" ")[0].replace("NA", "0"))
							dict_cdr3 = trip[line[0]][19]
							dict_phasing_pct = 0 if trip[line[0]][14] == "NA" else float(trip[line[0]][14].split("/")[0])/float(trip[line[0]][14].split(" ")[0].split("/")[1])*100 if trip[line[0]][14].split(" ")[0].split("/")[1] != "0" else 100
							dict_muts_low_confidence = 300 if trip[line[0]][14] == "NA" else float(trip[line[0]][14].split(" - ")[1])
							if line[20] != "NA" and dict_cdr3 == "NA": # if the same VDJ, keep the one with info in CDR3 aa seq
								del trip[line[0]]
							elif line[20] == "NA" and dict_cdr3 != "NA":
								pr = 1
							elif spl_ins >= dict_spl_ins*0.75 and spl_ins <= dict_spl_ins*1.25: # if similar scores
								if phasing_pct > dict_phasing_pct: # based on phasing
									del trip[line[0]]
								elif dict_phasing_pct > phasing_pct:
									pr = 1
								elif muts_low_confidence < 3 and dict_muts_low_confidence > 3 and muts_low_confidence < dict_muts_low_confidence*2: # based on number of potential mutations not correctly phased, keep lowest number
									del trip[line[0]]
								elif muts_low_confidence > 3 and dict_muts_low_confidence < 3 and muts_low_confidence > dict_muts_low_confidence*2:
									pr = 1
								elif mq > 50 and dict_mq < 10: # based on map qual
									del trip[line[0]]
								elif mq < 10 and dict_mq > 50:
									pr = 1
								elif spl_ins > dict_spl_ins: # different scores, keep highest score
									del trip[line[0]]
								else:
									pr = 1
							elif seqDepth == "high" and phasing_pct == 100 and phasing_pct > dict_phasing_pct and spl_ins >= dict_spl_ins*0.45: # if high-depth, prioritize phasing over score
								del trip[line[0]]
							elif spl_ins > dict_spl_ins: # different scores, keep highest score
								del trip[line[0]]
							else:
								pr = 1
						
						else: # force no exact match to enter here...
							### d4) check if a "not exact VDJ" is already annotated...
							for keys in [k for k in trip]: # make list of keys to avoid dictionary changed size during iteration
								if trip[keys][18] == "Partial rearrangement": continue
								ts = keys.split(" - ") # annotated values
								nw = line[0].split(" - ") # new value
								common = set(ts).intersection(nw) # intersection between values in dict and value analysed
								ts2 = trip[keys][20].split(" - ") # annotated values based on IgBlast
								ts2 = [g.split("*")[0] for g in ts2] # remove alleles and keep first gene 
								common2 = set(ts2).intersection(nw2) # intersection between values in dict and value analysed
								dict_spl_ins = trip[keys][21]
								dict_spl_ins_phased = len(trip[keys][23].split(","))
								dict_original_spl_ins = trip[keys][1]+trip[keys][2]
								dict_mq = float(trip[keys][22].split(" ")[0].replace("NA", "0"))
								dict_cdr3 = trip[keys][19]
								dict_phasing_pct = 0 if trip[keys][14] == "NA" else float(trip[keys][14].split("/")[0])/float(trip[keys][14].split(" ")[0].split("/")[1])*100 if trip[keys][14].split(" ")[0].split("/")[1] != "0" else 100
								dict_muts_low_confidence = 300 if trip[keys][14] == "NA" else float(trip[keys][14].split(" - ")[1])
								#### Specific conditions after manual review of multiple analyses:
								##### A: 2 genes in common based on breaks + CDR3 inside
								##### B: 2 genes in common + highly similary CDR3s
								##### C: same IGHJ breaks + (CDR3 inside or highly similary CDR3s)
								##### D: same J-V, one without split-reads
								##### E: 2 genes in common -also considering IgBlast annotation- and one with CDR3 = NA
								##### F: same as D but considering V-gene family in the comparison instead of V-gene
								##### G: 2 genes in common -also considering IgBlast annotation- and same V seq or V seq within or V seq same length and similar
								##### H: 2 genes in common -also considering IgBlast annotation-, one with 3x spl_ins
								##### I: 2 genes in common, one with 2x spl_ins_phased and 2x original_spl_ins
								##### J: for IGH: same J gene or IGH J break or V break, same V family, and highly similar CDR3
								condiA = len(common) == 2 and line[20] != "NA" and dict_cdr3 != "NA" and (line[20] in dict_cdr3 or dict_cdr3 in line[20])
								condiB = (len(common) == 2 or len(common2) >= 2) and line[20] != "NA" and dict_cdr3 != "NA" and abs(len(line[20])-len(dict_cdr3)) <= 1 and SequenceMatcher(None, line[20], dict_cdr3).ratio() > 0.8
								condiC = nw[0].startswith("IGHJ") and line[4] == trip[keys][3] and line[5] == trip[keys][4] and line[20] != "NA" and dict_cdr3 != "NA" and ((line[20] in dict_cdr3 or dict_cdr3 in line[20]) or (abs(len(line[20])-len(dict_cdr3)) <= 1 and SequenceMatcher(None, line[20], dict_cdr3).ratio() > 0.9))
								condiD = ((nw[0] == ts[0] and nw[-1] == ts[-1]) or (nw2[0] == ts2[0] and nw2[-1] == ts2[-1])) and ((line[2] == 0 and line[6] == 0 and line[9] == 0) or (trip[keys][1] == 0 and trip[keys][5] == 0 and trip[keys][8] == 0))
								condiE = (len(common) == 2 or len(common2) >= 2) and (line[20] == "NA" or dict_cdr3 == "NA")
								condiF = ((nw[0] == ts[0] and nw[-1].split("-")[0] == ts[-1].split("-")[0]) or (nw2[0] == ts2[0] and nw2[-1].split("-")[0] == ts2[-1].split("-")[0])) and ((line[2] == 0 and line[6] == 0 and line[9] == 0) or (trip[keys][1] == 0 and trip[keys][5] == 0 and trip[keys][8] == 0))
								condiG = (len(common) == 2 or len(common2) >= 2) and (line[13] in trip[keys][12] or trip[keys][12] in line[13] or CustomSequenceMatcher(line[13], trip[keys][12], "nt") > collapseSequencesSimilarity)
								condiH = ((nw[0] == ts[0] and nw[-1] == ts[-1]) or (nw2[0] == ts2[0] and nw2[-1] == ts2[-1])) and (dict_spl_ins*3 < spl_ins or spl_ins*3 < dict_spl_ins)
								condiI = (len(common) == 2 or len(common2) >= 2) and ( (dict_spl_ins_phased*2 < spl_ins_phased and dict_original_spl_ins*2 < original_spl_ins) or (spl_ins_phased*2 < dict_spl_ins_phased and original_spl_ins*2 < dict_original_spl_ins))
								condiJ = nw[0].startswith("IGHJ") and (nw[0] == ts[0] or nw2[0] == ts2[0] or line[5] == trip[keys][4] or line[7] == trip[keys][6]) and (nw[-1].split("-")[0] == ts[-1].split("-")[0] or nw2[-1].split("-")[0] == ts2[-1].split("-")[0]) and (line[20] in dict_cdr3 or dict_cdr3 in line[20] or CustomSequenceMatcher(line[20], dict_cdr3, "aa") > 0.6)
								if condiA or condiB or condiC or condiD or condiE or condiF or condiG or condiH or condiI or condiJ:
									if dict_spl_ins*3 < spl_ins: # check if one has more split/insert reads
										del trip[keys]
									elif spl_ins*3 < dict_spl_ins:
										pr = 1
									elif dict_spl_ins_phased*2 < spl_ins_phased and dict_original_spl_ins*2 < original_spl_ins: # check if one has more reads (including phased reads) and more initially mapped split/insert reads
										del trip[keys]
									elif spl_ins_phased*2 < dict_spl_ins_phased and original_spl_ins*2 < dict_original_spl_ins:
										pr = 1
									elif line[20] != "NA" and dict_cdr3 == "NA": # keep the one with info in CDR3 aa seq
										del trip[keys]
									elif line[20] == "NA" and dict_cdr3 != "NA":
										pr = 1
									elif len(ts) > len(nw) and len(ts2) > len(nw2): # if the one annotated has len=3 (VDJ) and the new one 2 (VJ), keep the one annotated
										pr = 1
									elif len(ts) < len(nw) and len(ts2) < len(nw2): # if the other way around... keep the new one
										del trip[keys]
									elif (line[13] in trip[keys][12] or trip[keys][12] in line[13]) and dict_spl_ins < spl_ins: # if same V-seq, keep only the one with higher score
										del trip[keys]
									elif (line[13] in trip[keys][12] or trip[keys][12] in line[13]) and spl_ins < dict_spl_ins:
										pr = 1
									elif line[5] != trip[keys][4] and line[7] != trip[keys][6]: # if different breakpoints, keep both if similar score
										if spl_ins*0.5 > dict_spl_ins:
											del trip[keys]
										elif spl_ins*1.5 < dict_spl_ins:
											pr = 1
									elif spl_ins >= dict_spl_ins*0.75 and spl_ins <= dict_spl_ins*1.25: # if similar scores
										if phasing_pct > dict_phasing_pct: # based on phasing
											del trip[keys]
										elif dict_phasing_pct > phasing_pct:
											pr = 1
										elif muts_low_confidence < 3 and dict_muts_low_confidence > 3 and muts_low_confidence < dict_muts_low_confidence*2: # based on number of potential mutations not correctly phased, keep lowest number
											del trip[keys]
										elif muts_low_confidence > 3 and dict_muts_low_confidence < 3 and muts_low_confidence > dict_muts_low_confidence*2:
											pr = 1
										elif mq > 50 and dict_mq < 10: # based on map qual
											del trip[keys]
										elif mq < 10 and dict_mq > 50:
											pr = 1
										elif line[5] != trip[keys][4] or line[7] != trip[keys][6]: # if similar score, same phasing, similar mq, and at least one different breakpoint, keep both 
											pr = 0
										elif spl_ins > dict_spl_ins: # different scores, keep highest score
											del trip[keys]
										else:
											pr = 1
									elif seqDepth == "high" and phasing_pct == 100 and phasing_pct > dict_phasing_pct and spl_ins >= dict_spl_ins*0.65: # if high-depth, prioritize phasing over score
										del trip[keys]
									elif spl_ins > dict_spl_ins: # different scores, keep highest score
										del trip[keys]					
									else:
										pr = 1
									break # force only one match in condis
								
								#### 1 gene in common (IGLV/IGKV/TRAV/TRGV genes)
								elif len(common) == 1 and ( nw[1].startswith("IGLV") or nw[1].replace("D", "").startswith("IGKV") or nw[1].startswith("TRAV") or nw[1].startswith("TRGV") ):
									if nw[1].replace("D", "") == ts[1].replace("D", "") and ts[0] != "IGKKde":
										if spl_ins > dict_spl_ins:
											del trip[keys]
										else:
											pr = 1
										break # force only one match
					
				### d5) add if needed
				if pr == 0:
					trip[line[0]] = line[1:]
	
	## Report Only productive
	if reportOnlyProductive == "no":
		return(trip)
	else:
		prodTrip = {}
		for ky in trip:
			line = trip[ky]
			if "Productive" in line[18] or "Potentially productive" in line[18] or line[18] == "Rearrangement without junction coverage" or line[18] == "Complete rearrangement (no junction found)":
				prodTrip[ky] = line
			elif reportOnlyProductive == "oof" and line[18] == "Unproductive (stop codons)":
				prodTrip[ky] = line
		return(prodTrip)

def classSwitchAnalysis(wkDir, data, annot_table_JV, bedFile, baseq, chromGene, bamT, bamN, pathToSamtools, tumorPurity, seq, scoreCutoffCSR, errLogMpileup):
	class_switch = []
	class_switch_filt = []
	reductionMeans  = []
	
	for key in data:
		kGenes = key.split(" - ")[0]+" - "+key.split(" - ")[1]
		kReadtype = key.split(" - ")[2]
		
		score = round(data[key] / tumorPurity, 1)
		
		# hard filter: keep only CSR if supported by >= 5 insertSize reads (by default) and by "Deletion"
		hardCutoff = 4 if scoreCutoffCSR > 4 else scoreCutoffCSR
		if score < hardCutoff or kReadtype != "Deletion": continue

		# add readNames and maqQual
		gene1 = kGenes.split(" - ")[0]
		gene2 = kGenes.split(" - ")[1]
		readNames = []
		mapQual = []
		ANNOT_TABLE_JV = open(annot_table_JV, "r")
		for j in ANNOT_TABLE_JV:
			w = j.rstrip("\n").split("\t")
			if gene1 == w[18] and gene2 == w[19]:
				readNames.append(w[0])
				mapQual.append(int(w[4]))
		readNames = ",".join(set(readNames))
		mapQual = str(round(mean(mapQual),1))+" ("+str(min(mapQual))+"-"+str(max(mapQual))+")"
		numReads = len(readNames.split(","))

		# study coverage and soft filter if seq == wgs:
		if seq == "wgs":
			isotypye = kGenes.split(" - ")[0]

			VDJ = open(bedFile, "r")
			for k in VDJ:
				v = k.rstrip("\n").split("\t")
				if isotypye == v[3]:
					startA = int(v[1])-1500
					endA = int(v[1])
					startB = int(v[2])
					endB = int(v[2])+1500
					break
			VDJ.close()		
			
			covs = {}
			for i in ["A", "B"]:
				if i == "A":
					st = startA
					en = endA							
				else:
					st = startB
					en = endB
				
				mpileupFile = wkDir+"/tmp/"+bamT.split("/")[-1].replace(".bam", "_output_mpileup.tsv")
				subprocess.call(pathToSamtools+"samtools mpileup -d 0 -B -Q "+baseq+" -r "+chromGene+":"+str(st)+"-"+str(en)+" "+bamT+ " > "+mpileupFile+" 2>> "+errLogMpileup, shell=True)					
				
				pos = st
				lst = []
				seq = open(mpileupFile, "r")
				for sLine in seq:
					sList = sLine.rstrip("\n").split("\t")
					while pos < int(sList[1]):
						lst.append(0)
						pos += 1
					lst.append(int(sList[3]))
					pos += 1
				seq.close()
				
				while pos <= en:
					lst.append(0)
					pos += 1	
					
				covs[i] = lst
			
			meanA = round(mean(covs["A"]), 3)
			meanB = round(mean(covs["B"]), 3)
			
			meanNormA = 1
			meanNormB = 1

			# if normal bam file available, substract depth in normal to correct for different distribution of coverage and correct means to positive
			if bamN is not None:
						
				for i in ["A", "B"]:
					
					if i == "A":
						st = startA
						en = endA							
					else:
						st = startB
						en = endB
					
					mpileupFile = wkDir+"/tmp/"+bamN.split("/")[-1].replace(".bam", "_output_mpileup.tsv")
					subprocess.call(pathToSamtools+"samtools mpileup -d 0 -B -Q "+baseq+" -r "+chromGene+":"+str(st)+"-"+str(en)+" "+bamN+ " > "+mpileupFile+" 2>> "+errLogMpileup, shell=True)					
					
					pos = st
					idx = 0
					forMeanNorm = []
					seq = open(mpileupFile, "r")
					for sLine in seq:
						sList = sLine.rstrip("\n").split("\t")
						while pos < int(sList[1]):
							covs[i][idx] = covs[i][idx] - 0
							forMeanNorm.append(0)
							idx += 1
							pos += 1
						covs[i][idx] = covs[i][idx] - int(sList[3])
						forMeanNorm.append(int(sList[3]))
						idx += 1
						pos += 1
					seq.close()
					
					if i == "A":
						meanNormA = mean(forMeanNorm)
					else:
						meanNormB = mean(forMeanNorm)
			
			
			# adjusted meanA to ratio means in normal
			meanA = round(meanA/(meanNormA/meanNormB), 3)
			
			# wilcoxon
			statistic, pvalue = stats.wilcoxon(covs["A"], covs["B"])
			
			# reduction adjusted means adjusted by tumor purity
			reductionMean = round(100-(meanB/meanA*100), 3) / tumorPurity
					
			# hard filter:
			if reductionMean > 0:
				class_switch.append([ kGenes.split(" - ")[0], kReadtype, score, mapQual, numReads, meanA, meanB, pvalue, reductionMean, readNames ])
			
			# pre-defined soft filer:
			if meanA > 8 and reductionMean >= 30 and pvalue < 0.0000000001: # 1e-10
				if (score >= scoreCutoffCSR and reductionMean >= 60) or (score >= (scoreCutoffCSR*1.5) and reductionMean >= 30):
					class_switch_filt.append([ kGenes.split(" - ")[0], kReadtype, score, mapQual, numReads, meanA, meanB, pvalue, reductionMean, readNames ])
					reductionMeans.append(reductionMean)

		# just filter based on score if seq != wgs:
		else:
			class_switch.append([kGenes.split(" - ")[0], kReadtype, score, mapQual, numReads, "NA", "NA", "NA", "NA", readNames])
			if score >= scoreCutoffCSR:
				class_switch_filt.append([kGenes.split(" - ")[0], kReadtype, score, mapQual, numReads, "NA", "NA", "NA", "NA", readNames])
				reductionMeans.append(score)
	
	class_switch = sorted(class_switch, key=operator.itemgetter(1), reverse=True)
	
	return(class_switch, class_switch_filt, reductionMeans)

def getIgTranslocations(wkDir, genomeVersion, inputsFolder, pathToSamtools, threadsForSamtools, coordsToSubset, bamT, bamN, pairedMode, chrom, geneToAnalyze, tumorPurity, mntonco, mntoncoPass, vafOnco, mnnonco, mapqOnco, mncPoN, genesOncoIg, customGenesOncoIg, genesOncoIgDistance, customGenesOncoIgDistance, reportReadNames, errLogMpileup):
	
	if genomeVersion == "hg19":
		chrom14_IGH = [106052774, 107288051] # IGH region 
		chrom2_IGK = [89131589, 90274600] # IGK region
		chrom22_IGL = [22380000, 23266000] # IGL region
		chrom14_TRA_TRD = [22069991, 23034042] # TRA region (TRD within TRA)
		chrom7_TRB = [141979017, 142531084] # TRB region
		chrom7_TRG = [38272981, 38427770] # TRG region
	else:
		chrom14_IGH = [105576937, 106889844] # IGH region 
		chrom2_IGK = [88822278, 90245370] # IGK region
		chrom22_IGL = [22016076, 22932913] # IGL region
		chrom14_TRA_TRD = [21611904, 22562132] # TRA region (TRD within TRA)
		chrom7_TRB = [142289011, 142823287] # TRB region
		chrom7_TRG = [38230024, 38378055] # TRG region

	chroms = [chrom+str(i) for i in range(1,22)]+[chrom+"X", chrom+"Y", "="] # chroms considered ("=" to consider deletions/inversions/gains within the same chromosome)
	
	samT = wkDir+"/tmp/"+bamT.split("/")[-1].replace(".bam", ".sam")
	comms = pathToSamtools+"samtools view -@ "+threadsForSamtools+" -q "+mapqOnco+" "+bamT+" "+coordsToSubset+" > "+samT
	subprocess.call(comms, shell=True)
	
	if bamN is not None and pairedMode == "paired":
		samN = wkDir+"/tmp/"+bamN.split("/")[-1].replace(".bam", ".sam")
		comms = pathToSamtools+"samtools view -@ "+threadsForSamtools+" -q "+mapqOnco+" "+bamN+" "+coordsToSubset+" > "+samN
		subprocess.call(comms, shell=True)

	# 1. annotate potential 1-read translocations
	dicForTranslocations = {} 
	dicForTranslocations[chrom+"2"] = {}
	dicForTranslocations[chrom+"7"] = {}
	dicForTranslocations[chrom+"14"] = {}
	dicForTranslocations[chrom+"22"] = {}

	samfile = open(samT, "r")
	for i in samfile:
		
		w = i.rstrip("\n").split("\t")
		sa = []
		for x in (w[11:]): # get SA:... after qualities
			if x.startswith("SA:Z"):
				sa.append(x)
		if len(sa) != 0: w = w[:10]+[sa[0]]
		else: w = w[:10]+["NA"]

		if w[5] != "*" and w[6] in chroms:
			if w[6] == "=" and w[10] == "NA":
				if abs(int(w[3]) - int(w[7])) > 10000: # for inversion, deletions, gains
					if w[2] == chrom+"14" and int(w[7]) >= int(chrom14_IGH[0]) and int(w[7]) <= int(chrom14_IGH[1]): continue
					elif w[2] == chrom+"22" and int(w[7]) >= int(chrom22_IGL[0]) and int(w[7]) <= int(chrom22_IGL[1]): continue
					elif w[2] == chrom+"2" and int(w[7]) >= int(chrom2_IGK[0]) and int(w[7]) <= int(chrom2_IGK[1]): continue
					elif w[2] == chrom+"14" and int(w[7]) >= int(chrom14_TRA_TRD[0]) and int(w[7]) <= int(chrom14_TRA_TRD[1]): continue
					elif w[2] == chrom+"7" and int(w[7]) >= int(chrom7_TRB[0]) and int(w[7]) <= int(chrom7_TRB[1]): continue
					elif w[2] == chrom+"7" and int(w[7]) >= int(chrom7_TRG[0]) and int(w[7]) <= int(chrom7_TRG[1]): continue
				else: continue
				
			flgBin = flagToCustomBinary(w[1])
			strands = "regular" 
			
			# break 1
			inChrom = w[2]
			split1 = re.findall(r'[A-Za-z]|[0-9]+', w[5])
			two1 = [split1[x:x+2] for x in range(0, len(split1),2)]
			if any("H" in sub for sub in two1): continue # remove hard-clipped reads
			posInChrom = w[3]
			strandInChrom = "-" if "1" == flgBin[4] else "+"
			if any("S" in sub for sub in two1): # adjust strandInChrom if needed
				if strandInChrom == "+" and next(i for i, val in enumerate(two1, 0) if "S" in val) < next(i for i, val in enumerate(two1, 0) if "M" in val):
					strandInChrom = "-"
					strands = "strand1_pos_is_neg"
				if strandInChrom == "-" and next(i for i, val in enumerate(two1, 0) if "M" in val) < next(i for i, val in enumerate(two1, 0) if "S" in val):
					strandInChrom = "+"
					strands = "strand1_neg_is_pos"
			if strandInChrom == "+": # adjust posInChrom if needed
				if not any("S" in sub for sub in two1):
					posInChrom = str(int(posInChrom) + sum([int(i[0]) for i in two1 if "M" in i or "D" in i]) - 1)
				elif next(i for i, val in enumerate(two1, 0) if "M" in val) < next(i for i, val in enumerate(two1, 0) if "S" in val): ## check if M before S
					posInChrom = str(int(posInChrom) + sum([int(i[0]) for i in two1 if "M" in i or "D" in i]) - 1)

			# break 2
			outChrom = w[6].replace("=", inChrom)
			posOutChrom = w[7]
			readType = "paired"
			strandOutChrom = "-" if "1" == flgBin[5] else "+" # strand from insert size to get orientation of the translocation
			strandOutChromSA = ""
			if w[10].startswith("SA:Z") and w[10].split(":")[2].split(",")[0] in chroms: # get from splits if any
				readType = "split"
				outChrom = w[10].split(":")[2].split(",")[0]
				posOutChrom = w[10].split(":")[2].split(",")[1]
				strandOutChromSA = w[10].split(":")[2].split(",")[2]
				split2 = re.findall(r'[A-Za-z]|[0-9]+', w[10].split(":")[2].split(",")[3])
				two2 = [split2[x:x+2] for x in range(0, len(split2),2)]
				if strandOutChromSA == "+" and next(i for i, val in enumerate(two2, 0) if "S" in val) < next(i for i, val in enumerate(two2, 0) if "M" in val):
					strandOutChromSA = "-"
					strands = "strand2_pos_is_neg"
				if strandOutChromSA == "-" and next(i for i, val in enumerate(two2, 0) if "M" in val) < next(i for i, val in enumerate(two2, 0) if "S" in val):
					strandOutChromSA = "+"
					strands = "strand2_neg_is_pos"
				if next(i for i, val in enumerate(two2, 0) if "M" in val) < next(i for i, val in enumerate(two2, 0) if "S" in val):  ## check if M before S
					posOutChrom = str(int(posOutChrom) + sum([int(i[0]) for i in two2 if "M" in i or "D" in i]) - 1)
			
			# get N-nucleotides
			nNucleotides = "NA"
			if w[10].startswith("SA:Z") and w[10].split(":")[2].split(",")[0] in chroms:
				if strandInChrom == "+" and strandOutChromSA == "+" and any("M" in sub for sub in two1):
					lastM1  = len(two1) - next(i for i, val in enumerate(reversed(two1), 0) if "M" in val) - 1
					firstM2 = next(i for i, val in enumerate(two2, 0) if "M" in val)
					matches1 = sum([int(val[0]) for i, val in enumerate(two1, 0) if i <= lastM1])
					if strands == "regular": 
						nonMatches2 = sum([int(val[0]) for i, val in enumerate(two2, 0) if i < firstM2])
						nNucleotides = w[9][matches1:nonMatches2] if matches1 < nonMatches2 else "None"
					elif strands == "strand1_neg_is_pos" or strands == "strand2_neg_is_pos":
						nonMatches2 = sum([int(val[0]) for i, val in enumerate(two2, 0) if i > firstM2])
						nNucleotides = w[9][matches1:nonMatches2] if matches1 < nonMatches2 else "None"
					else: nNucleotides = "NA"
				elif strandInChrom == "-" and strandOutChromSA == "-" and any("S" in sub for sub in two1):
					firstM1 = next(i for i, val in enumerate(two1, 0) if "M" in val)
					firstS2 = next(i for i, val in enumerate(two2, 0) if "S" in val)
					nonmatches1 = sum([int(val[0]) for i, val in enumerate(two1, 0) if i < firstM1])
					if strands == "regular": 
						matches2 = sum([int(val[0]) for i, val in enumerate(two2, 0) if i < firstS2])
						nNucleotides = w[9][matches2:nonmatches1] if matches2 < nonmatches1 else "None"
					elif strands == "strand1_pos_is_neg" or strands == "strand2_pos_is_neg":
						matches2 = sum([int(val[0]) for i, val in enumerate(two2, 0) if i > firstS2])
						nNucleotides = w[9][matches2:nonmatches1] if matches2 < nonmatches1 else "None"
					else: nNucleotides = "NA"
				elif strandInChrom == "+" and strandOutChromSA == "-" and any("M" in sub for sub in two1):
					lastM1  = len(two1) - next(i for i, val in enumerate(reversed(two1), 0) if "M" in val) - 1
					lastM2 = len(two2) - next(i for i, val in enumerate(reversed(two2), 0) if "M" in val) - 1
					matches1 = sum([int(val[0]) for i, val in enumerate(two1, 0) if i <= lastM1])
					if strands == "regular":
						nonMatches2 = sum([int(val[0]) for i, val in enumerate(two2, 0) if i > lastM2])
						nNucleotides = w[9][matches1:nonMatches2] if matches1 < nonMatches2 else "None"
					elif strands == "strand1_neg_is_pos" or strands == "strand2_pos_is_neg":
						nonMatches2 = sum([int(val[0]) for i, val in enumerate(two2, 0) if i < lastM2])
						nNucleotides = w[9][matches1:nonMatches2] if matches1 < nonMatches2 else "None"
					else: nNucleotides = "NA"
				elif strandInChrom == "-" and strandOutChromSA == "+" and any("S" in sub for sub in two1):
					lastS1  = len(two1) - next(i for i, val in enumerate(reversed(two1), 0) if "S" in val) - 1
					lastS2  = len(two2) - next(i for i, val in enumerate(reversed(two2), 0) if "S" in val) - 1
					nonmatches1 = sum([int(val[0]) for i, val in enumerate(two1, 0) if i <= lastS1])
					if strands == "regular": 
						matches2 = sum([int(val[0]) for i, val in enumerate(two2, 0) if i > lastS2])
						nNucleotides = w[9][matches2:nonmatches1] if matches2 < nonmatches1 else "None"
					elif strands == "strand1_pos_is_neg" or strands == "strand2_neg_is_pos":
						matches2 = sum([int(val[0]) for i, val in enumerate(two2, 0) if i < lastS2])
						nNucleotides = w[9][matches2:nonmatches1] if matches2 < nonmatches1 else "None"
					else: nNucleotides = "NA"
				else:
					nNucleotides = "NA"

				if len(nNucleotides) > len(w[9])*0.75: continue # potential artefact or situation not considered
			
			if strandOutChromSA != "": strandOutChrom = strandOutChromSA

			if outChrom in dicForTranslocations[inChrom]:
				dicForTranslocations[inChrom][outChrom].append([posInChrom, strandInChrom, posOutChrom, strandOutChrom, nNucleotides, readType, w[4], w[0]])
			else:
				dicForTranslocations[inChrom][outChrom] = [[posInChrom, strandInChrom, posOutChrom, strandOutChrom, nNucleotides, readType, w[4], w[0]]]

	samfile.close()
	
	# 2. Merge individual one-read translocations into potential translocations (kep only if number of reads (ie score) > mntonco)
	translocations = {}
	translocations[chrom+"2"] = {}
	translocations[chrom+"7"] = {}
	translocations[chrom+"14"] = {}
	translocations[chrom+"22"] = {}

	position1 = list()
	strand1 = ""
	position2 = list()
	strand2 = ""
	nNucleotidesList = list()
	readTypeList = list()
	mapQualList = list()
	readNameList = list()
	for key1 in dicForTranslocations:
		for key2 in dicForTranslocations[key1]:
			for item in dicForTranslocations[key1][key2]:
				alreadyConsidered = "no"
				
				# check if new one-read translocation could be added to an already merged potential translocation:
				if key2 in translocations[key1]:
					for item2 in translocations[key1][key2]:
						if ( abs(int(item[0]) - int(item2[1])) < 200 or abs(int(item[0]) - int(item2[2])) < 200 ) and item[1] == item2[3] and ( abs(int(item[2]) - int(item2[5])) < 1000 or abs(int(item[2]) - int(item2[6])) < 1000 )  and item[3] == item2[7]:
							# check if same readName and readType already considered (overlaping R1 and R2)
							if item[7] in item2[12].split(","):
								indexOfRead = item2[12].split(",").index(item[7])
								if item2[10].split("-")[indexOfRead] == item[5]:
									alreadyConsidered = "yes"
									break
							item2[1] = str(min([int(item2[1]), int(item[0])]))
							item2[2] = str(max([int(item2[2]), int(item[0])]))
							item2[5] = str(min([int(item2[5]), int(item[2])]))
							item2[6] = str(max([int(item2[6]), int(item[2])]))
							item2[8] = item2[8]+1
							item2[9] = item2[9]+"-"+item[4]
							item2[10] = item2[10]+"-"+item[5]
							item2[11] = item2[11]+"-"+item[6]
							item2[12] = item2[12]+","+item[7]
							alreadyConsidered = "yes"
							break
				
				# if not added before, initialize new potential translocation:
				if alreadyConsidered == "no":
					if position1 == list():
						position1.append(int(item[0]))
						strand1 = item[1]
						position2.append(int(item[2]))
						strand2 = item[3]
						nNucleotidesList.append(item[4])
						readTypeList.append(item[5])
						mapQualList.append(item[6])
						readNameList.append(item[7])

					elif min(abs(p1 - int(item[0])) for p1 in position1) < 200 and strand1 == item[1] and min(abs(p2 - int(item[2])) for p2 in position2) < 1000 and strand2 == item[3]:
						# check if same readName and readType already considered (overlaping R1 and R2 as split reads or two insertSize/paired reads)
						if item[7] in readNameList:
							indexOfRead = readNameList.index(item[7])
							if readTypeList[indexOfRead] == item[5]:
								alreadyConsidered = "yes"
						if alreadyConsidered == "no":
							position1.append(int(item[0]))
							position2.append(int(item[2]))
							nNucleotidesList.append(item[4])
							readTypeList.append(item[5])
							mapQualList.append(item[6])
							readNameList.append(item[7])

					else:
						if key2 in translocations[key1]: translocations[key1][key2].append([key1, str(min(position1)), str(max(position1)), strand1, key2, str(min(position2)), str(max(position2)), strand2, len(readNameList), "-".join(nNucleotidesList), "-".join(readTypeList), "-".join(mapQualList), ",".join(readNameList), 0 if bamN is not None and pairedMode == "paired" else "NA"]) # 0 will be the count in normal
						else: translocations[key1][key2] = [ [key1, str(min(position1)), str(max(position1)), strand1, key2, str(min(position2)), str(max(position2)), strand2, len(readNameList), "-".join(nNucleotidesList), "-".join(readTypeList), "-".join(mapQualList), ",".join(readNameList), 0 if bamN is not None and pairedMode == "paired" else "NA"] ]

						position1 = [int(item[0])]
						strand1 = item[1]
						position2 = [int(item[2])]
						strand2 = item[3]
						nNucleotidesList = [item[4]]
						readTypeList = [item[5]]
						mapQualList = [item[6]]
						readNameList = [item[7]]
						
			# if no more positions in second chrom, end iteration and reset:
			if key2 in translocations[key1]: translocations[key1][key2].append([key1, str(min(position1)), str(max(position1)), strand1, key2, str(min(position2)), str(max(position2)), strand2, len(readNameList), "-".join(nNucleotidesList), "-".join(readTypeList), "-".join(mapQualList), ",".join(readNameList), 0 if bamN is not None and pairedMode == "paired" else "NA"])
			else: translocations[key1][key2] = [ [key1, str(min(position1)), str(max(position1)), strand1, key2, str(min(position2)), str(max(position2)), strand2, len(readNameList), "-".join(nNucleotidesList), "-".join(readTypeList), "-".join(mapQualList), ",".join(readNameList), 0 if bamN is not None and pairedMode == "paired" else "NA"] ]
			
			position1 = list()
			strand1 = ""
			position2 = list()
			strand2 = ""
			nNucleotidesList = list()
			readTypeList = list()
			mapQualList = list()
			readNameList = list()
	
	# filter based on mntonco cutoff
	translocationsFiltered = {}
	translocationsFiltered[chrom+"2"] = {}
	translocationsFiltered[chrom+"7"] = {}
	translocationsFiltered[chrom+"14"] = {}
	translocationsFiltered[chrom+"22"] = {}
	for key1 in translocations:
		for key2 in translocations[key1]:
			for item in translocations[key1][key2]: 				
				if item[8] >= mntonco:
					if key2 in translocationsFiltered[key1]: translocationsFiltered[key1][key2].append(item)
					else: translocationsFiltered[key1][key2] = [ item ]
	
	# 3. Annotate in normal
	if bamN is not None and pairedMode == "paired":
		readNamesUsedInPoN = [] # to avoid counting R1 and R2 twice
		samfile = open(samN, "r")
		for i in samfile:
			
			w = i.rstrip("\n").split("\t")
			if w[0] in readNamesUsedInPoN: continue
			
			sa = []
			for x in (w[11:]): # get SA:... after qualities
				if x.startswith("SA:Z"):
					sa.append(x)
			if len(sa) != 0: w = w[:10]+[sa[0]]
			else: w = w[:10]+["NA"]
			
			if w[5] != "*" and w[6] in chroms: 
				if w[6] == "=" and w[10] == "NA":
					if abs(int(w[3]) - int(w[7])) > 8000: # 8000 instead of 10000 just to be more permessive in the normal...
						if w[2] == chrom+"14" and int(w[7]) >= int(chrom14_IGH[0]) and int(w[7]) <= int(chrom14_IGH[1]): continue
						elif w[2] == chrom+"22" and int(w[7]) >= int(chrom22_IGL[0]) and int(w[7]) <= int(chrom22_IGL[1]): continue
						elif w[2] == chrom+"2" and int(w[7]) >= int(chrom2_IGK[0]) and int(w[7]) <= int(chrom2_IGK[1]): continue
						elif w[2] == chrom+"14" and int(w[7]) >= int(chrom14_TRA_TRD[0]) and int(w[7]) <= int(chrom14_TRA_TRD[1]): continue
						elif w[2] == chrom+"7" and int(w[7]) >= int(chrom7_TRB[0]) and int(w[7]) <= int(chrom7_TRB[1]): continue
						elif w[2] == chrom+"7" and int(w[7]) >= int(chrom7_TRG[0]) and int(w[7]) <= int(chrom7_TRG[1]): continue
					else: continue
				
				flgBin = flagToCustomBinary(w[1])
				split1 = re.findall(r'[A-Za-z]|[0-9]+', w[5])
				two1 = [split1[x:x+2] for x in range(0, len(split1),2)]
				if any("H" in sub for sub in two1): continue # remove hard-clipped reads

				# break 1
				inChrom =  w[2]
				posInChrom = int(w[3])
				strandInChrom = "-" if "1" == flgBin[4] else "+"					
				if any("S" in sub for sub in two1): # adjust strandInChrom if needed
					if strandInChrom == "+" and next(i for i, val in enumerate(two1, 0) if "S" in val) < next(i for i, val in enumerate(two1, 0) if "M" in val):
						strandInChrom = "-"
					if strandInChrom == "-" and next(i for i, val in enumerate(two1, 0) if "M" in val) < next(i for i, val in enumerate(two1, 0) if "S" in val):
						strandInChrom = "+"
				if strandInChrom == "+": # adjust posInChrom if needed
					if not any("S" in sub for sub in two1):
						posInChrom = int(posInChrom) + sum([int(i[0]) for i in two1 if "M" in i or "D" in i]) - 1
					elif next(i for i, val in enumerate(two1, 0) if "M" in val) < next(i for i, val in enumerate(two1, 0) if "S" in val): ## check if M before S
						posInChrom = int(posInChrom) + sum([int(i[0]) for i in two1 if "M" in i or "D" in i]) - 1		

				# break 2
				outChrom = w[6].replace("=", inChrom)
				posOutChrom = int(w[7])
				strandOutChrom = "-" if "1" == flgBin[5] else "+" # strand from insert size to get orientation of the translocation
				strandOutChromSA = ""
				if w[10].startswith("SA:Z") and w[10].split(":")[2].split(",")[0] in chroms: # get from splits if any
					outChrom = w[10].split(":")[2].split(",")[0]
					posOutChrom = int(w[10].split(":")[2].split(",")[1])
					strandOutChromSA = w[10].split(":")[2].split(",")[2]
					split2 = re.findall(r'[A-Za-z]|[0-9]+', w[10].split(":")[2].split(",")[3])
					two2 = [split2[x:x+2] for x in range(0, len(split2),2)]
					if strandOutChromSA == "+" and next(i for i, val in enumerate(two2, 0) if "S" in val) < next(i for i, val in enumerate(two2, 0) if "M" in val):
						strandOutChromSA = "-"
					if strandOutChromSA == "-" and next(i for i, val in enumerate(two2, 0) if "M" in val) < next(i for i, val in enumerate(two2, 0) if "S" in val):
						strandOutChromSA = "+"
					if next(i for i, val in enumerate(two2, 0) if "M" in val) < next(i for i, val in enumerate(two2, 0) if "S" in val):  ## check if M before S
						posOutChrom = int(posOutChrom) + sum([int(i[0]) for i in two2 if "M" in i or "D" in i]) - 1
				if strandOutChromSA != "": strandOutChrom = strandOutChromSA
				
				# add PoN count
				if outChrom not in translocationsFiltered[inChrom]: continue
				for trans in translocationsFiltered[inChrom][outChrom]:
					if trans[0] == inChrom and int(trans[1])-200 <= posInChrom and int(trans[2])+200 >= posInChrom and trans[3] == strandInChrom and trans[4] == outChrom and int(trans[5])-1000 <= posOutChrom and int(trans[6])+1000 >= posOutChrom and trans[7] == strandOutChrom:
						trans[13] = trans[13]+1
						readNamesUsedInPoN.append(w[0])
		
		samfile.close()

	# 4. Prepare output, annotate RepeatMasker and GeneID, and return
	mask_expand = 20
	if genomeVersion == "hg19":
		RepeatMasker_dicti = pickle.load(gzip.open(inputsFolder+'/hg19/dicts/RepeatMasker_rmsk_hg19_dictionary.pkl.gz', 'rb'))
		if genesOncoIg == "all":
			GeneID_dicti = pickle.load(gzip.open(inputsFolder+'/hg19/dicts/GRCh37.p13_Biomart_HGNCsymbol.pkl.gz', 'rb'))
		else:
			GeneID_dicti = pickle.load(gzip.open(inputsFolder+'/hg19/dicts/GRCh37.p13_Biomart_HGNCsymbol_proteinCoding.pkl.gz', 'rb'))
	else:
		RepeatMasker_dicti = pickle.load(gzip.open(inputsFolder+'/hg38/dicts/RepeatMasker_rmsk_hg38_dictionary.pkl.gz', 'rb'))
		if genesOncoIg == "all":
			GeneID_dicti = pickle.load(gzip.open(inputsFolder+'/hg38/dicts/GRCh38.p14_Biomart_HGNCsymbol.pkl.gz', 'rb'))
		elif genesOncoIg == "protein_coding":
			GeneID_dicti = pickle.load(gzip.open(inputsFolder+'/hg38/dicts/GRCh38.p14_Biomart_HGNCsymbol_proteinCoding.pkl.gz', 'rb'))
		else: # protein_coding_canonical
			GeneID_dicti = pickle.load(gzip.open(inputsFolder+'/hg38/dicts/GRCh38.p14_Biomart_HGNCsymbol_proteinCodingCanonical.pkl.gz', 'rb'))

	translocationsList = []
	for key1 in translocationsFiltered:
		for key2 in translocationsFiltered[key1]:
			for item in translocationsFiltered[key1][key2]: 
				translocationsList.append(item)

	translocationsList = sorted(translocationsList, key=operator.itemgetter(8), reverse=True)
	translocationsALL = list()
	translocationsPASS = list()
	translocationsALL.append("\t".join(["Rearrangement", "Mechanism", "Score", "MQ", "Num_reads", "Read_types", "Depths_and_VAF", "Reads_in_normal", "Count_in_PoN", "Repeat_masker", "Chr_A", "Position_A", "Strand_A", "Chr_B", "Position_B", "Strand_B", "N_nucleotides", "Gene_ID", "Distance_to_gene"])+("" if reportReadNames == "no" else "\tRead_names"))

	for i in translocationsList:

		if i[0] == i[4]:
			if int((i[1] if i[3] == "-" else i[2])) < int((i[5] if i[7] == "-" else i[6])):
				idxA = 0
				idxB = 4
			else:
				idxA = 4
				idxB = 0
			
			chrA = i[idxA]
			strandA = i[idxA+3]
			positionA = i[idxA + (2 if strandA == "+" else 1)]
			chrB = i[idxB]
			strandB = i[idxB+3]
			positionB = i[idxB + (2 if strandB == "+" else 1)]
			
			if strandA == "+" and strandB == "-": 
				mechanism = "Deletion"
				traAnnot = "del("+chrA+":"+positionA+"-"+positionB+")"
			
			elif strandA == "-" and strandB == "+": 
				mechanism = "Gain"
				traAnnot = "gain("+chrA+":"+positionA+"-"+positionB+")"
			else: 
				mechanism = "Inversion"
				traAnnot = "inv("+chrA+":"+positionA+"-"+positionB+")"
		
		else:
			mechanism = "Translocation"
			
			traAnnot = ("t("+str(min([int(i[0].replace("chr", "").replace("X", "23").replace("Y", "24")), int(i[4].replace("chr", "").replace("X", "23").replace("Y", "24"))]))+";"+str(max([int(i[0].replace("chr", "").replace("X", "23").replace("Y", "24")), int(i[4].replace("chr", "").replace("X", "23").replace("Y", "24"))]))+")").replace("23", "X").replace("24", "Y")
			
			if i[0].replace("chr", "") == traAnnot.replace("t(", "").split(";")[0]: 
				idxA = 0
				idxB = 4
			else: 
				idxA = 4
				idxB = 0
			
			chrA = i[idxA]
			strandA = i[idxA+3]
			positionA = i[idxA + (2 if strandA == "+" else 1)]
			chrB = i[idxB]
			strandB = i[idxB+3]
			positionB = i[idxB + (2 if strandB == "+" else 1)]

		score = round( i[8] / tumorPurity, 1 )
		nNucleotidesFinal = "NA" if all(i == "NA" for i in i[9].split("-")) else Counter(i for i in i[9].split("-") if i != "NA").most_common(1)[0][0]
		nSplits = i[10].split("-").count("split")
		nPaired = i[10].split("-").count("paired")
		readTypeFinal = str(nSplits)+" split + "+str(nPaired)+" paired"
		quals = [int(q) for q in i[11].split("-")]
		mapQualReport = str(round(mean(quals),1))+" ("+str(min(quals))+"-"+str(max(quals))+")"
		readNamesReport = ",".join(set(i[12].split(",")))
		numReads = len(readNamesReport.split(","))
		scoreNormal = i[13]
		
		## RepeatMasker and GeneID:
		minDistance = genesOncoIgDistance
		repeatMasker = "none"
		if genomeVersion == "hg19":
			AllGenesBedToOpen = inputsFolder+"/hg19/dicts/AllRegionsAndGenes_hg19.bed"
		elif genomeVersion == "hg38":
			AllGenesBedToOpen = inputsFolder+"/hg38/dicts/AllRegionsAndGenes_hg38.bed"
		
		## Iterate breakA
		breakAisIG = "no"
		if geneToAnalyze != "tcr" and chrA == chrom+"14" and int(positionA) >= int(chrom14_IGH[0]) and int(positionA) <= int(chrom14_IGH[1]): geneID = "IGH"; locusID = "IGH"; breakAisIG = "IG"
		elif geneToAnalyze != "tcr" and chrA == chrom+"22" and int(positionA) >= int(chrom22_IGL[0]) and int(positionA) <= int(chrom22_IGL[1]): geneID = "IGL"; locusID = "IGL"; breakAisIG = "IG"
		elif geneToAnalyze != "tcr" and chrA == chrom+"2" and int(positionA) >= int(chrom2_IGK[0]) and int(positionA) <= int(chrom2_IGK[1]): geneID = "IGK"; locusID = "IGK"; breakAisIG = "IG"
		elif geneToAnalyze != "ig" and chrA == chrom+"14" and int(positionA) >= int(chrom14_TRA_TRD[0]) and int(positionA) <= int(chrom14_TRA_TRD[1]): geneID = "TRA_TRD"; locusID = "TRA_TRD"; breakAisIG = "TCR"
		elif geneToAnalyze != "ig" and chrA == chrom+"7" and int(positionA) >= int(chrom7_TRB[0]) and int(positionA) <= int(chrom7_TRB[1]): geneID = "TRB"; locusID = "TRB"; breakAisIG = "TCR"
		elif geneToAnalyze != "ig" and chrA == chrom+"7" and int(positionA) >= int(chrom7_TRG[0]) and int(positionA) <= int(chrom7_TRG[1]): geneID = "TRG"; locusID = "TRG"; breakAisIG = "TCR"

		if breakAisIG != "no":
			AllGenesBed = open(AllGenesBedToOpen, "r")
			for AllGenesBedLine in AllGenesBed:
				AllGenesBedList = AllGenesBedLine.rstrip("\n").split("\t")
				if chrA.replace("chr", "") == AllGenesBedList[0] and int(positionA) >= (int(AllGenesBedList[1])-5) and int(positionA) <= (int(AllGenesBedList[2])+5):
					geneID = AllGenesBedList[3]
					break
			AllGenesBed.close()
		else:
			gene = ""
			for element in GeneID_dicti[chrA.replace("chr","")]:
				if element[2] in customGenesOncoIg and ( ( int(positionA) >= int(element[0]) and int(positionA) <= int(element[1]) ) or abs(int(positionA) - int(element[0])) < customGenesOncoIgDistance or abs(int(positionA) - int(element[1])) < customGenesOncoIgDistance ):
					gene = element[2]
					if int(positionA) >= int(element[0]) and int(positionA) <= int(element[1]): 
						minDistance = 0
					else:
						minDistance = abs(int(positionA) - int(element[0])) if abs(int(positionA) - int(element[0])) < abs(int(positionA) - int(element[1])) else abs(int(positionA) - int(element[1]))
					break
				elif ( int(positionA) >= int(element[0]) and int(positionA) <= int(element[1]) ) or abs(int(positionA) - int(element[0])) < minDistance or abs(int(positionA) - int(element[1])) < minDistance:
					gene = element[2]
					if int(positionA) >= int(element[0]) and int(positionA) <= int(element[1]): 
						minDistance = 0
						if customGenesOncoIg == [""]: break
					else: 
						minDistance = abs(int(positionA) - int(element[0])) if abs(int(positionA) - int(element[0])) < abs(int(positionA) - int(element[1])) else abs(int(positionA) - int(element[1]))
			if gene.startswith("IGHV"): repeatMasker = "IGHV_pseudogene"
			if gene.startswith("IGHD"): repeatMasker = "IGHD_pseudogene"
			geneID = gene if gene != "" else "none"
			locusID = "NoIG"
			minDistance = minDistance if minDistance < genesOncoIgDistance else "NA"
			
			for element in RepeatMasker_dicti[chrA.replace("chr","")]:
				if int(positionA) >= int(element[0]) - mask_expand and int(positionA) <= int(element[1]) + mask_expand:
					repeatMasker = element[2]
					break
		
		## Iterate breakB
		breakBisIG = "no"
		if geneToAnalyze != "tcr" and chrB == chrom+"14" and int(positionB) >= int(chrom14_IGH[0]) and int(positionB) <= int(chrom14_IGH[1]): geneID = geneID+" - IGH"; locusID = locusID+" - IGH"; breakBisIG = "IG"
		elif geneToAnalyze != "tcr" and chrB == chrom+"22" and int(positionB) >= int(chrom22_IGL[0]) and int(positionB) <= int(chrom22_IGL[1]): geneID = geneID+" - IGL"; locusID = locusID+" - IGL"; breakBisIG = "IG"
		elif geneToAnalyze != "tcr" and chrB == chrom+"2" and int(positionB) >= int(chrom2_IGK[0]) and int(positionB) <= int(chrom2_IGK[1]): geneID = geneID+" - IGK"; locusID = locusID+" - IGK"; breakBisIG = "IG"
		elif geneToAnalyze != "ig" and chrB == chrom+"14" and int(positionB) >= int(chrom14_TRA_TRD[0]) and int(positionB) <= int(chrom14_TRA_TRD[1]): geneID = geneID+" - TRA_TRD"; locusID = locusID+" - TRA_TRD"; breakBisIG = "TCR"
		elif geneToAnalyze != "ig" and chrB == chrom+"7" and int(positionB) >= int(chrom7_TRB[0]) and int(positionB) <= int(chrom7_TRB[1]): geneID = geneID+" - TRB"; locusID = locusID+" - TRB"; breakBisIG = "TCR"
		elif geneToAnalyze != "ig" and chrB == chrom+"7" and int(positionB) >= int(chrom7_TRG[0]) and int(positionB) <= int(chrom7_TRG[1]): geneID = geneID+" - TRG"; locusID = locusID+" - TRG"; breakBisIG = "TCR"
		
		if breakAisIG == "no" and breakBisIG == "no": continue # Remove if none of the breaks are within an IG/TCR locus

		if breakBisIG != "no":
			AllGenesBed = open(AllGenesBedToOpen, "r")
			for AllGenesBedLine in AllGenesBed:
				AllGenesBedList = AllGenesBedLine.rstrip("\n").split("\t")
				if chrB.replace("chr", "") == AllGenesBedList[0] and int(positionB) >= (int(AllGenesBedList[1])-5) and int(positionB) <= (int(AllGenesBedList[2])+5):
					geneID = geneID.split(" - ")[0]+" - "+AllGenesBedList[3]
					break
			AllGenesBed.close()
		else:
			gene = ""
			for element in GeneID_dicti[chrB.replace("chr","")]:
				if element[2] in customGenesOncoIg and ( ( int(positionB) >= int(element[0]) and int(positionB) <= int(element[1]) ) or abs(int(positionB) - int(element[0])) < customGenesOncoIgDistance or abs(int(positionB) - int(element[1])) < customGenesOncoIgDistance ):
					gene = element[2]
					if int(positionB) >= int(element[0]) and int(positionB) <= int(element[1]):
						minDistance = 0
					else:
						minDistance = abs(int(positionB) - int(element[0])) if abs(int(positionB) - int(element[0])) < abs(int(positionB) - int(element[1])) else abs(int(positionB) - int(element[1]))
					break
				elif ( int(positionB) >= int(element[0]) and int(positionB) <= int(element[1]) ) or abs(int(positionB) - int(element[0])) < minDistance or abs(int(positionB) - int(element[1])) < minDistance:
					gene = element[2]
					if int(positionB) >= int(element[0]) and int(positionB) <= int(element[1]):
						minDistance = 0
						if customGenesOncoIg == [""]: break
					else:
						minDistance = abs(int(positionB) - int(element[0])) if abs(int(positionB) - int(element[0])) < abs(int(positionB) - int(element[1])) else abs(int(positionB) - int(element[1]))
			gene = gene if gene != "" else "none"
			if gene.startswith("IGHV"): repeatMasker = "IGHV_pseudogene"
			if gene.startswith("IGHD"): repeatMasker = "IGHD_pseudogene"
			geneID = geneID+" - "+gene
			locusID = locusID+" - "+"NoIG"
			minDistance = minDistance if minDistance < genesOncoIgDistance else "NA"
			
			for element in RepeatMasker_dicti[chrB.replace("chr","")]:
				if int(positionB) >= int(element[0]) - mask_expand and int(positionB) <= int(element[1]) + mask_expand:
					repeatMasker = element[2]
					break
		
		## Remove sv within the same IG/TCR locus
		if locusID  in ["IGH - IGH", "IGL - IGL", "IGK - IGK", "TRA_TRD - TRA_TRD", "TRB - TRB", "TRG - TRG"]: continue
		
		## Calculate VAF of the IG/TCR breakpoint
		if breakAisIG != "no": 
			if strandA == "+": region = chrA+":"+str(int(positionA)-19)+"-"+str(int(positionA))
			else: region = chrA+":"+str(int(positionA))+"-"+str(int(positionA)+19)
		else:
			if strandB == "+": region = chrB+":"+str(int(positionB)-19)+"-"+str(int(positionB))
			else: region = chrB+":"+str(int(positionB))+"-"+str(int(positionB)+19)

		# depth
		comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapqOnco+" -r "+region+" "+bamT+" 2>> "+errLogMpileup+" | cut -f 4"
		process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		stdout, stderr = process.communicate()
		depth = int(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/20, 0))

		# altDepth
		readsToFile = readNamesReport.split(",")
		readFile = wkDir+"/tmp/"+bamT.split("/")[-1].replace(".bam", "_readnames.txt")
		readBam = wkDir+"/tmp/"+bamT.split("/")[-1].replace(".bam", "_readnames.bam")
		READFILE = open(readFile, "w")
		READFILE.write("\n".join(readsToFile))
		READFILE.close()
		comms = pathToSamtools+"samtools view -@ "+threadsForSamtools+" -b -h -N "+readFile+" -o "+readBam+" "+bamT
		subprocess.call(comms, shell=True)
		comms = pathToSamtools+"samtools index "+readBam
		subprocess.call(comms, shell=True)
		comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapqOnco+" -r "+region+" "+readBam+" 2>> "+errLogMpileup+" | cut -f 4"
		process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		stdout, stderr = process.communicate()
		altDepth = int(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/20, 0))
		if altDepth > len(readsToFile): altDepth = len(readsToFile)
		if altDepth > 0:
			vaf = round(altDepth/depth*100, 2)
			vafAdj = round((altDepth/depth*100)/tumorPurity, 2)
			if vafAdj >= 100: vafAdj = 100
		else:
			vaf = 0
			vafAdj = 0
		vafString = str(altDepth)+"/"+str(depth)+" ("+str(vaf)+"% ["+str(vafAdj)+"%])"

		## PoN:
		ponCount = 0
		
		if breakAisIG != "no": 
			igLocus = locusID.split(" - ")[0]
			igOrder = 0
		else:
			igLocus = locusID.split(" - ")[1]
			igOrder = 1
		
		if genomeVersion == "hg19": PoN = open(inputsFolder+'/hg19/PoN/hg19_PoN.tsv', 'r')
		else: PoN = open(inputsFolder+'/hg38/PoN/hg38_PoN.tsv', 'r')
		
		for ponLine in PoN:
			if ponLine.startswith("Rearrangement"): continue
			ponList = ponLine.rstrip("\n").split("\t")
			if igLocus in ponList[15].split("::")[igOrder]:
				if igOrder == 0:
					if chrB.replace("chr", "") == ponList[11].replace("chr", "") and int(positionB) >= int(ponList[12])-1000 and int(positionB) <= int(ponList[12])+1000 and strandB == ponList[13]:
						ponCount += 1
				else:
					if chrA.replace("chr", "") == ponList[8].replace("chr", "") and int(positionA) >= int(ponList[9])-1000 and int(positionA) <= int(ponList[9])+1000 and strandA == ponList[10]:
						ponCount += 1
		
		PoN.close()
		
		## Convert geneID to ::
		geneID = geneID.replace(" - ", "::")		

		## Return all
		if vafAdj > 0:
			translocationsALL.append("\t".join([traAnnot, mechanism, str(score), mapQualReport, str(numReads), readTypeFinal, vafString, str(scoreNormal), str(ponCount), repeatMasker, chrA, positionA, strandA, chrB, positionB, strandB, nNucleotidesFinal, geneID, str(minDistance)])+("" if reportReadNames == "no" else "\t"+readNamesReport))
		
		## Return pass
		if score >= mntoncoPass and vafAdj >= vafOnco*100 and ( scoreNormal == "NA" or scoreNormal <= mnnonco ) and ponCount <= mncPoN:
			if mechanism == "Translocation": traAnnot = traAnnot+" ["+chrA+":"+positionA+":"+strandA+";"+chrB+":"+positionB+":"+strandB+"] ["+nNucleotidesFinal+"] ["+geneID+"] ["+str(vafAdj)+"%]"
			else: traAnnot = traAnnot+" ["+strandA+"/"+strandB+"] ["+nNucleotidesFinal+"] ["+geneID+"] ["+str(vafAdj)+"%]"
			translocationsPASS.append("\t".join(["Oncogenic "+("IG" if geneToAnalyze == "ig" else "TCR" if geneToAnalyze == "tcr" else "IG/TCR")+" rearrangement", traAnnot, mechanism, str(score)+" ("+str(scoreNormal)+") ["+str(ponCount)+"] ["+repeatMasker+"]", mapQualReport, str(numReads)]+["NA"]*6)+("" if reportReadNames == "no" else "\t"+readNamesReport))
		
	return(translocationsALL, translocationsPASS)

def getPurity(wkDir, seq, chrom, genomeVersion, inputsFolder, chrAnnot, filterOutputFile, listGenes, estimatePurityCoverage, bamT, bamN, seqDepth, pathToSamtools, mapq, scoreCutoffPurity, plotPurityCoverage, reportReadNames, errLogMpileup):
	
	if genomeVersion == "hg19":	
		bedFile = inputsFolder+"/hg19/"+chrAnnot+"/wgEncodeGencodeBasicV19_hg19_JgenesForPurity.bed"
		IGKKdePos = 89132285
		IGKRSSPos = 89159112
	else:
		bedFile = inputsFolder+"/hg38/"+chrAnnot+"/GencodeV40_hg38_JgenesForPurity.bed"
		IGKKdePos = 88832765
		IGKRSSPos = 88859600
	
	CovReductionAll = list()
	CovReductionSelected = list()
	puritySampleList = list()
	purityCutoff = 0.15 if bamN is not None and seqDepth == "high" else 0.20 if seqDepth == "high" or bamN is not None else 0.25
	covReductionCutoff = 0.15 if bamN is not None and seqDepth == "high" else 0.25 if seqDepth == "high" or bamN is not None else 0.3
	windowDeletedBreak = 50 if seq == "wgs" else 10
	listRegionsPlot = [] # list of list with regions to plot (chr, start, end, GENE, gene-name)
	
	# Iterate each locus analyzed
	for GENE in listGenes:

		if GENE == "CSR" or GENE == "TRD": continue
		if GENE == "TRA": 
			GENEtoStore = "TRA_TRD"
			genesToMatch = tuple(["TRAJ", "TRDJ"])
		else: 
			GENEtoStore = GENE
			genesToMatch = tuple([GENE+"J"])

		countGeneRearrangement = 0
		puritySpecificRearrangement = 0
		purityGeneRearrangement = list()
		CovReductionSelectedRearrangement = list()
		purityGene = 0
		locusCompleted = "no"
		
		# Create list of IGHJ genes found rearranged in GENE
		JgenesRearrangedList = []
		SUMM = open(filterOutputFile, "r")
		for sLine in SUMM:
			sList = sLine.rstrip("\n").split("\t")
			if any([True if i in sList[1] else False for i in genesToMatch]):
				score = float(sList[3].split(" ")[0])
				if score < scoreCutoffPurity: continue
				if sList[0].startswith("Oncogenic"):
					if sList[2] != "Translocation": continue
					allGenes = sList[1].split(" ")[3].replace("[", "").replace("]", "").split("::")
					for g in allGenes:
						if g.startswith(genesToMatch):
							JgenesRearrangedList.append([g, "geneOnly", "geneOnly", score])
							break
				else:
					if sList[2] != "Deletion": continue
					if sList[9] == "Partial rearrangement":
						g = sList[1].split(" - ")[0]
						JgenesRearrangedList.append([g, "geneOnly", "geneOnly", score])
					else:
						jGenes = ",".join([j.split("*")[0] for j in sList[1].split(" - ")[0].split(",")])
						locus = jGenes.split("J")[0]
						seq = sList[11]
						CHAIN = open(filterOutputFile.replace("filtered.tsv", locus+".tsv"), "r")
						for cLine in CHAIN:
							cList = cLine.rstrip("\n").split("\t")
							if seq == cList[15]:
								JgenesRearrangedList.append([jGenes, cList[4], cList[5], score])
								break
						CHAIN.close()
		SUMM.close()
		if len(JgenesRearrangedList) > 2: # limit to two J genes per GENE
			JgenesRearrangedList.sort(key=lambda p: p[3], reverse=True) # order based on score
			JgenesRearrangedList = JgenesRearrangedList[:2]
		
		# Open bed file and iterate over J genes
		jCount = 0
		B = open(bedFile, "r")
		for bLine in B:
			bList = bLine.rstrip("\n").split("\t")
			if bList[3].startswith(genesToMatch):
				jCount += 1
				chr = bList[0]
				posStart = bList[1]
				posEnd = bList[2]
				gene = bList[3]
				strand = bList[5]

				# region deleted after break
				posBreak = int(bList[2]) if strand == "-" else int(bList[1])
				if strand == "-":
					regionDeleted_1 = chr+":"+str(posBreak+10)+"-"+str(posBreak+10+windowDeletedBreak-1)
				else:
					regionDeleted_1 = chr+":"+str(posBreak-10-windowDeletedBreak+1)+"-"+str(posBreak-10)
				
				# region deleted J gene
				adjustPos = 0 if gene.endswith("_DEL") else 15
				regionDeletedGene = chr+":"+posStart+"-"+str(int(posEnd) - adjustPos) if strand == "-" else chr+":"+str(int(posStart) + adjustPos)+"-"+posEnd
				windowDeletedGene = int(posEnd) - int(posStart) - adjustPos + 1
				
				# region to store for ploting
				if not gene.endswith("_DEL"):
					if strand == "-":
						pos1 = int(bList[1])
						pos2 = int(bList[2])+10+windowDeletedBreak-1
					else:
						pos1 = int(bList[1])-10-windowDeletedBreak+1
						pos2 = int(bList[2])
					locusForPlot = "TRA_TRD" if GENE in ["TRA", "TRD"] else GENE
					listRegionsPlot.append([chr, pos1, pos2, locusForPlot, gene])

				# store or calculate coverage
				if jCount == 1:
					geneRearranged = gene
					posStartGeneRearranged = posStart
					posEndGeneRearranged = posEnd
					regionNormal = regionDeletedGene
					windowNormal = windowDeletedGene
					regionDeletedBreak = regionDeleted_1					

				else:
					# Get mean depths tumor bam
					comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -r "+regionNormal+" "+bamT+" 2>> "+errLogMpileup+" | cut -f 4"
					process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
					stdout, stderr = process.communicate()
					depthNormal = float(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/windowNormal, 0))
					if ( seqDepth == "high" and depthNormal <= 50 ) or ( seqDepth == "int" and depthNormal <= 15 ) or ( seqDepth == "low" and depthNormal <= 7 ): flagCov = "LowCoverage"
					else: flagCov = "PASS"

					comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -r "+regionDeletedBreak+" "+bamT+" 2>> "+errLogMpileup+" | cut -f 4"
					process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
					stdout, stderr = process.communicate()
					depthDeletedBreak = float(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/windowDeletedBreak, 0))
					
					comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -r "+regionDeletedGene+" "+bamT+" 2>> "+errLogMpileup+" | cut -f 4"
					process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
					stdout, stderr = process.communicate()
					depthDeletedGene = float(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/windowDeletedGene, 0))

					# Get mean depths normal bam
					factorDeletedBreak = 1
					factorDeletedGene = 1
					if bamN is not None:
						comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -r "+regionNormal+" "+bamN+" 2>> "+errLogMpileup+" | cut -f 4"
						process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
						stdout, stderr = process.communicate()
						depthNormalNormal = float(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/windowNormal, 0))
						
						comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -r "+regionDeletedBreak+" "+bamN+" 2>> "+errLogMpileup+" | cut -f 4"
						process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
						stdout, stderr = process.communicate()
						depthDeletedBreakNormal = float(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/windowDeletedBreak, 0))
						
						comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -r "+regionDeletedGene+" "+bamN+" 2>> "+errLogMpileup+" | cut -f 4"
						process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
						stdout, stderr = process.communicate()
						depthDeletedGeneNormal = float(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/windowDeletedGene, 0))

						if depthNormalNormal > 10:
							if depthDeletedBreakNormal > 5:
								factorDeletedBreak = depthDeletedBreakNormal / depthNormalNormal
							if depthDeletedGeneNormal > 5: 
								factorDeletedGene = depthDeletedGeneNormal / depthNormalNormal
					
					# Adjust depthDeleted based on normal-BAM-derived factorDeleted and calculate covReduction
					depthDeletedBreak = round(depthDeletedBreak / factorDeletedBreak, 0)
					depthDeletedGene = round(depthDeletedGene / factorDeletedGene, 0)
					if depthNormal > 0:
						covReductionBreak = round(1 - depthDeletedBreak / depthNormal, 3)
						covReductionGene = round(1 - depthDeletedGene / depthNormal, 3)
						if covReductionBreak < 0: covReductionBreak = 0
						if covReductionGene < 0: covReductionGene = 0
					else:
						covReductionBreak = 0
						covReductionGene = 0

					# Check if IGHJ gene found in a rearranged allele to estimate purity
					countGeneSpecificRearrangement = 0
					countGeneSpecificRearrangementDetection = "no"
					for jgrl in JgenesRearrangedList:
						if jgrl[1] == "geneOnly":
							if geneRearranged == jgrl[0]:
								countGeneSpecificRearrangement += 1
						elif (strand == "-" and int(posStartGeneRearranged) == int(jgrl[1])) or (strand == "+" and int(posEndGeneRearranged) == int(jgrl[2])):
							countGeneSpecificRearrangement += 1
					if countGeneSpecificRearrangement > 0:
						countGeneSpecificRearrangementDetection = "rearrangement"
						countGeneRearrangement = countGeneRearrangement + countGeneSpecificRearrangement
					
					# Potentially missed rearrangements, based on coverage if specified
					if ( estimatePurityCoverage == "yes" or (estimatePurityCoverage == "igh" and GENE == "IGH") ) and flagCov == "PASS" and countGeneSpecificRearrangement == 0 and covReductionBreak > covReductionCutoff and covReductionGene > covReductionCutoff: 
						countGeneSpecificRearrangement = 1
						countGeneSpecificRearrangementDetection = "coverage"
						countGeneRearrangement = countGeneRearrangement + countGeneSpecificRearrangement

					# Select and summarize covReduction
					if ( covReductionBreak > covReductionCutoff and covReductionGene > covReductionCutoff and flagCov == "PASS") or countGeneSpecificRearrangement > 0:
						if covReductionBreak >= 0.9 and covReductionGene <= 0.7:
							covReduction = round(covReductionBreak, 3)
						elif covReductionGene >= 0.9 and covReductionBreak <= 0.7:
							covReduction = round(covReductionGene, 3)
						elif covReductionBreak >= 0.5 and covReductionGene <= 0.1:
							covReduction = round(covReductionBreak, 3)
						else:
							covReduction = round(mean([covReductionBreak, covReductionGene]), 3)
					else:
						covReduction = min([covReductionBreak, covReductionGene])

					# Adjust countGeneRearrangement and countGeneSpecificRearrangement if covReduction of biallelic rearrangement
					if covReduction >= 0.80 and countGeneRearrangement < 2 and flagCov == "PASS":
						while countGeneRearrangement < 2:
							countGeneRearrangement += 1
							countGeneSpecificRearrangement += 1

					# Calculate purity
					if countGeneRearrangement < 2: purity = covReduction * 2
					else: purity = covReduction
					if purity >= 1: purity = 1
					if purity < 0: purity = 0

					# Adjust flag if found rearranged but with no covReduction
					if countGeneSpecificRearrangement > 0 and purity < 0.08: flagCov = "NoDepthReduction"

					# Keep all info
					listToStore = [GENEtoStore, geneRearranged, regionNormal, regionDeletedBreak, regionDeletedGene, str(int(depthNormal)), str(int(depthDeletedBreak)), str(int(depthDeletedGene)),str(covReductionBreak), str(covReductionGene), str(covReduction), str(countGeneSpecificRearrangement), str(countGeneRearrangement), str(purity), flagCov]
					CovReductionAll.append("\t".join(listToStore))

					# Keep highest purities per gene (ie locus)
					if flagCov == "PASS" and (locusCompleted == "no" or countGeneSpecificRearrangementDetection == "rearrangement"):
						if countGeneSpecificRearrangement > 0:
							purityGeneRearrangement.append(purity)
							CovReductionSelectedRearrangement.append("\t".join(listToStore))

						elif ( estimatePurityCoverage == "yes" or (estimatePurityCoverage == "igh" and GENE == "IGH") ) and ( (seqDepth == "high" and purity > purityGene and purity > purityCutoff) or (seqDepth != "high" and purity > purityGene and purity > purityCutoff and covReduction > covReductionCutoff) ):				
							purityGene = purity
							CovReductionGeneInfo = "\t".join(listToStore)

					# Store regionDeleted and current gene for next iteration
					geneRearranged = gene
					posStartGeneRearranged = posStart
					posEndGeneRearranged = posEnd
					regionNormal = regionDeletedGene
					windowNormal = windowDeletedGene
					regionDeletedBreak = regionDeleted_1

					# check locusCompleted
					if countGeneRearrangement >= 2:
						locusCompleted = "yes"

		B.close()

		if purityGeneRearrangement:
			cutoff = sorted(purityGeneRearrangement)[-2] if len(purityGeneRearrangement) >= 2 else max(purityGeneRearrangement)
			puritySpecificRearrangement = max(purityGeneRearrangement)
			if seqDepth == "high" and puritySpecificRearrangement > 0.9 and cutoff < 0.45: cutoff = puritySpecificRearrangement
			if seqDepth == "int" and puritySpecificRearrangement > 0.85 and cutoff < 0.15: cutoff = puritySpecificRearrangement
			for indx in range(0, len(purityGeneRearrangement), 1):
				if purityGeneRearrangement[indx] >= cutoff:
					CovReductionSelected.append(CovReductionSelectedRearrangement[indx])
					puritySampleList.append(purityGeneRearrangement[indx])
		
		if locusCompleted == "no":
			if purityGene > puritySpecificRearrangement:
				CovReductionSelected.append(CovReductionGeneInfo)
				puritySampleList.append(purityGene)
	
	# Check IGKKde/IGKRSS deletions (if IGK analyzed and if IGKKde/IGKRSS deleteions identified)
	if "IGK" in listGenes:
		
		kdes = 0
		flag_kdes = "PASS"
		rsss = 0
		flag_rsss = "PASS"
		kde_rss_s = 0
		SUMM = open(filterOutputFile, "r")
		for sLine in SUMM:
			sList = sLine.rstrip("\n").split("\t")
			if "IGK" == sList[0] and "IGKKde" in sList[1]:
				if float(sList[3].split(" ")[0]) < scoreCutoffPurity: continue
				if sList[2] == "Deletion": kdes += 1
				else: flag_kdes = "PotentialInversion"
			if "IGK" == sList[0] and "IGKRSS" in sList[1]:
				if float(sList[3].split(" ")[0]) < scoreCutoffPurity: continue
				if sList[2] == "Deletion": rsss += 1
				else: flag_rsss = "PotentialInversion"
			if "IGK" == sList[0] and "IGKKde" in sList[1] and "IGKRSS" in sList[1]:
				kde_rss_s += 1
		SUMM.close()
				
		if kdes == 0 or flag_kdes == "PotentialInversion":
			flagToPrint = "NoRearranged" if kdes == 0 else "PotentialInversion"
			CovReductionGeneInfo = "".join(["IGKKde", "\tNA"*13, "\t"+flagToPrint])
			CovReductionAll.append(CovReductionGeneInfo)
		
		if rsss == 0 or flag_rsss == "PotentialInversion":
			flagToPrint = "NoRearranged" if rsss == 0 else "PotentialInversion"
			CovReductionGeneInfo = "".join(["IGKRSS", "\tNA"*13, "\t"+flagToPrint])
			CovReductionAll.append(CovReductionGeneInfo)
		
		regionsToIterate = ["IGKKde", "IGKRSS"] if kdes > 0 and rsss > 0 and flag_kdes == "PASS" and flag_rsss == "PASS" else ["IGKKde"] if kdes > 0 and flag_kdes == "PASS" else ["IGKRSS"] if rsss > 0 and flag_rsss == "PASS" else list()
		
		if regionsToIterate:
			for igkRegion in regionsToIterate:

				if igkRegion == "IGKKde": 
					igkPositionRearrangements = kdes
					regionNormal = chrom+"2:"+str(IGKKdePos - 15 - windowDeletedBreak)+"-"+str(IGKKdePos - 15 - 1)
					regionDeleted = chrom+"2:"+str(IGKKdePos + 15 + 1)+"-"+str(IGKKdePos + 15 + windowDeletedBreak)
					listRegionsPlot.append([chrom+"2", IGKKdePos - 15 - windowDeletedBreak, IGKKdePos + 15 + windowDeletedBreak, "IGKKde", "IGKKde"])
				elif igkRegion == "IGKRSS":
					igkPositionRearrangements = rsss
					regionNormal = chrom+"2:"+str(IGKRSSPos + 15 + 1)+"-"+str(IGKRSSPos + 15 + windowDeletedBreak)
					regionDeleted = chrom+"2:"+str(IGKRSSPos - 15 - windowDeletedBreak)+"-"+str(IGKRSSPos - 15 - 1)
					listRegionsPlot.append([chrom+"2", IGKRSSPos - 15 - windowDeletedBreak, IGKRSSPos + 15 + windowDeletedBreak, "IGKRSS", "IGKRSS"])

				# Get mean depths tumor windowDeletedBreak
				comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -r "+regionNormal+" "+bamT+" 2>> "+errLogMpileup+" | cut -f 4"
				process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
				stdout, stderr = process.communicate()
				depthNormal = float(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/windowDeletedBreak, 0))
				
				if ( seqDepth == "high" and depthNormal <= 50 ) or ( seqDepth == "int" and depthNormal <= 15 ) or ( seqDepth == "low" and depthNormal <= 7 ): 
					CovReductionGeneInfo = "".join([igkRegion, "\tNA"*13, "LowCoverage"])
					CovReductionAll.append(CovReductionGeneInfo)
				
				else:
					flagCov = "PASS"
					comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -r "+regionDeleted+" "+bamT+" 2>> "+errLogMpileup+" | cut -f 4"
					process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
					stdout, stderr = process.communicate()
					depthDeleted = float(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/windowDeletedBreak, 0))

					# Get mean depths normal bam
					factorDeleted = 1
					if bamN is not None:
						comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -r "+regionNormal+" "+bamN+" 2>> "+errLogMpileup+" | cut -f 4"
						process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
						stdout, stderr = process.communicate()
						depthNormalNormal = float(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/windowDeletedBreak, 0))

						comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -r "+regionDeleted+" "+bamN+" 2>> "+errLogMpileup+" | cut -f 4"
						process = subprocess.Popen(comms, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
						stdout, stderr = process.communicate()
						depthDeletedNormal = float(round(float(sum([int(cov) for cov in stdout.decode("utf-8").split("\n") if cov != ""]))/windowDeletedBreak, 0))
						
						if depthNormalNormal > 10: factorDeleted = depthDeletedNormal / depthNormalNormal
				
					# Calculate covReduction
					depthDeleted = round(depthDeleted / factorDeleted, 0)
					covReduction = round(1 - depthDeleted / depthNormal, 3)
					if covReduction < 0: covReduction = 0

					# Adjust kdes/rsss if covReduction of biallelic rearrangement
					if covReduction > 0.80 and igkPositionRearrangements < 2: igkPositionRearrangements = 2

					# Calculate purity
					if igkPositionRearrangements < 2: purity = covReduction * 2
					else: purity = covReduction
					if purity >= 1: purity = 1
					if purity < 0: purity = 0

					# Adjust, if necessary, purity if only one IGKKde-IGKRSS identified but it looks potentially biallelic considering also all other purities found
					if purity > 0.80 and kdes == 1 and rsss == 1 and kde_rss_s == 1 and all(purs < purity*0.7 for purs in puritySampleList): purity = round(purity/2, 3)

					# Adjust flag if found rearranged but with no covReduction
					if purity < 0.08: flagCov = "NoDepthReduction"
					
					# Store info
					CovReductionGeneInfo = "\t".join([igkRegion, igkRegion, regionNormal, regionDeleted, "NA", str(int(depthNormal)), str(int(depthDeleted)), "NA", str(covReduction), "NA", str(covReduction), str(igkPositionRearrangements), str(igkPositionRearrangements), str(purity), flagCov])
					if flagCov == "PASS":
						CovReductionSelected.append(CovReductionGeneInfo)
						puritySampleList.append(purity)
					else:
						CovReductionAll.append(CovReductionGeneInfo)

	# Print output...
	if len(puritySampleList) == 0:
		medPurity = "<"+str(purityCutoff)
		puritySample = "<"+str(purityCutoff)+" (NA)"
	else: 
		puritySampleList.sort()
		medPurity = round(median(puritySampleList), 3)
		if medPurity < purityCutoff: medPurity = "<"+str(purityCutoff)
		puritySample = str(medPurity)+" ("+";".join([str(p) for p in puritySampleList])+")"
	
	## ...filtered
	with open(filterOutputFile, 'a') as file: file.write("Purity\t"+str(puritySample)+"\t"+"\t".join(["NA"]*(10 if reportReadNames == "no" else 11))+"\n")
	
	## ...purity file
	locusToPlot = []
	locusCounted = []
	genesCounted = [] 
	O = open(filterOutputFile.replace("_filtered.tsv", "_purity.tsv"), "w")
	O.write("Analysis\tJ_gene\tNormal_region\tDeleted_region_break\tDeleted_region_next_gene\tDepth_normal_region\tDepth_deleted_region_break\tDepth_deleted_region_next_gene\tDepth_reduction_break\tDepth_reduction_next_gene\tDepth_reduction\tNumber_of_times_J_gene_specific_rearrangments\tNumber_of_times_J_gene_rearranged\tPurity\tFlag\tSelected\n")
	for eToPrint in CovReductionSelected:
		O.write(eToPrint+"\tYes\n")
		locusToPlot.append(eToPrint.split("\t")[0])
		locusCounted.append(eToPrint.split("\t")[0])
		genesCounted.extend([eToPrint.split("\t")[1]] * (int(eToPrint.split("\t")[12]) if int(eToPrint.split("\t")[12]) > 0 else 1))
	for eToPrint in CovReductionAll:
		if eToPrint not in CovReductionSelected:
			O.write(eToPrint+"\tNo\n")
			if eToPrint.split("\t")[-1] != "LowCoverage":
				locusToPlot.append(eToPrint.split("\t")[0])
	O.close()
	locusToPlot = list(set(locusToPlot))

	## ...plot coverage along J genes
	if plotPurityCoverage == "yes":
		
		# prepare DataFrame and bed file for mpileup
		chroms = []
		positions = []
		locus = []
		gene = []
		bedForMpileup = wkDir+"/tmp/bedFileForMpileup.bed"
		BEDFORMPILEUP = open(bedForMpileup, "w")
		for regionsForMpileup in listRegionsPlot:
			if regionsForMpileup[3] in locusToPlot:
				# for bed file
				BEDFORMPILEUP.write(str(regionsForMpileup[0])+"\t"+str(int(regionsForMpileup[1])-1)+"\t"+str(regionsForMpileup[2])+"\n")
				# for dataframe
				lengthRegion = int(regionsForMpileup[2])-int(regionsForMpileup[1])+1
				chroms.extend([str(regionsForMpileup[0])] * lengthRegion)
				positions.extend(list(range(int(regionsForMpileup[1]), int(regionsForMpileup[2])+1)))
				locus.extend([str(regionsForMpileup[3])] * lengthRegion)
				gene.extend([str(regionsForMpileup[4])] * lengthRegion)
		BEDFORMPILEUP.close()
		covDF = pd.DataFrame({ "chrom" : chroms, "position" : positions, "locus" : locus, "gene" : gene })

		# get coverage tumor
		covCountsTumor = wkDir+"/tmp/covCountsTumor.tsv"
		comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -l "+bedForMpileup+" "+bamT+" > "+covCountsTumor+" 2>> "+errLogMpileup
		subprocess.call(comms, shell=True)
		COV = open(covCountsTumor, "r")
		chroms = []
		positions = []
		cov = []
		for covLine in COV:
			covList = covLine.rstrip("\n").split("\t")
			chroms.append(str(covList[0]))
			positions.append(int(covList[1]))
			cov.append(int(covList[3]))
		COV.close()
		tumorCV = pd.DataFrame({"chrom" : chroms, "position" : positions, "tumorCov" : cov})
		covDF = pd.merge(covDF, tumorCV, on=['chrom', 'position'], how='left')
		covDF["tumorCov"] = covDF["tumorCov"].fillna(0)

		# get coverage normal
		if bamN is not None:
			covCountsNormal = wkDir+"/tmp/covCountsNormal.tsv"
			comms = pathToSamtools+"samtools mpileup -d 0 -a -A -B -q "+mapq+" -l "+bedForMpileup+" "+bamN+" > "+covCountsNormal+" 2>> "+errLogMpileup
			subprocess.call(comms, shell=True)
			COV = open(covCountsNormal, "r")
			chroms = []
			positions = []
			cov = []
			for covLine in COV:
				covList = covLine.rstrip("\n").split("\t")
				chroms.append(str(covList[0]))
				positions.append(int(covList[1]))
				cov.append(int(covList[3]))
			COV.close()
			normalCV = pd.DataFrame({"chrom" : chroms, "position" : positions, "normalCov" : cov})
			covDF = pd.merge(covDF, normalCV, on=['chrom', 'position'], how='left')
			covDF["normalCov"] = covDF["normalCov"].fillna(1)

		# normalize coverages
		covDF["tumorCovNormalized"] = covDF["tumorCov"] / covDF["tumorCov"].mean()
		if bamN is not None:
			covDF["normalCovNormalized"] = covDF["normalCov"] / covDF["normalCov"].mean()
			covDF["tumorCovNormalized"] = covDF["tumorCovNormalized"] / covDF["normalCovNormalized"]
		
		# re-name genes to display number of breaks counted
		for gene in covDF["gene"].unique():
			newGeneLabel = gene+" (n="+str(genesCounted.count(gene))+")"
			covDF.loc[covDF["gene"] == gene, "gene"] = newGeneLabel
		
		# plot
		covDF["gene"] = covDF["gene"].astype(str)
		covDF["locus"] = covDF["locus"].astype(str)
		with PdfPages(bamT.replace(".bam", "_output_purity_RawTumorCoverage.pdf")) as pdf:
			for locus in covDF["locus"].unique():
				covLocus = covDF[covDF["locus"] == locus]
				covLocus = covLocus.sort_values(by="position", ascending=True)
				nCol = covLocus["gene"].nunique() if covLocus["gene"].nunique() < 13 else 13
				nRow = int(covLocus["gene"].nunique()/nCol)
				g = sns.FacetGrid(covLocus, col="gene", sharex=False, margin_titles=True, col_wrap=nCol)
				g.map_dataframe(sns.lineplot, x="position", y="tumorCov", color="#32B39E")
				for ax in g.axes.flat: ax.set_ylim(bottom=0)
				g.set_axis_labels("Position", "Raw tumor coverage")
				g.set_titles(col_template="{col_name}")
				g.figure.subplots_adjust(top=0.85 if nRow == 1 else 0.95)
				g.figure.suptitle("Purity: "+str(medPurity)+" | Locus: "+locus+" (n= "+str(locusCounted.count(locus))+")", fontsize=11)
				g.figure.set_size_inches(3 * nCol, 4 * nRow)
				pdf.savefig(g.figure)
				plt.close()
		
		with PdfPages(bamT.replace(".bam", "_output_purity_NormalizedTumorCoverage.pdf")) as pdf:
			for locus in covDF["locus"].unique():
				covLocus = covDF[covDF["locus"] == locus]
				covLocus = covLocus.sort_values(by="position", ascending=True)
				nCol = covLocus["gene"].nunique() if covLocus["gene"].nunique() < 13 else 13
				nRow = int(covLocus["gene"].nunique()/nCol)
				g = sns.FacetGrid(covLocus, col="gene", sharex=False, margin_titles=True, col_wrap=nCol)
				g.map_dataframe(sns.lineplot, x="position", y="tumorCovNormalized", color="#32B39E")
				for ax in g.axes.flat: ax.set_ylim(bottom=0)
				g.set_axis_labels("Position", "Norm. tumor coverage")
				g.set_titles(col_template="{col_name}")
				g.figure.subplots_adjust(top=0.85 if nRow == 1 else 0.95)
				g.figure.suptitle("Purity: "+str(medPurity)+" | Locus: "+locus+" (n= "+str(locusCounted.count(locus))+")", fontsize=11)
				g.figure.set_size_inches(3 * nCol, 4 * nRow)
				pdf.savefig(g.figure)
				plt.close()
