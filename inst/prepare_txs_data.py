#!/usr/bin/env python

#	prepare_txs_data.py

#	Prepare training and test data at gene level for MaLTE.

#	Copyright (C) 2013  Paul K. Korir

#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.

#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.

#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
import sys
import tempfile
import argparse
import time
import os
import os.path
import gzip
import re
import multiprocessing

def get_data( f_fn, g_fn, ma_L, platform="ma", column=1 ):
	"""
	create a dictionary that maps feature IDs to expression values
	"""
	unwanted = [ '#' ]
	
	if platform == "ma":
		index = 2
		if re.search( r"gz$", f_fn ):
			f = gzip.open( f_fn )
		elif re.search( r"txt$", f_fn ):
			f = open( f_fn )
	elif platform == "hts":
		index = 1
		if re.search( r"gz$", f_fn ):
			f = gzip.open( f_fn )
		elif re.search( r"txt$", f_fn ):
			f = open( f_fn )
	
	g = open( g_fn, 'w+r' )
	
	names_L = list()
	genes_L = list()
	
	row_count = 0
	for row in f:
		if row[0] in unwanted: continue
		L = row.strip().split( '\t' )
		
		# make and print the header		
		if row[0] == 'p' or row[0] == 't' and row_count == 0:
			names_L = L[index:]
			if platform == "ma":
				print >> g, "\t".join( [ L[0], L[1] ] + ma_L )
			elif platform == "hts":
				if len( ma_L ) == 0:
					print >> g, "tx_id\nsamples"
				else:
					print >> g, "\t".join( [ L[0] ] + ma_L )
				row_count += 1
			continue
		elif row[0] != 'p' and row[0] != 't' and row_count == 0:
			print >> sys.stderr, "Error: missing header in HTS/microarray data file. Please amend and retry."
			sys.exit( 1 )

		expr_D = dict( zip( names_L, L[index:] ) )
		
		# test to see that the sample names of interest (ma_L) are contained in
		# the sample names in the header (names_L)
		try:
			assert len( set( names_L ).intersection( set( ma_L ))) == len( ma_L )
		except:
			print >> sys.stderr, "\nError: possible sample name mismatch. Please ensure that all sample names in '%s' file header are present in 'samples.txt'." % platform
			sys.exit( 1 )
		
		# make and print the data
		if platform == "ma":
			print >> g, "\t".join( [ L[0], L[column] ] + [ expr_D[k] for k in ma_L ] )
		elif platform == "hts":
			if len( ma_L ) == 0:
				print >> g, "%s\tNA" % L[0]
			else:
				print >> g, "\t".join( [ L[0] ] + [ expr_D[k] for k in ma_L ] )
		row_count += 1
	
	f.close()
	g.close
	
	return

def get_sample_names( f ):
	"""
	process the file 'samples.txt' to extract sample names
	"""
	train_hts_L = list()
	test_hts_L = list()
	train_ma_L = list()
	test_ma_L = list()

	# file should have two columns
	row_count = 0
	for row in f:
		L = row.strip().split( '\t' )
		if row[0] == 'h' and row_count == 0: # skip the header
			row_count += 1
			continue
		elif row[0] != 'h' and row_count == 0:
			print >> sys.stderr, "Error: missing header in sample names file. Please add a header 'hts<tab>ma'."
			sys.exit( 1 )
		if L[0] == 'NA':
			test_ma_L += [L[1]]
		elif L[0][0] == '*':
			test_hts_L += [L[0].lstrip( '*' )]
			test_ma_L += [L[1]]
		else:
			train_hts_L += [L[0]]
			train_ma_L += [L[1]]
		row_count += 1
	
	return train_hts_L, test_hts_L, train_ma_L, test_ma_L
	
def get_gene_to_tx( f ):
	"""
	create a dictionary of gene to transcripts
	"""
	gene2txs_D = dict()
	for row in f:
		if row[0] == 'g' or row[:3] == 'Ens':
			continue
		L = row.strip().split( '\t' )
		if len( L ) == 1: # just in case
			continue
		if L[0] not in gene2txs_D:
			gene2txs_D[L[0]] = set( [ L[1] ] ) # enforce uniqueness
		else:
			gene2txs_D[L[0]].add( L[1] )
	return gene2txs_D

def make_txs_data( f, g_f, g2tx_D, design="train", ofn="_txs_data.txt.gz", waste=True ):
	"""
	f - training/test data for genes
	g_f - temporary file holding the training/test transcript expression data
	g2tx_D - dictionary for gene-to-transcripts
	"""	
	# waste: a file to dump all genes without probesets and probesets without intensities
	w = gzip.open( design + "_transcripts_missing.txt.gz", 'wb' )
	
	no_genes = 0
	missing_genes = 0
	no_txs = 0
	missing_txs = 0
	
	# put the hts data into a hash table
	rs_data = dict()
	for row in g_f:
		if row[0] == 'g': continue
		L = row.strip().split( '\t' )
		rs_data[ L[0] ] = ",".join( L[1:] )
		
	# combine the data into a pack
	h = gzip.open( design + ofn, 'wb' )
	c = 0
	for row in f: # for each row in the gene data
		no_genes += 1
		if c > 100: break
		if row[0] == 't': continue
		L = row.strip().split( '\t' )
		g = L[0]
		
		# get the associated transcripts
		txs = g2tx_D[ L[0] ] # assumes that all genes have transcripts
		
		# get the tx expr
		rs = list()
		present_txs = list()
		for t in txs:
			try:
				rs += [ rs_data[t] ]
			except KeyError:
				print >> w, t
				missing_txs += 1
				continue
			present_txs += [ t ]
		
		no_txs += len( present_txs )
		
		if len( rs ) == 0:
			missing_genes += 1
			continue
		else:
			print >> h, "\t".join( [ g, ",".join( present_txs ), str( len( present_txs )), L[1], L[2], L[3], ",".join( rs ), L[5] ] )	
		c += 0
	
	return no_genes, missing_genes, no_txs, missing_txs

if __name__ == '__main__':
	unwanted = [ '#' ]
	
	# initialise the parser
	parser = argparse.ArgumentParser( description="Produce training and testing data for MaLTE gene expression models." )
	
	# build the options
	parser.add_argument( '-s', "--samples", default="samples.txt", help="the name of the file containing the map of sample names between HTS and microarray [default: samples.txt]" )
	parser.add_argument( '-r', "--train-data", default="train_data.txt", help="the name of the file containing training data for genes [default: train__data.txt.gz]" )
	parser.add_argument( '-a', "--test-data", default="test_data.txt", help="the name of the file containing test data for genes [default: test__data.txt.gz]" )
	parser.add_argument( '-e', "--txs-expr-data", default="hts_txs_data.txt", help="the name of the file containing a matrix of transcript isoform expression estimates (e.g. from Cufflinks) with first column being EnsEMBL transcript ID and other columns being a sample; the file MUST contain a header of sample names [default: hts_txs_data.txt]" )
	parser.add_argument( '-t', "--gene-transcripts", default="gene_transcripts.txt", help="the name of the file containig the map of gene IDs to transcript IDs [default: gene_transcripts.txt]" )
	
	# if no args given
	if len( sys.argv ) == 1:
		print >> sys.stderr, "Error: missing input files."
		sys.exit( 1 )

	# build vars
	args = parser.parse_args()
	id_fn = args.samples
	train_fn = args.train_data
	test_fn = args.test_data
	hts_txs_fn = args.txs_expr_data
	g2tx_fn = args.gene_transcripts
	
	# check if the files exist and has proper name
	if os.path.exists( id_fn ):
		if not re.search( r"gz$", id_fn ) and not re.search( r"txt$", id_fn ):
			print >> sys.stderr, "Error: Invalid input file '%s': should be .txt or .gz file." % id_fn
			sys.exit( 1 )
		else:
			pass
	else:
		print >> sys.stderr, "Error: Missing sample names file (e.g. '-s samples.txt')."
		sys.exit( 1 )
	
	if os.path.exists( train_fn ):
		if not re.search( r"gz$", train_fn ) and not re.search( r"txt$", train_fn ):
			print >> sys.stderr, "Error: Invalid input file '%s': should be .txt or .gz file." % train_fn
			sys.exit( 1 )
		else:
			pass
	else:
		print >> sys.stderr, "Error: Missing training data file (e.g. '-r train_data.txt')."
		sys.exit( 1 )	
	
	if os.path.exists( test_fn ):
		if not re.search( r"gz$", test_fn ) and not re.search( r"txt$", test_fn ):
			print >> sys.stderr, "Error: Invalid input file '%s': should be .txt or .gz file." % test_fn
			sys.exit( 1 )
		else:
			pass
	else:
		print >> sys.stderr, "Error: Missing test data file (e.g. '-a test_data.txt')."
		sys.exit( 1 )		
	
	if os.path.exists( hts_txs_fn ):
		if not re.search( r"gz$", hts_txs_fn ) and not re.search( r"txt$", hts_txs_fn ):
			print >> sys.stderr, "Error: Invalid input file '%s': should be .txt or .gz file." % hts_txs_fn
			sys.exit( 1 )
		else:
			pass
	else:
		print >> sys.stderr, "Error: Missing transcript expression file (e.g. '-e hts_txs_data.txt')."
		sys.exit( 1 )
	
	if os.path.exists( g2tx_fn ):
		if not re.search( r"gz$", g2tx_fn ) and not re.search( r"txt$", g2tx_fn ):
			print >> sys.stderr, "Error: Invalid input file '%s': should be .txt or .gz file." % g2tx_fn
			sys.exit( 1 )
		else:
			pass
	else:
		print >> sys.stderr, "Error: Missing gene-to-transcript file (e.g. '-t gene_transcripts.txt')."
		sys.exit( 1 )
	
	# raise a warning if the named directory already exists
	try:
		os.mkdir( "temp" )
	except OSError: # and overwrite...
		print >> sys.stderr, "Warning: Found directory 'temp'. Overwriting..."
	
	# get sample names
	print >> sys.stderr, "Creating a map of sample names between datasets from '%s'..." % id_fn,
	if re.search( r"gz$", id_fn ):
		f = gzip.open( id_fn )
	elif re.search( r"txt$", id_fn ):
		f = open( id_fn )
	train_hts_L, test_hts_L, train_ma_L, test_ma_L = get_sample_names( f )
	f.close()
	print >> sys.stderr, "OK"
	print >> sys.stderr, "No. of samples: HTS (train/test); microarray (train/test)... (%s/%s) vs. (%s/%s)" % ( len( train_hts_L ), len( test_hts_L ), len( train_ma_L ), len( test_ma_L ) )
	
	# create a map of gene to transcripts
	if re.search( r"gz$", g2tx_fn ):
		f = gzip.open( g2tx_fn )
	elif re.search( r"txt$", g2tx_fn ):
		f = open( g2tx_fn )
	print >> sys.stderr, "Creating a map of genes to transcripts...",
	g2tx_D = get_gene_to_tx( f )
	f.close()	
	print >> sys.stderr, "OK"
	
	# create the training hts txs data
	p1 = multiprocessing.Process( target=get_data, args=( hts_txs_fn, "temp/train_txs_hts", train_hts_L ), kwargs={ "platform": "hts" })
	p1.start()
	
	p2 = multiprocessing.Process( target=get_data, args=( hts_txs_fn, "temp/test_txs_hts", test_hts_L ),
	kwargs={ "platform": "hts" })
	p2.start() 
	
	p1.join()
	p2.join()
	
	g_train_hts_f = open( "temp/train_txs_hts" )
	g_test_hts_f = open( "temp/test_txs_hts" )
	
	# combine train txs and probe data
	if re.search( r"gz$", train_fn ):
		f = gzip.open( train_fn )
	elif re.search( r"txt$", train_fn ):
		f = open( train_fn )
	print >> sys.stderr, "Compiling training data...",
	train_no_genes, train_missing_genes, train_no_txs, train_missing_txs = make_txs_data( f, g_train_hts_f, g2tx_D, design="train" )
	f.close()
	print >> sys.stderr, "DONE"
	
	# combine test txs and probe data
	if re.search( r"gz$", test_fn ):
		f = gzip.open( test_fn )
	elif re.search( r"txt$", test_fn ):
		f = open( test_fn )
	print >> sys.stderr, "Compiling test data...",
	test_no_genes, test_missing_genes, test_no_txs, test_missing_txs = make_txs_data( f, g_test_hts_f, g2tx_D, design="test" )
	f.close()
	print >> sys.stderr, "DONE"
	
	print >> sys.stderr
	print >> sys.stderr, "Stats: Training data"
	print >> sys.stderr, "Total number of genes       = %s" % train_no_genes
	print >> sys.stderr, "Total genes included        = %s" % ( train_no_genes - train_missing_genes )
	print >> sys.stderr, "Total number of transcripts = %s" % train_no_txs
	print >> sys.stderr, "Total transcripts included  = %s" % ( train_no_txs - train_missing_txs )
	print >> sys.stderr
	print >> sys.stderr, "Stats: Testing data"
	print >> sys.stderr, "Total number of genes       = %s" % test_no_genes
	print >> sys.stderr, "Total genes included        = %s" % ( test_no_genes - test_missing_genes )
	print >> sys.stderr, "Total number of transcripts = %s" % test_no_txs
	print >> sys.stderr, "Total transcripts included  = %s" % ( test_no_txs - test_missing_txs )


	# cleanup
	g_train_hts_f.close()
	g_test_hts_f.close()
	
	# remove temporary directory
	print >> sys.stderr, "Cleaning up...",
	os.unlink( "temp/train_txs_hts" )
	os.unlink( "temp/test_txs_hts" )
	os.rmdir( "temp" )
	print >> sys.stderr, "DONE"
	
	# return exit success
	sys.exit( 0 )









