#!/usr/bin/env python

#	prepare_data.py

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
		if (row[0] == 'p' or row[0] == 'g') and row_count == 0:
			names_L = L[index:]
			if platform == "ma":
				print >> g, "\t".join( [ L[0], L[column] ] + ma_L )
			elif platform == "hts":
				if len( ma_L ) == 0:
					print >> g,"gene_id\tsamples"
				else:
					print >> g, "\t".join( [ L[0] ] + ma_L )
			row_count += 1
			continue
		elif row[0] != 'p' and row[0] != 'g' and row_count == 0:
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
				print >> g,"%s\tNA" % L[0]
			else:
				print >> g, "\t".join( [ L[0] ] + [ expr_D[k] for k in ma_L ] )
		row_count += 1
	
	f.close()
	g.close()
	
	return

def get_sample_names( f, extract_PCs=False ):
	"""
	process the file 'samples.txt' to extract sample names
	"""
	train_hts_L = list()
	test_hts_L = list()
	train_ma_L = list()
	test_ma_L = list()
	
	if extract_PCs:
		train_PCs = open( "train_PCs.txt", 'w' )
		test_PCs = open( "test_PCs.txt", 'w' )
	
	# file should have two columns
	row_count = 0
	for row in f:
		L = row.strip().split( '\t' )
		if row[0] == 'h' and row_count == 0: # skip the header
			row_count += 1
			if len( L ) > 2:
				PCs_present = True 
			else:
				PCs_present = False
			continue
		elif row[0] != 'h' and row_count == 0:
			print >> sys.stderr, "Error: missing header in sample names file. Please add a header 'hts<tab>ma'."
			sys.exit( 1 )
		if L[0] == "*NA":
			test_ma_L += [L[1]]
			if PCs_present and extract_PCs:	
				print >> test_PCs, "\t".join( L[2:] )
		elif L[0][0] == '*':
			test_hts_L += [L[0].lstrip( '*' )]
			test_ma_L += [L[1]]
			if PCs_present and extract_PCs:
				print >> test_PCs, "\t".join( L[2:] )
		else:
			train_hts_L += [L[0]]
			train_ma_L += [L[1]]
			if PCs_present and extract_PCs:
				print >> train_PCs, "\t".join( L[2:] )
		row_count += 1
	
	return train_hts_L, test_hts_L, train_ma_L, test_ma_L

def get_gene_to_ps( f ):
	"""
	create a dictionary of gene to probesets
	"""
	g2ps_D = dict()
	for row in f:
		if row[0] == 'g' or row[:3] == "Ens":
			continue
		L = row.strip().split( '\t' )
		if len( L ) == 1:
			continue
		if L[1] == '.':
			continue
		g = L[0]; ps = int( L[1] )
		if g not in g2ps_D:
			g2ps_D[g] = set([ps])
		else:
			g2ps_D[g].add(ps)
	
	return g2ps_D

def make_data( f, g_f, g2ps_D, design="train", ofn="_data.txt.gz", waste=True ):
	"""
	take raw data and create training and test data
	f - ma data
	g_f - hts data
	g2ps_D - map of gene to probesets
	"""
	unwanted = [ '#' ]
	no_samples = 0
	
	# some stats
	no_genes = 0
	no_genes_included = 0
	no_probesets = 0
	no_probesets_included = 0
	no_probes = 0
	no_probes_included = 0
	
	# waste: a file to dump all genes without probesets and probesets without intensities
	w1 = gzip.open( design + "_genes_missing.txt.gz", 'wb' )
	w2 = gzip.open( design + "_probesets_missing.txt.gz", 'wb' )
	
	# put the probe data into a hash table
	ps2int = dict()
	for row in f:
		if row[0] in unwanted: continue
		L = row.strip().split( '\t' )
		if row[0] == 'p':
			sample_names = L[2:]
			continue
		ps = int(L[1])
		p = L[0]
		if no_samples == 0:
			no_samples = len( L[2:] )
	
		ints = ",".join( L[2:] )
	
		if ps not in ps2int: # ps are unique
			ps2int[ps] = dict()
	
		ps2int[ps][p] = ints
		no_probes += 1
	
	# put the hts data into a hash table
	rs_data = dict()
	for row in g_f:
		if row[0] == 'g': continue
		L = row.strip().split( '\t' )
		rs_data[ L[0] ] = ",".join( L[1:] )
	
	# combine the data into a pack
	h = gzip.open( design + ofn, 'wb' )
	
	c = 0
	for g in g2ps_D:
		no_genes += 1
		if c > 100: break
	
		# get the RNA-Seq values
		try:
			rs = rs_data[g]
		except KeyError:
			print >> w1, g
			continue
			
		# get probe intensities
		probes = ""
		intensities = ""
		for ps in g2ps_D[g]:
			no_probesets += 1
			try:
				probe_ids = ps2int[ps].keys()
			except KeyError:
				print >> w2, ps
				continue
			for p in probe_ids:
				no_probes_included += 1
				if probes == "":
					probes = p
					intensities = ps2int[ps][p]
				else:
					probes += "," + p
					intensities += "," + ps2int[ps][p]
			
			if probes != "":
				no_probesets_included += 1

		if probes == "":
			continue
#			print >> h, "\t".join( [ g, str( no_samples ), "0", "NA", rs, "NA" ] )
		else:
			print >> h, "\t".join( [ g, str( no_samples ), str( len( probes.split( ',' )) ), probes, rs, intensities ] )
			no_genes_included += 1
			c += 0
	
	h.close()
	
	return no_genes, no_genes_included, no_probesets, no_probesets_included, no_probes, no_probes_included

if __name__ == '__main__':
	unwanted = [ '#' ]
	
	# initialise the parser
	parser = argparse.ArgumentParser( description="Produce training and testing data for MaLTE gene expression models." )
	
	# build the options
	parser.add_argument( '-s', "--samples", default="samples.txt", help="the name of the file containing the map of sample names between HTS and microarray [default: samples.txt]" )
	parser.add_argument( '-t', "--hts-data", default="hts_data.txt", help="the name of the file containing a matrix of HTS expression estimates [default: hts_data.txt]" )
	parser.add_argument( '-m', "--ma-data", default="ma_data.txt", help="the name of the file containing a matrix of microarray probe intensities [default: ma_data.txt]" )
	parser.add_argument( '-g', "--gene-probesets", default="gene_probesets.txt", help="the name of the file containig the map of gene names to probeset names [default: gene_probesets.txt]" )
	parser.add_argument( '-r', "--raw", action='store_true', default=False, help="should you used output directly from APT apt-cel-extract? [default: False]" )
	parser.add_argument( '-p', "--principal-components", action='store_true', default=False, help="should you extract principal components? (if present: see documentation on this) [defaultt: False]" )
	
	# if no args given
	if len( sys.argv ) == 1:
		print >> sys.stderr, "Error: missing input files."
		sys.exit( 1 )
	
	# build vars
	args = parser.parse_args()
	id_fn = args.samples
	ma_fn = args.ma_data
	hts_fn = args.hts_data
	g2ps_fn = args.gene_probesets
	raw = args.raw
	extract_PCs = args.principal_components
	
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
	
	if os.path.exists( ma_fn ):
		if not re.search( r"gz$", ma_fn ) and not re.search( r"txt$", ma_fn ):
			print >> sys.stderr, "Error: Invalid input file '%s': should be .txt or .gz file." % ma_fn
			sys.exit( 1 )
		else:
			pass
	else:
		print >> sys.stderr, "Error: Missing probe intensity file (e.g. '-m ma_data.txt' or  '-r -m raw_ma_data.txt')."
		sys.exit( 1 )	
	
	if os.path.exists( hts_fn ):
		if not re.search( r"gz$", hts_fn ) and not re.search( r"txt$", hts_fn ):
			print >> sys.stderr, "Error: Invalid input file '%s': should be .txt or .gz file." % hts_fn
			sys.exit( 1 )
		else:
			pass
	else:
		print >> sys.stderr, "Error: Missing HTS data file (e.g. '-t hts_data.txt')."
		sys.exit( 1 )		
	
	if os.path.exists( g2ps_fn ):
		if not re.search( r"gz$", g2ps_fn ) and not re.search( r"txt$", g2ps_fn ):
			print >> sys.stderr, "Error: Invalid input file '%s': should be .txt or .gz file." % g2ps_fn
			sys.exit( 1 )
		else:
			pass
	else:
		print >> sys.stderr, "Error: Missing gene-to-(probe set) file (e.g. '-g gene_probesets.txt')."
		sys.exit( 1 )
	
	# raw microarray file (direct from APT)?
	if raw:
		column = 4
	else:
		column = 1
	
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
	train_hts_L, test_hts_L, train_ma_L, test_ma_L = get_sample_names( f, extract_PCs )
	f.close()
	print >> sys.stderr, "OK"
	print >> sys.stderr, "No. of samples: HTS (train/test); microarray (train/test)... (%s/%s) vs. (%s/%s)" % ( len( train_hts_L ), len( test_hts_L ), len( train_ma_L ), len( test_ma_L ) )
	if extract_PCs:
		print >> sys.stderr, "Extracted principal components present to files 'train_PCs.txt' and 'test_PCs.txt'."
	
	# get map of genes to ps
	print >> sys.stderr, "Creating a map of genes to probe sets from '%s'..." % g2ps_fn,
	if re.search( r"gz$", g2ps_fn ):
		f = gzip.open( g2ps_fn )
	elif re.search( r"txt$", g2ps_fn ):
		f = open( g2ps_fn )
	g2ps_D = get_gene_to_ps( f )
	f.close()
	print >> sys.stderr, "OK"

	print >> sys.stderr, "Creating temporary training and test data for HTS and microarray...",
#	# create the training probe data
	p1 = multiprocessing.Process( target=get_data, args=( ma_fn, "temp/train_ma", train_ma_L ), kwargs={ "platform": "ma", "column": column } )
	p1.start()
	
	# create the test probe data
	p2 = multiprocessing.Process( target=get_data, args=( ma_fn, "temp/test_ma", test_ma_L ),  kwargs={ "platform": "ma", "column": column } )
	p2.start()

	p3 = multiprocessing.Process( target=get_data, args=( hts_fn, "temp/train_hts", train_hts_L ), kwargs={ "platform": "hts" } )
	p3.start()
	
	# create the test hts data
	p4 = multiprocessing.Process( target=get_data, args=( hts_fn, "temp/test_hts", test_hts_L ), kwargs={ "platform": "hts" } )
	p4.start()
	
	# wait until you're done
	p1.join()
	p2.join()
	p3.join()
	p4.join()
	
	print >> sys.stderr, "OK"
#	sys.exit( 1 )
	
	print >> sys.stderr, "Extracted microarray and HTS data"

	g_train_ma_f = open( "temp/train_ma" )
	g_test_ma_f = open( "temp/test_ma" )
	g_train_hts_f = open( "temp/train_hts" )
	g_test_hts_f = open( "temp/test_hts" )

	# combine training data
	print >> sys.stderr, "Compiling training data...",
	tr_genes_i, tr_genes_present_i, tr_probesets_i, tr_probesets_present_i, tr_probes_i, tr_probes_present_i = make_data( g_train_ma_f, g_train_hts_f, g2ps_D, design="train", waste=True )
	print >> sys.stderr, "DONE"
	
	print >> sys.stderr
	print >> sys.stderr, "Stats: Training data"
	print >> sys.stderr, "Total number of genes     = %s" % tr_genes_i
	print >> sys.stderr, "Total genes included      = %s" % tr_genes_present_i
	print >> sys.stderr, "Total number of probesets = %s" % tr_probesets_i
	print >> sys.stderr, "Total probesets included  = %s" % tr_probesets_present_i
	print >> sys.stderr, "Total number of probes    = %s" % tr_probes_i
	print >> sys.stderr, "Total probes included     = %s"	% tr_probes_present_i
	print >> sys.stderr

	# combine testing data
	print >> sys.stderr, "Compiling test data...",
	te_genes_i, te_genes_present_i, te_probesets_i, te_probesets_present_i, te_probes_i, te_probes_present_i = make_data( g_test_ma_f, g_test_hts_f, g2ps_D, design="test", waste=True )
	print >> sys.stderr, "DONE"
	
	print >> sys.stderr, "Stats: Testing data"
	print >> sys.stderr, "Total number of genes     = %s" % te_genes_i
	print >> sys.stderr, "Total genes included      = %s" % te_genes_present_i
	print >> sys.stderr, "Total number of probesets = %s" % te_probesets_i
	print >> sys.stderr, "Total probesets included  = %s" % te_probesets_present_i
	print >> sys.stderr, "Total number of probes    = %s"	% te_probes_i
	print >> sys.stderr, "Total probes included     = %s"	% te_probes_present_i
	
	# clean up
	g_train_ma_f.close()
	g_test_ma_f.close()
	g_train_hts_f.close()
	g_test_hts_f.close()
		
#	# remove temporary directory
	print >> sys.stderr, "Cleaning up...",
	os.unlink( "temp/train_ma" )
	os.unlink( "temp/test_ma" )
	os.unlink( "temp/train_hts" )
	os.unlink( "temp/test_hts" )
	os.rmdir( "temp" )
	print >> sys.stderr, "DONE"
	
	# return exit success
	sys.exit( 0 )
