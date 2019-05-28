#!/usr/bin/python2

import math
import getopt
import sys
import os
import re
from time import strftime, localtime
from fractions import Fraction
from itertools import product, permutations

platforms = ["comet", "stampede", "bridges"]
one_dim_perms = map(lambda a: a.replace('(', '').replace(')', '').replace(',', '')
										,[str(p) for p in product(permutations([0,1,2]), repeat=2)])
current_time = strftime("%d-%m-%Y-%H%M%S", localtime())
# job size
NUMNODES = 1		# --nodes
TASKSPERNODE = 16   # --ntasks-per-node

MT_NUMTHREADS = 2	  # env var OMP_NUM_THREADS
MT_RANKSPERNODE = 8	# used for dims

assert(MT_NUMTHREADS * MT_RANKSPERNODE == TASKSPERNODE)

def usage_exit(msg):
	print msg
	print "USAGE: ./configure.py -s comet|gordon|edison|cori|stampede [-m] [-p MINCORES MAXCORES [NUMTHREADS if -m is used] [MINGRID MAXGRID]]"
	print "Make sure to run this script from one level above your p3dfft.3 source directory!"
	print "-h displays usage information"
	print "-d source directory created from configure.py"
	print "-s  specifies which platform"
	print "-m  to build -mt branch"
	print "-p  to build performance test"
	sys.exit(1)

# Return list of all tests of this specific configuration
def gettests(sd):
	p3dfft_dirs = next(os.walk(sd))[1]
	pattern = re.compile('p3dfft\+\+_compiled\S+')
	p3dfft_dirs = sorted(filter(pattern.match, p3dfft_dirs))
	pattern = re.compile('test.+_([cf]|cpp)$')
	all_tests = []
	for p3dir in p3dfft_dirs:
		for root, _, files in os.walk(os.path.join(sd, p3dir, 'sample')):
			if root == os.path.join(p3dir, 'sample'):
				continue
			all_tests += [os.path.abspath(os.path.join(root, f)) for f in filter(pattern.match, files)]
	return all_tests


# Get dimensions based on whether mt version or not
def getdims(mt):
	if mt:
		n = MT_RANKSPERNODE
	else:
		n = TASKSPERNODE
	dims = []
	facs = [i for i in range(1, n+1) if n%i == 0]
	if (len(facs) % 2 == 0):
					# take the two factors in the middle
		dims.append("'" + str(facs[len(facs)/2-1]) + " " + str(facs[len(facs)/2]) + "'")
	else:
		# perfect square, take the factor in the middle
		dims.append("'" + str(facs[len(facs)/2]) + " " + str(facs[len(facs)/2]) + "'")
	dims.append("'1 " + str(n) + "'")
	dims.append("'" + str(n) + " 1'")
	return dims

# Standard test run line
def runline(platform, mt, output_dir, test):
	r = ''
	if platform == "comet":
		if mt:
			r = "ibrun -n " + str(MT_RANKSPERNODE) + " " + test
		else:
			r = "ibrun -n " + str(TASKSPERNODE * NUMNODES) + " " + test
	elif platform == "stampede":
		if mt:
			r = "ibrun -n " + str(MT_RANKSPERNODE) + " -o 0 tacc_affinity " + test
		else:
			r = "ibrun -n " + str(TASKSPERNODE * NUMNODES) + " -o 0 " + test
	elif platform == "bridges":
		if mt:
			r = "mpirun -n " + str(MT_RANKSPERNODE) + " " + test
		else:
			r = "mpirun -n " + str(TASKSPERNODE * NUMNODES) + " " + test
	return r + " &>> " + os.path.join(output_dir, "output_" + os.path.basename(test)) + "\n"

# Test for 1x1 dims
def onebyone(platform, mt, output_dir, test):
	r = ''
	if platform == "comet":
		r = "ibrun --npernode 1 " + test
	elif platform == "stampede":
		if mt:
			r = "ibrun -n 1 -o 0 tacc_affinity " + test
		else:
			r = "ibrun -n 1 -o 0 " + test
	elif platform == "bridges":
		if mt:
			rt = "mpirun -n 1 " + test
		else:
			r = "mpirun -n 1" + test
	return r + " &>> " + os.path.join(output_dir, "output_" + os.path.basename(test)) + "\n"

# Write all tests for all dims
def buildall(platform, mt, all_tests, all_dims, batchf, output_dir):
	for test in all_tests:
		#if "cheby" in test:
		#	batchf.write("echo '32 32 33 2 1' > stdin\n")
		#elif "many" in test:
		#	batchf.write("echo '32 32 32 2 5 1' > stdin\n")
		#elif "pruned" in test:
		#	batchf.write("echo '64 64 64 32 32 32 2 1' > stdin\n")
		#else:
		#	batchf.write("echo '128 128 128 2 1' > stdin\n")
		basename = os.path.basename(test)
		if '3D' in basename:
			batchf.write("echo '128 128 128 2 1' > stdin\n")
			for dims in all_dims:
				batchf.write("echo " + dims + " > dims\n")
				batchf.write(runline(platform,mt,output_dir,test))
			batchf.write("echo '1 1' > dims\n")
			batchf.write(onebyone(platform, mt, output_dir, test))
		elif '1D' in basename:
			batchf.write("rm -f dims\n")
			for dims in all_dims:
				batchf.write("echo " + dims + " > dims\n")
				for perm in one_dim_perms:
					dim_in = perm.find('0')/2
					dim_out = perm.find('0', 5)/2 - 3
					batchf.write("echo -e '128 128 128 " + str(dim_in) + ' 1\\n' + perm[:5] + '\\n' + perm[6:] + "' > trans.in\n")
					batchf.write(runline(platform, mt, output_dir, test))
					batchf.write("echo -e '128 128 128 " + str(dim_out) + ' 1\\n' + perm[:5] + '\\n' + perm[6:] + "' > trans.in\n")
					batchf.write(runline(platform, mt, output_dir, test))
		elif 'IDIR' in basename:
			batchf.write("echo '128 128 128 2 1 3' > stdin\n")
			for dims in all_dims:
				batchf.write("echo " + dims + " > dims\n")
				batchf.write(runline(platform, mt, output_dir, test))
			batchf.write("echo '1 1' > dims\n")
			batchf.write(onebyone(platform, mt, output_dir, test))
		# 1x1 dims test

def unevengrid(platform, mt, all_tests, all_dims, batchf, output_dir):
	for test in all_tests:
		if "cheby" in test:
			batchf.write("echo '14 26 38 2 1' > stdin\n")
		elif "many" in test:
			batchf.write("echo '14 26 38 2 5 1' > stdin\n")
		elif "pruned" in test:
			batchf.write("echo '64 64 64 32 32 32 2 1' > stdin\n")
		else:
			batchf.write("echo '14 26 38 2 1' > stdin\n")
		# write dims
		batchf.write("echo " + all_dims[0] + " > dims\n")
		# run test
		batchf.write(runline(platform, mt, output_dir, test))

# Test for performance
def perftest(platform, mt, test, curr_numcores, num_threads):
	if platform == "comet":
		if mt:
			return "ibrun -n " + str(curr_numcores/num_threads) + " " + test + "\n"
		else:
			return "ibrun -N " + str(curr_numcores) + " " + test + "\n"
	elif platform == "stampede":
		if mt:
			return "ibrun -n " + str(curr_numcores/num_threads) + " -o 0 tacc_affinity " + test + "\n"
		else:
			return "ibrun -n " + str(curr_numcores) + " -o 0 " + test + "\n"
	elif platform == "bridges":
		if mt:
			return "mpirun -n " + str(curr_numcores/num_threads) + " " + test + "\n"
		else:
			return "mpirun -n " + str(curr_numcores) + " " + test + "\n"

# Write all tests for performance testing
#TODO NOT WORKING
def runperf(platform, mt, perf, batchf, MINCORES, MAXCORES, MINGRID, MAXGRID, PERF_NUMTHREADS):
	# Get test_sine
	basedir = os.getcwd()
	p3dfft_dir = os.path.join(basedir, "p3dfft++_compiled_p_")
	f_test = os.path.join(p3dfft_dir, 'sample/C++/test_sin_cpp')

	# Run test_sine for all cores arranged in all dims
	curr_numcores = MINCORES
	while curr_numcores <= MAXCORES:
		# Calculate maximum grid size based on cores, and initialise grid size.
		NMAX = int(math.floor(math.pow(curr_numcores*4*math.pow(10,9)/48,Fraction(1,3))))
		if MINGRID > NMAX:
			print("MINGRID > NMAX")
			sys.exit(-1)

		curr_gridsize = MINGRID

		#while (curr_gridsize < NMAX)
		while curr_gridsize <= MAXGRID: # TODO: comment me out later
			batchf.write("echo '" + str(curr_gridsize) + " " + str(curr_gridsize) + " " + str(curr_gridsize) + " 2 5' > stdin\n")

			# Calculate dims
			all_dims = []
			all_dims.append("'16 " + str(int(curr_numcores/(16*PERF_NUMTHREADS))) + "'")
			all_dims.append("'32 " + str(int(curr_numcores/(32*PERF_NUMTHREADS))) + "'")
			all_dims.append("'64 " + str(int(curr_numcores/(64*PERF_NUMTHREADS))) + "'")
			#all_dims = ["'16 NUMCORES/16'", "'32 NUMCORES/32'", "'64 NUMCORES/64'"]

			for dims in all_dims:
				# write dims
				batchf.write("echo " + dims + " > dims\n")
				# run test
				batchf.write(perftest(platform, mt, f_test, curr_numcores, PERF_NUMTHREADS))
			curr_gridsize *= 2
		curr_numcores *= 2


def script_header(platform, batchf, mt, perf, email, output_dir,sd):
	batchf.write('#!/bin/bash\n')
	if platform == "comet":
		batchf.write('#SBATCH --job-name="' + "p3dfft++_compiled" + '"\n')
		batchf.write('#SBATCH --output="' + os.path.join(output_dir,'out.%j') + '"\n')
		batchf.write('#SBATCH --partition=compute\n')
		if perf:
			batchf.write('#SBATCH --nodes=' + str(int(MAXCORES/(32*PERF_NUMTHREADS))) + '\n')
			batchf.write('#SBATCH --ntasks-per-node=32\n')
		else:
			batchf.write('#SBATCH --nodes=' + str(NUMNODES) + '\n')
			batchf.write('#SBATCH --ntasks-per-node=' + str(TASKSPERNODE) + '\n')
		batchf.write('#SBATCH --export=ALL\n')
		batchf.write('#SBATCH --switches=1\n')
		if email:
			batchf.write('#SBATCH --mail-user=' + email + '\n')
			batchf.write('#SBATCH --mail-type=ALL\n')
		batchf.write('#SBATCH -t 01:00:00\n')
	elif platform == "stampede":
		batchf.write('#SBATCH -J ' + "p3dfft++_compiled" + '\n')
		batchf.write('#SBATCH -o out/out.%j\n')
		batchf.write('#SBATCH -e out/out.%j\n')
		batchf.write('#SBATCH -p normal\n')
		if perf:
			batchf.write('#SBATCH -N ' + str(int(MAXCORES/(16*PERF_NUMTHREADS))) + '\n')
		else:
			batchf.write('#SBATCH -N ' + str(NUMNODES) + '\n')
		if perf:
			batchf.write('#SBATCH --ntasks-per-node ' + str(MAXCORES/PERF_NUMTHREADS) + '\n')
		elif mt:
			batchf.write('#SBATCH --ntasks-per-node ' + str(MT_RANKSPERNODE) + '\n')
		else:
			batchf.write('#SBATCH --ntasks-per-node ' + str(TASKSPERNODE) + '\n')
		if email:
			batchf.write('#SBATCH --mail-user=' + email + '\n')
			batchf.write('#SBATCH --mail-type=all\n')
		batchf.write('#SBATCH -t 01:00:00\n')
	elif platform == "bridges":
		batchf.write('#SBATCH --job-name="' + "p3dfft++_compiled" + '"\n')
		batchf.write('#SBATCH --output="' + os.path.join(output_dir,'out.%j') + '"\n')
		batchf.write('#SBATCH --partition=RM\n')
		if perf:
			batchf.write('#SBATCH --nodes=' + str(int(MAXCORES/(32*PERF_NUMTHREADS))) + '\n')
			batchf.write('#SBATCH --ntasks-per-node=32\n')
		else:
			batchf.write('#SBATCH --nodes=' + str(NUMNODES) + '\n')
			batchf.write('#SBATCH --ntasks-per-node=' + str(TASKSPERNODE) + '\n')
		batchf.write('#SBATCH --export=ALL\n')
		batchf.write('#SBATCH --switches=1\n')
		if email:
			batchf.write('#SBATCH --mail-user=' + email + '\n')
			batchf.write('#SBATCH --mail-type=ALL\n')
		batchf.write('#SBATCH -t 01:00:00\n')
	batchf.write('cd ' + sd + '\n')
	batchf.write('\n')

def main():
	platform = None
	source_dir = None
	mt = False
	perf = False
	uneven = False
	email = ''
	perfopts = ""

	# parse command line options
	try:
		opts = getopt.getopt(sys.argv[1:], 'd:s:e:hmp:u')
	except getopt.GetoptError as err:
		usage_exit(str(err))
	for o, a in opts[0]:
		if o == '-s':
			platform = a
		elif o == '-d':
			source_dir = a
		elif o == '-e':
			email = a
		elif o == '-h':
			usage_exit("**Help Message**")
		elif o == '-m':
			pass
			#mt = True
		elif o == '-p':
			pass
			#perf = True
			#perfopts = a
		elif o == '-u':
			uneven = True
		else:
			assert False, "unhandled option"
	if platform == None:
		usage_exit("no platform specified")
	elif platform not in platforms:
		usage_exit("invalid platform specified")
	elif source_dir == None:
		usage_exit("no source directory supplied")
	elif not os.path.isdir(source_dir):
		usage_exit("source directory doesn't exist")

	if uneven and perf:
		usage_exit("no uneven grid performance tests")

	if perf:
		perfopts = perfopts.split()
		MINCORES = int(perfopts[0])
		MAXCORES = int(perfopts[1])
		if mt:
			if len(perfopts) != 3 and len(perfopts) != 5:
				usage_exit("num threads not specified")
			if len(perfopts) < 5:
				print("Using 1024^3 as default grid size.")
				MINGRID = 1024
				MAXGRID = 1024
			else:
				MINGRID = int(perfopts[3])
				MAXGRID = int(perfopts[4])
			PERF_NUMTHREADS = int(perfopts[2])
		else:
			if len(perfopts) < 4:
				print("Using 1024^3 as default grid size.")
				MINGRID = 1024
				MAXGRID = 1024
			else:
				MINGRID = int(perfopts[2])
				MAXGRID = int(perfopts[3])
			PERF_NUMTHREADS = 1

		if MINCORES > MAXCORES:
			print("MINCORES > MAXCORES")
			sys.exit(-1)
		if MINGRID > MAXGRID:
			print("MINGRID > MAXGRID")
			sys.exit(-1)

	# start write job	jobs/comet-p3dfft-mt.sh
	jobs_dir = os.path.join(os.getcwd(), "jobs_" + current_time)
	fname = os.path.join(jobs_dir, platform + "-" + "p3dfft++_tests" + os.path.basename(source_dir) + ".sh")
	try:
		os.mkdir(jobs_dir)
	except IOError:
		pass

	batchf = open(fname, 'w')

	output_dir = os.path.join(jobs_dir, 'out')
	try:
		os.mkdir(output_dir)
	except IOError:
		pass

	# write header
	script_header(platform, batchf, mt, perf, email, output_dir,jobs_dir)
	if mt:
		if perf:
			batchf.write('export OMP_NUM_THREADS=' + str(PERF_NUMTHREADS) + '\n')
		else:
			batchf.write('export OMP_NUM_THREADS=' + str(MT_NUMTHREADS) + '\n')

	if perf:
		runperf(platform, mt, perf, batchf, MINCORES, MAXCORES, MINGRID, MAXGRID, PERF_NUMTHREADS)
	else:
		all_tests = gettests(source_dir)
		all_dims = getdims(mt)
		if uneven:
			unevengrid(platform, mt, all_tests, all_dims, batchf, output_dir)
		else:
			buildall(platform, mt, all_tests, all_dims, batchf, output_dir)
	batchf.write("grep -rcwE * -e 'Error|incorrect|BAD'\n")
	# Close the file. Done.
	batchf.close()

	print "Wrote to " + fname

if __name__ == '__main__':
	main()
