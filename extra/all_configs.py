#!/usr/bin/python2
#This script essentially builds all permuations for P3DFFT++ based on the configure options available
#It does this by copying all the files from the source directory into a new directory for each build

import getopt
import sys
import os
import math
from subprocess import call
from itertools import combinations
from time import strftime, localtime


platforms = ["comet", "bridges","stampede"]
compilers = ["intel", "gnu", "pgi", "cray", "ibm"]
#Add configure options to this list once there are some
options = []
configs = { "comet": './configure --enable-fftw',
			"stampede": './configure --enable-fftw --with-fftw-lib=$TACC_FFTW3_LIB --with-fftw-inc=$TACC_FFTW3_INC',
			"bridges": './configure --enable-fftw CXX=mpiicpc CC=mpiicc FC=mpiifort'
			}
sourcedir = "p3dfft.3"
destdir = "p3dfft++_configs_" + strftime("%d-%m-%Y-%H%M%S", localtime())

def usage_exit(msg):
	print msg
	print "USAGE: ./all_configs.py -s comet|bridges|stampede [-c intel|gnu|pgi|cray|ibm] [-p] [-f extra flags]"
	print "Make sure to run this script from one level above your p3dfft.3 source directory!"
	print "-h displays usage information"
	print "-s specifies which platform"
	print "-c to specify non-default compiler"
	#print "-m to build -mt branch"
	#print "-p to build performance test"
	print "-f extra configure flags"
	sys.exit(1)

def main():
	cflags = []
	platform = ''
	comp = ''
	extra = ''
	#mt = False
	perf = False

	# parse command line options
	try:
		opts = getopt.getopt(sys.argv[1:], 'hs:c:pf:')
	except getopt.GetoptError as err:
		usage_exit(str(err))
	for o, a in opts[0]:
		if o == '-h':
			usage_exit("**Help Message**")
		if o == '-s':
			platform = a
		elif o == '-c':
			comp = a
		elif o == '-p':
			#perf = True
			#cflags += ['-O3']
			pass
		elif o == '-f':
			extra = a
		else:
			usage_exit( "unhandled option")
	if platform == None:
		usage_exit("no platform specified")
	elif platform not in platforms:
		usage_exit("invalid platform specified")

	# make configline according to compiler
	configline = configs[platform]
	if comp and comp not in compilers:
		usage_exit("invalid compiler specified")
	if comp:
		configline += " --enable-" + comp
	else:
		comp = "intel"
	if extra != None:
		configline += " " + extra
	if cflags:
		configline += ' CFLAGS=\" '
		for flag in cflags:
			configline += flag + ' '
		configline += '\"'

	# ensure that the source dir exists
	source = sourcedir
	dest = destdir
	dest = dest + "_" + comp
	if perf:
		dest += "_p"
	cwd = os.getcwd()
	if not os.path.isdir(cwd + '/' + source):
		usage_exit(source + " dir does not exist. Make sure you are at the right directory level")

	# start build
	print configline
	print "Source Directory: " + source
	print "Destination Directory: " + dest
	print "********** Starting build... **********"

	# Perf is not working!
	if perf:
		d = os.path.join(cwd,dest)
		dd = os.path.join(d, 'p3dfft++_compiled_p')
		try:
			os.mkdir(d)
			os.mkdir(dd)
		except Exception as e:
			print e
			sys.exit(1)
		call('cp -r ' + cwd + '/' + source + '/* ' + dd, shell=True)
		os.chdir(dd)
		c = configline
		call(c, shell=True)
		call('make', shell=True)
	#TODO Modify once options are available
	else:
		combos = []
		d = os.path.join(cwd,dest)
		try:
			os.mkdir(d)
		except Exception as e:
			print e
			sys.exit(1)
		# Generates a list of combinations of the options list
		# Then iterates through that list and builds p3dfft in a new directory for each combination
		for i in range(len(options)+1):
			combos += list(combinations(options, i))
		for combo in combos:
			dd = os.path.join(d, 'p3dfft++_compiled_' + '-'.join(combo))
			try:
				os.mkdir(dd)
			except Exception as e:
				print e
				sys.exit(1)
			call('cp -r ' + cwd + '/' + source + '/* ' + dd, shell=True)
			os.chdir(dd)
			c = configline
			for o in combo:
				c += ' --enable-' + o
			print "Configuring " + dd + " with "
			print "\t" + c
			c += " &> config_output"
			ret = call(c, shell=True)
			if ret != 0:
				print "CONFIG FAILED! CHECK config_output for log"
				continue
			print "Configured " + dd + " successfully"
			print "Making " + dd
			ret = call('make &> make_output', shell=True)
			if ret != 0:
				print "MAKE FAILED! CHECK make_output for log"
				continue
			print "Built " + dd + " successfully"

	print "********** Done. **********"

if __name__ == '__main__':
	main()
