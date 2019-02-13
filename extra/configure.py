#!/usr/bin/python2

import getopt
import sys
import os
import math
from subprocess import call

#TODO bridges
platforms = ["comet", "bridges","stampede"]
compilers = ["intel", "gnu", "pgi", "cray", "ibm"]
options = ['']
configs = { "comet": './configure --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc',
					  "stampede": './configure --enable-fftw --with-fftw=$TACC_FFTW3_DIR FC=mpif90 CC=mpicc'
			}
sourcedir = "p3dfft.3"
destdir = "p3dfft++_compiled"

def usage_exit(msg):
	print msg
	print "USAGE: ./configure.py -s comet|bridges|stampede [-c intel|gnu|pgi|cray|ibm] [-p] [-f extra flags]"
	print "Make sure to run this script from one level above your p3dfft.3 source directory!"
	print "-h displays usage information"
	print "-s specifies which platform"
	print "-c to specify non-default compiler"
	print "-m to build -mt branch"
	print "-p to build performance test"
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
			perf = True
			cflags += ['-O3']
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
	if perf:
		dest += "_p"
	dest = dest + "_" + comp
	cwd = os.getcwd()
	if not os.path.isdir(cwd + '/' + source):
		usage_exit(source + " dir does not exist. Make sure you are at the right directory level")

	# start build
	print configline
	print "Source Directory: " + source
	print "Destination Directory: " + dest
	print "********** Starting build... **********"

	if perf:
		d = cwd + '/' + dest
		try:
			os.mkdir(d)
		except:
			pass
		call('cp -r ' + cwd + '/' + source + '/* ' + d, shell=True)
		os.chdir(d)
		c = configline
		call(c, shell=True)
		call('make', shell=True)
	#TODO Modify once options are available
	else:
		for i in range(pow(2,len(options))):
			d = cwd + '/' + dest + str(i)
			try:
				os.mkdir(d)
			except:
				pass
			call('cp -r ' + cwd + '/' + source + '/* ' + d, shell=True)
			os.chdir(d)
			b = list(bin(i))[2:]
			b = map(int,['0']*(len(options)-len(b)) + b)
			c = configline
			for i in range(len(options)):
				if b[i]:
					c += ' --enable-' + options[i]
			print "Configuring " + d + " with "
			print "\t" + c
			c += " > config_output"
			ret = call(c, shell=True)
			if ret != 0:
				usage_exit("CONFIG FAILED! CHECK config_output for log")
			print "Configured " + d
			print "Making " + d
			ret = call('make > make_output', shell=True)
			if ret != 0:
				usage_exit("MAKE FAILED! CHECK make_output for log")
			print "Built " + d

	print "********** Done. **********"

if __name__ == '__main__':
	main()

