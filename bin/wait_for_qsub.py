#!/usr/bin/env python
import time
from subprocess import Popen, PIPE

if __name__ == '__main__' :

	# this is gross, but it works when you need to stall a pipeline until all your jobs are done
	done = False
	while not done :
		qstat_output = Popen('qstat',shell=True,stdout=PIPE).communicate()[0]
		if qstat_output == '' :
			done = True
		else :
			time.sleep(1)
