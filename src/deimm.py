import os
import psutil
import sys
import pickle
from subprocess import PIPE, STDOUT, Popen
import argparse
import socket
from argparse import RawTextHelpFormatter

def run(command, manager = False):
    command = ' '.join([str(i) for i in command])
    if manager:
    	return Popen(command,shell=True,stdout=PIPE,stderr=PIPE)
    else:
	return Popen(command,shell=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='deimm solve',
                                      usage='%(prog)s  [options] --port PORT  --output OUTPUT model_imm.pl model_en.pl',
                                       formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('model', 
                        help='model files (generated from preprocessing)', 
                        nargs=2,
                        type=str)
    parser.add_argument('--port','-p',
                      required = True,
                      type = int,
                      help = "Port")
    parser.add_argument('--output','-o',
                      required = True,
                      help = "Solution output file (ex: output.pcl)")
    parser.add_argument('--approximate','-a',
                      type = float,
                      default = 0.009,
                      required = False,
                      help = "Bound on approximation (default:0.009)")
    parser.add_argument('--relTol','-rel',
                      default = 0.00001,
                      type = float,
                      required = False,
                      help = "The relative tolerance for floating point comparison")
    parser.add_argument('--absTol','-abs',
                      default = 0.001,
                      type = float,
                      required = False,
                      help = "The absolute tolerance for floating point comparison, also used as epsilon in the model")
    parser.add_argument('--threads_per_worker',
                      required = False,
		       default = 4,
                      type = int,
                      help = "Number of threads per worker")
    parser.add_argument('-t',
                        '--threads', 
                        type = int,
                        default = 0,
                        required = False,
                        help = 'Number of threads (default: based on cpu count)')
    parser.add_argument('--verbose','-v',
                      default = 0,
                      type = int,
                      required = False,
                      help = "Verbosity (default 0)")

    args = parser.parse_args()
    
    script_dir = os.path.dirname(os.path.realpath(__file__))
    
    # Get number of threads
    t = args.threads
    if not args.threads:
        t = psutil.cpu_count()
       
    # Calculate number of workers
    if not args.threads_per_worker:
        workers = max(1,int(t/4))
    else:
        workers = max(1,int(t/args.threads_per_worker))
    threads_per_worker = int(t/workers)
    
    # Get host ip
    host_name = socket.gethostname() 
    host_ip = socket.gethostbyname(host_name) 
    
    # Model files
    imm, en = args.model
    
    # Start manager
    print "Starting manager"
    manager_command = script_dir + '/solver/rectangle_splitting_manager.py'
    manager = run(['python', '-u', manager_command,
         	   '-w', workers,
         	   '-p', args.port,
         	   '-a', args.approximate,
         	   '-k', 'rectangle',
         	   '-o', args.output,
         	   '-v', args.verbose], True)
    
    while True:
        line = manager.stdout.readline()
	if 'started' in line or not line:
            if args.verbose:
		print line
	    sys.stdout.flush()
	    break
    
    # Start workers
    print "Starting workers"
    worker_command = script_dir + '/solver/rectangle_splitting_worker.py'
    for i in range(workers):
        worker = run(['python ', worker_command,
                 '-i', imm, en,
                 '-t', threads_per_worker,
                 '-p', args.port,
                 '-m', host_ip,
                 '-a', 'rectangle',
                 '-v', args.verbose])
    worker.wait()
    manager.terminate()