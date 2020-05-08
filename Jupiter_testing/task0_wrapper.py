from ctypes import cdll
from ctypes import c_char_p
import ctypes 
import os
import time
import shutil
import math
import random
import subprocess

def task(filename,pathin,pathout):
     filename= "task0.c"
     subprocess.call(["gcc","-o","libtest.so","-shared","-fPIC","task0.c"])
    


     lib = cdll.LoadLibrary("./libtest.so")
     string_buffers = [ctypes.create_string_buffer(8) for i in range(4)]
     pointers = (ctypes.c_char_p*4)(*map(ctypes.addressof, string_buffers))
     lib.main(pointers)
     results = [s.value for s in string_buffers]
     print (results)

def main():
	filelist = '1botnet.ipsum'
	outpath = os.path.join(os.path.dirname(__file__), 'sample_input/')
	outfile = task(filelist, outpath, outpath)
	return outfile
    
if __name__ == "__main__":
    main()
