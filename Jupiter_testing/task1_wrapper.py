import os
import time
import shutil
import math
import random
import subprocess

def task():
     filename= "task1.c"
     subprocess.call(["gcc",filename])
     subprocess.call("./a.out")
        

def main():
    task()
    
if __name__ == "__main__":
    main()
