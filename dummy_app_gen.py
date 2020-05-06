"""
This code uses the output DAG of ``rand_task_gen.py``(by Diyi), which is used to generate a random DAG, to generate the corresponding dummy application working for Jupiter (version 3) , 

"""
__author__ = "Quynh Nguyen, Jiatong Wang, Diyi Hu and Bhaskar Krishnamachari"
__copyright__ = "Copyright (c) 2019, Autonomous Networks Research Group. All rights reserved."
__license__ = "GPL"
__version__ = "1.0"


import argparse
import numpy as np
import random
from functools import reduce
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import pylab as plt
import yaml
import os
import json
import shutil
from collections import defaultdict
import random
 
EPSILON = 1e-2

def parse_args(conf_file):
	parser = argparse.ArgumentParser(description='generate random task graphs')
	#parser.add_argument("--conf",required=True,type=str,help='yaml file specifying task_dag generation parameters')
	parser.add_argument('--conf', nargs='?', const=conf_file, type=str,help='yaml file specifying task_dag generation parameters')
	return parser.parse_args()

def random_list(depth,total_num,width_min,width_max):
    list_t = []
    if total_num>= (depth-2)*width_min+2 and total_num<=(depth-2)*width_max+2:
        list_t.append(1)
        for i in range(depth-2):
            list_t.append(2)
        list_t.append(1)
        #print(list_t)
        for i in range(total_num-sum(list_t)+2):
            while True:
                tmp = random.randint(1,len(list_t)-2)
                if list_t[tmp]+1 > 4:
                    pass
                else:
                    list_t[tmp] = list_t[tmp] + 1
                    break
    else:
        list_t.append(1)
        num_tries = 30
        for num in range(0,num_tries):
            t = random.randint(width_min,width_max)
            # print('-------')
            # print(t)
            # print(total_num)
            # print(list_t)
            a = sum(list_t)-1+t
            # print(a)
            # print(width_min)
            b = total_num -(sum(list_t)-1)
            # print(b)
            # print(width_max)

            if (sum(list_t)-1+t)<total_num:
                list_t.append(t)
            elif  total_num -(sum(list_t)-1) >=width_min and total_num -(sum(list_t)-1)<=width_max :
                list_t.append(total_num -(sum(list_t)-1))
                break
            else:
                print('something wrong')
                pass
        list_t.append(1)
    return list_t

def gen_task_nodes(depth,total_num,width_min,width_max):
	#num_levels = depth+2		# 2: 1 for entry task, 1 for exit task
	num_list = random_list(depth+2,total_num,width_min,width_max)
	print(num_list)
	num_levels = len(num_list)
	num_nodes_per_level = np.array(num_list)
	#num_nodes_per_level = np.array([random.randint(width_min,width_max) for i in range(num_levels)])
	num_nodes_per_level[0] = 1.
	num_nodes_per_level[-1] = 1.
	num_nodes = num_nodes_per_level.sum()
	level_per_task = reduce(lambda a,b:a+b, [[enum]*val for enum,val in enumerate(num_nodes_per_level)],[])
	#e.g. [0,1,2,2,3,3,3,3,4]
	level_per_task = {i:level_per_task[i] for i in range(num_nodes)}
	#level_per_task in the format of {task_i: level_of_i}
	task_per_level = {i:[] for i in range(num_levels)}
	for ti,li in level_per_task.items():
		task_per_level[li] += [ti]
		# task_per_level in the format of {level_i:[tasks_in_level_i]}
	return task_per_level, level_per_task

def gen_task_links(deg_mu,deg_sigma,task_per_level,level_per_task,delta_lvl=2):
	num_tasks = len(level_per_task)
	num_level = len(task_per_level)
	neighs_top_down = {t:np.array([]) for t in range(num_tasks)}
	neighs_down_top = {t:np.array([]) for t in range(num_tasks)}
	deg = np.random.normal(deg_mu,deg_sigma,num_tasks)
	deg2 = (deg/2.).astype(np.int)
	deg2 = np.clip(deg2,1,20)
	#add edges from top to down with deg2, then bottom-up with deg2
	edges = []
	# ---- top-down ----
	for ti in range(num_tasks):
		if level_per_task[ti] == num_level-1:	# exit task is a sink
			continue
		ti_lvl = level_per_task[ti]
		child_pool = []
		for li,tli in task_per_level.items():
			if li <= ti_lvl or li > ti_lvl+delta_lvl:
				continue
			child_pool += tli
		neighs_top_down[ti] = np.random.choice(child_pool,min(deg2[ti],len(child_pool)),replace=False)
		edges += [(str(ti),str(ci)) for ci in neighs_top_down[ti]]
	# ---- down-top ----
	for ti in reversed(range(num_tasks)):
		if level_per_task[ti] == 0:
			continue
		ti_lvl = level_per_task[ti]
		child_pool = []
		for li,tli in task_per_level.items():
			if li >= ti_lvl or li < ti_lvl-delta_lvl:
				continue
			child_pool += tli
		neighs_down_top[ti] = np.random.choice(child_pool,min(deg2[ti],len(child_pool)),replace=False)
		edges += [(str(ci),str(ti)) for ci in neighs_down_top[ti]]
	return list(set(edges)),neighs_top_down,neighs_down_top

def gen_attr(tasks,edges,ccr,comp_mu,comp_sigma,link_comm_sigma):
	task_comp = np.clip(np.random.normal(comp_mu,comp_sigma,len(tasks)), EPSILON, comp_mu+10*comp_sigma)
	link_comm = np.zeros(len(edges))
	link_comm_mu = comp_mu * ccr
	#link_comm is the data transmitted on links, comp is the computation workload. They both follow normal distribution. ccr is a constant
	link_comm = np.clip(np.random.normal(link_comm_mu,link_comm_sigma*link_comm_mu,len(edges)),EPSILON, link_comm_mu+10*link_comm_sigma*link_comm_mu)
	return task_comp,link_comm

def plot_dag(dag,dag_path_plot):
	pos = graphviz_layout(dag,prog='dot')
	node_labels = {n:'{}-{:3.1f}'.format(n,d['comp']) for n,d in dag.nodes(data=True)}
	edge_labels = {(e1,e2):'{:4.2f}'.format(d['data']) for e1,e2,d in dag.edges(data=True)}
	plt.clf()
	nx.draw(dag,pos=pos,labels=node_labels,font_size=8)
	nx.draw_networkx_edge_labels(dag,pos,edge_labels=edge_labels,label_pos=0.75,font_size=6)
	plt.savefig(dag_path_plot)



def prepare_task_dag(config_yml,dag_path_plot):
	with open(config_yml) as f_config:
		config = yaml.load(f_config)
	#--- generate task graph ---

	task_per_level,level_per_task = gen_task_nodes(config['depth'],config['total_num'],config['width_min'],config['width_max'])
	edges,adj_list_top_down,adj_list_down_top = gen_task_links(config['deg_mu'],config['deg_sigma'],task_per_level,level_per_task)
	task_comp,link_comm = gen_attr(np.arange(len(level_per_task)),edges,config['ccr'],config['comp_mu'],config['comp_sigma'],config['link_comm_sigma'])
	edge_d = [(e[0],e[1],{'data':link_comm[i]}) for i,e in enumerate(edges)]
	dag = nx.DiGraph()
	dag.add_edges_from(edge_d)
	for i,t in enumerate(task_comp):
		dag.node[str(i)]['comp'] = t
		##NOTE: actually it should not be called 'comp', but 'computation data amount'!!!
	if dag_path_plot is not None:
		plot_dag(dag,dag_path_plot)
	#print(dag.graph)
	
	return dag

#Generate configuration.txt
def generate_config(dag,app_path):
	f = open(app_path, 'w')
	total_node = len(dag.nodes())
	f.write(str(total_node) + "\n")
	task_dict = {}
	for j in dag.nodes():
		task_dict[j] = 0
	task_dict['0'] = 1
	for e0, e1 in dag.edges():
		if e1 in task_dict.keys():
			task_dict[e1] += 1

	data = dict()
	for i in dag.nodes():
		if i not in data.keys():
			data[i] = ""
		data[i] += str(task_dict[i]) + " "
		#data[i] += "true " ## send single output to all the children tasks
		data[i] += "false " ## send all children output to all the children tasks correspondingly
		for e0, e1 in dag.edges():
			if i == e0:
				data[i] +="task" + e1 + " "
		if int(i) == total_node - 1:
			data[i] += "home"
	for i in range(len(data)):
		f.write("task" + str(i) + " ")
		f.write(data[str(i)])
		f.write('\n')

	f.close()

def generate_communication(dag,app_path):
	f = open(app_path,'w')
	for i in dag.nodes():
		f.write("task"+i+ " ")
		for e0,e1,d in dag.edges(data=True):
			if i == e0:
				f.write('task'+e1+'-'+str(round(d['data'],2))+' ') #KB
		f.write("\n")
	f.close()


#Generate dummy scripts
def generate_scripts(dag,config_path,script_path,app_path,sample_path):
    

    print('------ Read input parameters from DAG  -------------')
    sorted_list = sorted(dag.nodes(data=True), key=lambda x: x[0], reverse=False)
    f = open(config_path, 'r')
    total_line = f.readlines()
    #dev = total_line[0]
    del total_line[0]
    f.close()


    for i in range(len(total_line)):
        tmp = total_line[i].split(' ')
        filename = tmp[0] + ".c"
        num = int(tmp[1])
        file_path = script_path + filename

        f = open(file_path, 'w')
        f.write("#include<stdio.h>\n")
        f.write("#include<stdlib.h>\n")
        f.write("#include<string.h>\n")
        f.write("#include<math.h>\n")
        f.write("#include<time.h>\n")
        f.write("#include<sys/stat.h>\n")
        f.write("#include<sys/types.h>\n")
        f.write("#define PATH_MAX 128\n")
        f.write("#include<libgen.h>\n")
        f.write("#define LSIZ 128\n")
        f.write("#define RSIZ 20\n")


        
        f.write('char* filename = "' + tmp[0] + '";\n')
        f.write("char file_path[RSIZ][LSIZ];\n")
        f.write("int execution_time;\n")
        f.write("char task_name[LSIZ];\n")
        f.write("char output1_list[LSIZ];\n")
        f.write("char input_file[LSIZ];\n")
        f.write("char output_name[LSIZ];\n")
        f.write("char new_path[LSIZ];\n")
        f.write("char bash_script[LSIZ];\n")
        f.write("int loop_var;\n")
        f.write("char* ptr;\n")
        f.write("char* base;\n")
        f.write("char rand_file[50];\n")
        f.write("char output_path[RSIZ][LSIZ];\n")
        f.write("char output_list[RSIZ][LSIZ];\n")
        f.write("char file_size[RSIZ][LSIZ];\n")
        f.write("char new_file[LSIZ];\n")

        f.write("char line[RSIZ][LSIZ];\n")
        f.write("char temp_comm[RSIZ][LSIZ];\n")
        f.write("char src[RSIZ][LSIZ];\n")
        f.write("char comm[RSIZ][LSIZ];\n")
        f.write("char dest_temp[RSIZ][LSIZ];\n")
        f.write("char dest[RSIZ][LSIZ];\n")
        f.write("char sizes[RSIZ][LSIZ];\n")
        f.write("char str[LSIZ];\n")
        f.write("char* fname;\n")
        f.write("FILE* fptr = NULL;\n")
        f.write("int i=0;\n")
        f.write("int k=0;\n")
        f.write("int flag=0;\n")
        f.write("static int tot=0;\n")
        f.write("int idx1=0;\n")
        f.write("int idx2=0;\n")

        f.write("\n")
        f.write("char* task(char* filename,char* pathin,char* pathout)\n")
        f.write("{ \n")
        f.write("\texecution_time = rand() % 10;\n")
        f.write('\tprintf("%d\\n",execution_time);\n')
        f.write("\ttime_t timeout = 0;\n")
        f.write("\ttimeout = time(&timeout) + execution_time;\n")

        f.write("\twhile (time(&timeout) < timeout)\n")
        f.write("{ \n")
        f.write("\t\t1+1 ;\n")
        f.write("} \n")

        f.write('\tchar *path = __FILE__;\n')
        f.write("\tchar *path_cpy = strdup(path);\n")
        f.write("\tbase=basename(path_cpy);\n")
        f.write("\tmemcpy(task_name,base,LSIZ);\n")

        f.write('\tprintf("-------------------------\\n");\n')
        f.write('\tprintf("%s\\n",task_name);\n')
        f.write('\tprintf("%s\\n",filename);\n')
        #f.write("\tprintf(pathin)\n")
        #f.write("\tprint(pathout)\n")
        f.write('\tprintf("-------------------------\\n");\n')
        f.write("\tstrcpy(output1_list,filename);\n")
        f.write('\tstrcat(output1_list,"_");\n')
        f.write("\tstrcat(output1_list,task_name);\n")
        f.write('\tstrcat(output1_list,".txt") ;\n')
        f.write("\tstrcpy(input_file,filename);\n")
        f.write('\tprintf("%s\\n",output1_list);\n')
        f.write('\tprintf("%s\\n",input_file);\n')

        f.write('\tstrcpy(output_name,input_file);\n')
        f.write('\tstrcat(output_name,"_");')
        f.write('\tstrcat(output_name,task_name);')
        f.write('\tprintf("%s\\n",output_name);\n')
        f.write('\tprintf("---------------------------\\n");\n')

        f.write("\tchar actualpath [PATH_MAX+1];\n")
        f.write('\tchar *file_comm = "communication.txt";\n')
        f.write("\tchar *ptr;\n")
        f.write("\tptr=realpath(file_comm,actualpath);\n")
        f.write('\tprintf("%s\\n",actualpath);\n')

        f.write('\tfname="communication.txt";\n')
        f.write('\tfptr=fopen(fname, "r" );\n')
        #f.write("\tfclose(f);\n")

        f.write("\twhile(fgets(line[i], LSIZ, fptr)) \n")
        f.write("{ \n")
        f.write("\t\tline[i][strlen(line[i]) - 1] = '\\0'; \n")
        f.write("\t\tstrcpy(temp_comm[i], line[i]);\n")
        f.write("\t\t i++;\n")
        f.write("} \n")
        f.write("\ttot=i;\n")
    
        f.write("\tfor(i=0; i < tot; ++i)\n")
        f.write("{ \n")
        f.write("\t\tstrncpy(src[i], temp_comm[i], 6);\n")
        f.write("\t\tchar *tmp = strchr(temp_comm[i], ' ');\n")
        f.write("\t\tif(tmp != NULL)\n")
        f.write("\t{ \n")
        f.write("\t\t\tstrcpy(dest_temp[i],tmp+1);\n")
        f.write("\t\t\tif(strlen(dest_temp[i])>0)\n")
        f.write("{ \n")
        f.write("\t\t\t\tstrcpy(comm[i] ,dest_temp[i]);\n")
        f.write("} \n")
        f.write("\t} \n")                                                      #IF TMP CLOSE
        f.write("} \n")                                                        #FOR I CLOSE                   
        f.write("\tfor(i=0; i < tot; ++i)\n")                                  #SECOND FOR I
        f.write("{ \n")
        f.write('\tif (strncmp(task_name,src[i],5)==0);\n')
        f.write("{ \n")
        f.write('\tstrcpy(str, comm[i]);\n')
        f.write("\tchar* pch; \n")
        f.write("\tint j=0;\n")
        f.write('\tpch = strtok (str," -");\n')
        f.write("\twhile (pch != NULL)\n")
        f.write("{ \n")
        f.write("\t\tif(j%2 == 0)\n")
        f.write("{ \n")
        f.write("\t\t\t\tstrcpy(dest[idx1], pch);\n")
        f.write("\t\t\t\tidx1++;\n")
        f.write("} \n")
        f.write("else\n")
        f.write("{  \n")
        f.write("\t\t\t\tstrcpy(sizes[idx2], pch);\n")
        f.write("\t\t\t\tidx2++;\n")
        f.write("} \n")
        f.write("\tj++; \n")
        f.write('\tpch = strtok (NULL, " -");\n')
        f.write("} \n")                                          #WHILE PCH close


        f.write("\tfor(int k = 0;k < idx1; k++)\n")
        f.write("{  \n")
        f.write('\t\t\t\tprintf("The neighor is:\\n");\n')
        f.write('\t\t\t\tprintf("%s\\n",dest[k]);\n')
        f.write('\t\t\t\tprintf("The IDX  is:\\n");\n')
        f.write('\t\t\t\tprintf("%d\\n",k);\n')
        f.write('\t\t\t\tstrcpy(new_file,output_name);\n')
        f.write('\t\t\t\tstrcat(new_file,"_");\n')
        f.write('\t\t\t\tstrcat(new_file,dest[k]);\n')
        f.write('\t\t\t\tstrcpy(output_list[k],new_file);\n')
        f.write('\t\t\t\tstrcpy(file_size[k],sizes[k]);\n')
        f.write('\t\t\t\tchar new_path[LSIZ];\n')
        f.write('\t\t\t\tstrcpy(new_path,pathout);\n')
        f.write('\t\t\t\tstrcpy(new_path,new_file);\n')
        f.write('\t\t\t\tstrcpy(output_path[k],new_path);\n')
        f.write('\t\t\t\tprintf("NEW PATH IS %s \\n",new_path);\n')
        f.write('\t\t\t\tstrcat(bash_script,"/centralized_scheduler/generate_random_files.sh");\n')
        f.write('\t\t\t\tstrcat(bash_script,new_path);\n')
        f.write('\t\t\t\tstrcat(bash_script,sizes[k]);\n')
        f.write('\t\t\t\tsystem(bash_script);\n')
        f.write("} \n")
        f.write("} \n")                                                     #END OF STRNCMP
        f.write('\t\telse if((strncmp(task_name,src[i],5)!=0))\n')
        f.write('\t\tflag++;\n')
        f.write("} \n")                                                     #END OF FOR



        f.write("\tif(flag == tot)\n")
        f.write("{  \n")
        f.write('\t\t\t\tint k=0;\n')
        f.write('\t\t\t\tstrcpy(new_file,output_name);\n')
        f.write('\t\t\t\tstrcat(new_file,"_");\n')
        f.write('\t\t\t\tstrcat(new_file,task_name);\n')
        f.write('\t\t\t\tchar new_path[LSIZ]; \n ')
        f.write('\t\t\t\tstrcpy(new_path,pathout);\n')
        f.write('\t\t\t\tstrcat(new_path,new_file);\n')
        f.write('\t\t\t\tstrcpy(output_path[k],new_path);\n')
        f.write('\t\t\t\tint x= (rand() + 1);\n')
        f.write('\t\t\t\tsprintf(rand_file,"%d",x);\n')
        f.write('\t\t\t\tstrcpy(new_file,output_name);\n')
        f.write('\t\t\t\tstrcat(bash_script,"/centralized_scheduler/generate_random_files.sh"); \n')
        f.write('\t\t\t\tstrcat(bash_script,new_path);\n')
        f.write('\t\t\t\tstrcat(bash_script,rand_file);\n')
        f.write('\t\t\t\tsystem(bash_script);\n')
        f.write("}  \n")

        f.write("\treturn (char *)output_path; \n")
        f.write("}  \n")

        f.write("\tchar* main()\n")
        f.write("{ \n")
        f.write('\tchar filelist[128] = "1botnet.ipsum"; \n')
        f.write('\tchar* outpath="/Desktop/testshraddha";\n')
        f.write('\tchar* final= task(filelist, outpath, outpath);\n')
        f.write('\tfor(int v=0;v < tot;v++)\n')
        f.write('\t{\n')
        f.write('\t\t\tprintf("OUTPUT PATH: %s\\n",output_path[v]);\n')
        f.write('\t}\n')
        f.write("}\n")                                          #END OF MAIN
                                                         
        

        

    shutil.copy('jupiterapp/app_config.ini',app_path)
    shutil.copy('jupiterapp/input_node.txt',app_path) #for WAVE
    shutil.copy('jupiterapp/generate_random_files.sh',script_path)
    os.mkdir(sample_path)
    shutil.copy('jupiterapp/1botnet.ipsum',sample_path)
    shutil.copy('jupiterapp/2botnet.ipsum',sample_path)



def generate_nameconvert(dag,app_path):

	f = open(app_path,'w')
	for i in dag.nodes():
		s = "task"+i+ " botnet botnet\n"
		f.write(s)
	f.write('input botnet botnet\n')
	f.close()


def generate_json(dag,app_path):
	f = open(app_path, 'w')
	print(app_path)
	taskname_map = dict()
	exec_profiler = dict()
	for i in dag.nodes():
		taskname_map["task"+i] = ["task"+i,True]
		exec_profiler["task"+i] = True
	final_json = {"taskname_map":taskname_map,"exec_profiler":exec_profiler}
	f.write(json.dumps(final_json,indent = 2))
	f.close()

def generate_dummy_app(dummy_conf_file,dummy_app_path):
	args = parse_args(dummy_conf_file)
	print(args)
	dummy_dag_plot = dummy_app_path + 'dag.png'
	dummy_config_path = dummy_app_path + 'configuration.txt'
	dummy_script_path = dummy_app_path + 'scripts/'
	dummy_json_path = dummy_script_path + 'config.json' 
	dummy_comm_path = dummy_script_path + 'communication.txt' 
	dummy_sample_path = dummy_app_path + 'sample_input/'
	dummy_name_path = dummy_app_path + "name_convert.txt"
	
	if os.path.isdir(dummy_app_path):
		shutil.rmtree(dummy_app_path)
	os.mkdir(dummy_app_path)

	print('Create dummy_app folder, generate DAG and plot corresponding dag.png')
	dag = prepare_task_dag(args.conf,dummy_dag_plot)

	print('Generate configuration.txt')
	generate_config(dag,dummy_config_path)

	

	os.mkdir(dummy_script_path)
	print('Generate dummy application scripts')
	generate_scripts(dag, dummy_config_path,dummy_script_path,dummy_app_path,dummy_sample_path)

	print('Generate communication.txt')
	generate_communication(dag,dummy_comm_path)

	print('Generate config.json file')
	generate_json(dag,dummy_json_path)

	print('Generate name_convert.txt')
	generate_nameconvert(dag,dummy_name_path)
	print('The dummy application is generated successfully!')

def generate_multiple_apps(dummy_app_root,dummy_conf_root,N,start,max_conf_num):
	for i in range(start,N+1):
		print(i)
		dummy_app_path= '%sdummy_app%d/'%(dummy_app_root,i)
		rand_config = random.randint(1,max_conf_num)
		dummy_conf_path = '%stask_config%d.yml'%(dummy_conf_root,rand_config)
		generate_dummy_app(dummy_conf_path, dummy_app_path)

if __name__ == '__main__':
	# generate multiple apps
	# dummy_app_root = 'dummy_app_list/'
	# dummy_conf_root = 'dummy_task_config/'
	# N = 100 # number of dags
	# start = 1
	# max_conf_num = 5
	# generate_multiple_apps(dummy_app_root,dummy_conf_root,N,start,max_conf_num)
	dummy_conf_file = 'task_config_100.yml'
	dummy_app_path = 'dummy_app_100/'
	generate_dummy_app(dummy_conf_file,dummy_app_path)
