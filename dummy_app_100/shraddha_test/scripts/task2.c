#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<sys/stat.h>
#include<sys/types.h>
#define PATH_MAX 128\n#include<libgen.h>
int tot=i;
char filename[128][128];
char file_path[128][128];
int execution_time;
char* task_name;
char output1_list[128];
char input_file[128];
FILE* f;
char comm [128][128];
char *src[128][20];
char* new_path[128];
char bash_script[128];
int loop_var;
char* ptr;

void task(filename,pathin,pathout)
{ 
	execution_time = rand() % 10;
	printf("%d",execution_time);
	time_t timeout;
	timeout = time(&timeout) + execution_time;
	while (time.time() < timeout)
{ 
		1+1 ;
} 
	char *path = '__FILE__';
	char *path_cpy = strdup(path);
	task_name=basename(path_cpy);
	printf("-------------------------");
	printf("%s",task_name);
	printf("%s",filename);
	printf("-------------------------" );
		strcpy(output1_list,&filename);
		strcat(output1_list,"."");
		strcat(output1_list,&task_name[0]);
		strcat(output1_list,".txt") ;
		strcpy(input_file, &filename);
	printf("%s",output1_list);
	printf("%s",input_file);
		char output_name[128]=strcat(input_file,"_");
		strcat(output_name,output_name);	print("%s",output_name);
	printf("---------------------------");
	char actualpath [PATH_MAX+1];
	char *file_comm [PATH_MAX+1]= "communication.txt";
	char *ptr;;
	ptr=realpath(file_comm,actualpath);
	printf("%s",ptr);
	FILE* f;
	fopen(/centralized_scheduler/communication.txt, "r" );
	fclose(f);
	char temp_comm [128][128];
	char comm [128][128];
	memset(temp_comm, " ", sizeof(temp_comm));
	memset(comm, ' ', sizeof(comm));;
	char line_read [128];
	memset(line_read, ' ', sizeof(line_read));
	int k = 0;
	while( fgets( line_read, sizeof(line_read), f ) != NULL )
{ 
	strcpy(temp_comm[k], line_read);
	k++;
} 
	int sum_k = k;
	char *src[128][20];
	char *dest_temp[128][20];
	for(int loop_var=0; loop_var < sum_k; loop_var++)
{ 
	strncpy(src[loop_var], temp_comm[loop_var], 6);
	char *tmp = strchr(temp_comm[loop_var], ' ');
	if(tmp != NULL)
	char *dest_temp[loop_var]=tmp+1;
	if(strlen(dest_temp[loop_var])>0)
	comm[loop_var] = dest_temp[loop_var];
	printf("---------------------------");
	printf("%s",comm[loop_var]);
	printf("%s",task_name);
	struct stat pathout_check;
		stat(pathout, &pathout_check);
		if (!(S_ISDIR(pathout_check.st_mode)));
		mkdir(pathout);
	char output_path=[];
	if strcmp(task_name,src[loop_var])
{ 
	printf("%s",comm[loop_var]);
	char output_list[128][128];
	char file_size[128][128];
		char new_file[128];
	printf("The neighor is:");
	printf("%s",comm[loop_var]);
	printf("The IDX  is:");
	printf("%s",src[loop_var]);
	strcpy(new_file,output_name);
	strcat(new_file,"_");
	strcat(new_file,src[loop_var]);
	output_list[j]=new_file;
	file_size[j]=comm[j];
	char* new_path[128]; 
	memset(new_path, ' ', sizeof(new_path));
	strcpy(new_path,pathout);
	strcat(new_path,new_file);
	strcpy(output_path,new_path);
	printf("%s",new_path);
	 bash_script="/centralized_scheduler/generate_random_files.sh"; 
	strcat(bash_script,new_path);
	strcat(bash_script,file_size);
	system(bash_script); 
} 
	else if !strcmp(task_name,src[loop_var]) 
	 	 {  
	strcpy(new_file,output_name); 
	strcat(new_file,"_"); 
	strcat(new_file,task_name); 
	strcpy(output_path,new_path);
	 file_size=itoa(rand() + 1);
	strcpy(new_file,output_name); 
	 bash_script="/centralized_scheduler/generate_random_files.sh"; 
	strcat(bash_script,new_path);
	strcat(bash_script,file_size);
	system(bash_script); 
	 	 } 
	return output_path; 
	 	 } 

	 int main() 
	 	 {  
	char filelist[128] = "1botnet.ipsum"; 
	char output_returned[128] ; 
	char outpath[128]=dirname(__FILE__); 
	strcat(outpath,"sample_input/"); 
 	output_returned = task(filelist, outpath, outpath); 
	return 0; 
	 	 } 
