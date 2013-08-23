
/*
 * docsmpi.c
 * testing version
 *  Created on: Nov 18, 2012
 *      Author: zellist (Tran Anh Phuong)
 *  Contiguous memory allocation (alloc_2d())
 *	Sending the whole chunk with MPI_iSEND()
 *	Freeing memory after use (free_2d())
 *	Output file name = "input file name" + "outx"
 *	Support huge file reading
 *	ubuntu - friendly
 *  time counting
 *
 *	Args: <input file name> <output file name> [cabinet]
 */



#include <string.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MASTER 0  /*taskID of the master's thread*/
#define ONEGB (128*1024*1024)  //currently support to 1GB of document data at each computer. Each computer suppose to be able to store >2x 1GB
#define SCALE 1		//change the value of SCALE in case of larger required memory, e.g: SCALE = 2 -> 2GB supported
#define NORMAL_FILE 0
#define HUGE_FILE 1

/*
 * the order of chunk,owner is reversed
 * task 0 will be in charge of the last chunk, while the taskid = (numtasks -1) will be in charge of the first chunk
 */
#define BLOCK_LOW(id,p,n)  ((p-1-id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(index,p,n) ((p-1)-(((p)*((index)+1)-1)/(n)))

/*Number of documents that can be stored in 1 chunk.*/
#define MAXCHUNKSIZE(num_subject) ((ONEGB*SCALE)/(num_subject))

/* distance calculating function*/
double distance(int size, double x[size],double y[size+1]){
	int i;
	double distance=0;
	for(i=0;i<size;i++){
		distance+=pow((x[i]-y[i]),2);
	}
	return distance;
}
/* contiguous memory allocating function*/
double **alloc_2d(int rows, int cols) {
    int i=0;
    double *data = (double *)calloc(rows*cols,sizeof(double));
    double **array= (double **)calloc(rows,sizeof(double*));
    for (i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}
/* memory freeing function*/
void free_2d(double **a){
		free(a[0]);
        free(a);
}


main(int argc, char **argv){


	/*MPI's related variables*/
		int numtasks,taskid,len;
		char hostname[MPI_MAX_PROCESSOR_NAME];
		MPI_Status status;
		MPI_Request myRequest1,myRequest2;

	/*Program's related variables*/
		int Number_Cabinets, Number_Documents, Number_Subjects;
		int chunksize;
		int i,j,k,temp;
		int dest;
		int repeatflag=1;
		int cabinet_num;
		int loop = 0;
		int position=0;
		/*int send;
		int sending_size;
		int sending_part,masterchunksize;*/
		int mode = 0;
		double start_time,computation_time;
		int doc_info[3];

	/*MPI Initialization*/
	MPI_Init(&argc, &argv);
	start_time = MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);  /*total number of tasks*/
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);  /*my TaskID*/
	MPI_Get_processor_name(hostname, &len);  /*my computer's name*/


	FILE *file;


	if(taskid == MASTER){     /*Only the master thread do this*/
			file = fopen ( argv[1], "r" );
			/*
			fscanf(file,"%d",&Number_Cabinets);
			fscanf(file,"%d",&Number_Documents);
			fscanf(file,"%d",&Number_Subjects);*/
			fscanf(file,"%d",&doc_info[0]);
			fscanf(file,"%d",&doc_info[1]);
			fscanf(file,"%d",&doc_info[2]);
			if(argv[3]!=NULL) {
				doc_info[0]=atoi(argv[3]);
			}

	}
		/*MASTER Broadcasts basic information about documents to other tasks*/
		MPI_Bcast(&doc_info[0],3,MPI_INT,MASTER,MPI_COMM_WORLD);


		Number_Cabinets=doc_info[0];
		Number_Documents=doc_info[1];
		Number_Subjects=doc_info[2];

		chunksize = BLOCK_SIZE(taskid,numtasks,Number_Documents);

		/*Preparing memory to store documents and send to others*/
		double**  doc = alloc_2d(chunksize,Number_Subjects);

		if(doc==NULL) {printf("Task %d is a bad boy. He doesn't have enough memory! \n",taskid);}

		/*choosing file reading mode*/
		if (BLOCK_SIZE(MASTER,numtasks,Number_Documents)<MAXCHUNKSIZE(Number_Subjects)){
			mode = NORMAL_FILE;
		}
		else {
			mode = HUGE_FILE;
			printf("Wow! it's a huge file\n", taskid);
		}

	if (taskid==MASTER){
		if (mode ==  NORMAL_FILE){
			/* The chunk size is smaller than [SCALE]GB  --> computer can store 2 chunks at the same time.*/
			double** tempdoc = alloc_2d(chunksize,Number_Subjects);
			if(tempdoc==NULL) {printf("Memory starving! \n");}
			dest=numtasks-1;
			while(dest>0){
				if((dest%2)==1){

					if(dest<(numtasks-2)){

						MPI_Wait(&myRequest1,&status);
					}
					int size=BLOCK_SIZE(dest,numtasks,Number_Documents);
					for(k=0;k<size;k++){

						fscanf(file,"%d",&temp);
						for (j=0;j<Number_Subjects;j++){
							fscanf(file, "%lf",&tempdoc[k][j]);
						}
					}

					MPI_Isend(&tempdoc[0][0],(Number_Subjects*size),MPI_DOUBLE,dest,dest,MPI_COMM_WORLD, &myRequest1);
					dest--;
				}
				else {

					if(dest<(numtasks-2)){

						MPI_Wait(&myRequest2,&status);
					}
					for(k=0;k<BLOCK_SIZE(dest,numtasks,Number_Documents);k++){
						fscanf(file,"%d",&temp);
						for (j=0;j<Number_Subjects;j++){
							fscanf(file, "%lf",&doc[k][j]);
						}
					}
					MPI_Isend(&doc[0][0],(Number_Subjects*BLOCK_SIZE(dest,numtasks,Number_Documents)),MPI_DOUBLE,dest,dest,MPI_COMM_WORLD, &myRequest2);
					dest--;
				}
			}
			if(numtasks>2){

				MPI_Wait(&myRequest2,&status);
			}



			for (i=0;i<chunksize;i++){
				fscanf(file,"%d",&temp);
				for (j=0;j<Number_Subjects;j++){
					fscanf(file, "%lf",&doc[i][j]);
				}
			}

			/* freeing memory*/
			if(numtasks>1){

				MPI_Wait(&myRequest1,&status);

						}
			if (tempdoc!=NULL) {free_2d(tempdoc);}
		}
		else {   /*HUGE_FILE mode*/
			/*actual chunk size is larger than [SCALE]GB. Computers are not supposed to be able to store 2 chunks at the same time*/
			dest=numtasks-1;
			while(dest>0){
				if(dest<(numtasks-1)){
					MPI_Wait(&myRequest1,&status);
				}
				for(k=0;k<(chunksize/2);k++){
					fscanf(file,"%d",&temp);
					for (j=0;j<Number_Subjects;j++){
						fscanf(file, "%lf",&doc[k][j]);
					}
				}
				MPI_Isend(&doc[0][0],(Number_Subjects*(chunksize/2)),MPI_DOUBLE,dest,1,MPI_COMM_WORLD, &myRequest1);
				if(dest<(numtasks-1)){
					MPI_Wait(&myRequest2,&status);
				}
				for(k=(chunksize/2);k<BLOCK_SIZE(dest,numtasks,Number_Documents);k++){
					fscanf(file,"%d",&temp);
					for (j=0;j<Number_Subjects;j++){
						fscanf(file, "%lf",&doc[k][j]);
					}
				}
				MPI_Isend(&doc[chunksize/2][0],(Number_Subjects*((BLOCK_SIZE(dest,numtasks,Number_Documents))-(chunksize/2))),MPI_DOUBLE,dest,2,MPI_COMM_WORLD, &myRequest2);
				dest--;
			}

			if(numtasks>1){
					MPI_Wait(&myRequest1,&status);
				}
			for (i=0;i<chunksize/2;i++){
				fscanf(file,"%d",&temp);
				for (j=0;j<Number_Subjects;j++){
					fscanf(file, "%lf",&doc[i][j]);
				}
			}
			if(numtasks>1){
				MPI_Wait(&myRequest2,&status);
			}
			for (i=chunksize/2;i<chunksize;i++){
					fscanf(file,"%d",&temp);
					for (j=0;j<Number_Subjects;j++){
						fscanf(file, "%lf",&doc[i][j]);
					}
				}
		}


		fclose(file); /*close the file*/
	}  /*end of master task*/

	else {   /* All the other tasks do this*/

		/*Slave task receives data from MASTER*/
		if(mode == NORMAL_FILE){
			MPI_Recv(&doc[0][0],(Number_Subjects*chunksize),MPI_DOUBLE,MASTER,taskid,MPI_COMM_WORLD,&status);
		}
		else {
			int masterchunksize = BLOCK_SIZE(MASTER,numtasks,Number_Documents);
			MPI_Recv(&doc[0][0],(Number_Subjects*(masterchunksize/2)),MPI_DOUBLE,MASTER,1,MPI_COMM_WORLD,&status);
			MPI_Recv(&doc[masterchunksize/2][0],(Number_Subjects*((BLOCK_SIZE(dest,numtasks,Number_Documents))-(masterchunksize/2))),MPI_DOUBLE,MASTER,2,MPI_COMM_WORLD,&status);
		}
	} /*end of slave task*/



	/*assigning the initial cabinet division*/
	int result[chunksize];
	position=BLOCK_LOW(taskid,numtasks,Number_Documents);

		for(i=0;i<chunksize;i++){
			result[i]=(position+i)%Number_Cabinets;
		}


		/*preparing memory for the cabinet*/

	double** cabinet=alloc_2d(Number_Cabinets,Number_Subjects+1);
	if(cabinet==NULL) {printf("No Cabinet memory for you today!\n");exit(0);}
	double** tempcabinet=alloc_2d(Number_Cabinets,Number_Subjects+1);
	if(tempcabinet==NULL) {printf("No memory for you today!\n");exit(0);}
	computation_time = MPI_Wtime();
		/*starting of the while*/
		while(repeatflag==1){
		  /*if(taskid==MASTER) {loop++;printf("%d\n",loop);}*/
			repeatflag=0;
			/* initialize cabinet */
			for(i=0;i<Number_Cabinets;i++){
					for(j=0;j<=Number_Subjects;j++){
						cabinet[i][j]=0;
					}
				}
			/*summing up the cabinet from chunk's data */
			for(i=0;i<chunksize;i++)
				{ 	cabinet_num=result[i];
					for(j=0;j<Number_Subjects;j++)
						{
						cabinet[cabinet_num][j]+=doc[i][j];
						}
					cabinet[cabinet_num][Number_Subjects]++;  //last element of cabinet store the number of Doc
				}


			/* Reducing values from cabinets */
			MPI_Allreduce(&cabinet[0][0],&tempcabinet[0][0],((Number_Cabinets)*(Number_Subjects+1)),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);


			/*Calculate the average*/
			for (i=0;i<Number_Cabinets;i++)
				{
					for(j=0;j<Number_Subjects;j++)
					{
								tempcabinet[i][j]=tempcabinet[i][j]/tempcabinet[i][Number_Subjects];
					}
				}




			for (i=0;i<chunksize;i++)  /*calculate the distance and decide whether we need to switch cabinets or not*/
				{
					double min;
					int switcher=0;double temp=0;
					min=distance(Number_Subjects,doc[i],tempcabinet[0]);
					for(j=1;j<Number_Cabinets;j++){
						temp=distance(Number_Subjects,doc[i],tempcabinet[j]);
						if (min>temp){
							min = temp;
							switcher=j;
						}
					}
					if (switcher != result[i]){
						result[i]=switcher;
						repeatflag=1;
					}

				}


			MPI_Allreduce(&repeatflag,&repeatflag,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

		}			/*end of while*/




	if (taskid==MASTER) printf("Computational time: %lf\n",(MPI_Wtime()-computation_time));

		/* finalresult[Number_Documents] to store at MASTER task */
		/*int finalresult[Number_Documents];*/

		/*int* displs = (int *)malloc(numtasks*sizeof(int));
		if (displs==NULL) printf("cannot allocate displs \n");
		int* rcounts = (int *)malloc(numtasks*sizeof(int));
		if (rcounts==NULL) printf("cannot allocate rcounts \n");*/
		int displs[numtasks];
		int rcounts[numtasks];
		/*root now has correct rcounts, using these we set displs[] so
					 * that data is placed contiguously (or concatenated) at receive end*/
		for (i=0; i<numtasks; ++i) {
			rcounts[i]= BLOCK_SIZE(i,numtasks,Number_Documents);

			}

		displs[numtasks-1]=0;
		for(i=1;i<numtasks;i++){
			displs[numtasks-1-i]=displs[numtasks-i]+rcounts[numtasks-i];

		}



		/*int* finalresult= (int *)malloc(Number_Documents*sizeof(int));
				if (finalresult==NULL){
					printf("too bad!\n");

				}*/


		int finalresult[Number_Documents];
		MPI_Gatherv(&result[0],chunksize,MPI_INT,&finalresult[0],rcounts, displs, MPI_INT, MASTER, MPI_COMM_WORLD);

			if(taskid==MASTER) {

				FILE *resultfile = fopen (argv[2], "w" );
					for (i=0;i<Number_Documents;i++){
							fprintf(resultfile,"%d %d\n",i,finalresult[i]);
						}
					fclose(resultfile);

					printf("Elapsed time: %lf\n",(MPI_Wtime()-start_time));

			}
		if(doc!=NULL)	{free_2d(doc);}
		if(tempcabinet!=NULL)	{free_2d(tempcabinet);}
		if(cabinet!=NULL)	{free_2d(cabinet);}
	MPI_Finalize();

}
/*
 * TODO: MPI + openMP
 * TODO:
 */
