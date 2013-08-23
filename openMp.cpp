#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>
#include<float.h>
#define NO_THDS 2
typedef struct {
  long int cabinetNum;
  long int id;
  double *sub_score;
  } doc;
  
 struct  d_ptr{
  doc *dct;
  struct d_ptr *next;
};
typedef struct d_ptr doc_ptr;
typedef struct {
  double *docs_sub_avg;
  double *docs_sub_sum;
  long int total_doc_cnt;
  long int *doc_cnt;
  doc_ptr **doc_ptr_arr;
  doc_ptr **last;
} cabinet;

int main(int argc ,char *argv[]){
  
  omp_set_num_threads(2);
	// Command line arguments
  double start_time=omp_get_wtime();
  long int num_Cabinets;
  long int num_Subjects;
  long int num_docs;
  FILE *in=fopen(argv[1],"r");;
   
    if(argc ==3){
    num_Cabinets=strtol(argv[2],NULL,10);
  //  printf("Number of cabinets=%ld \n",num_Cabinets);
    fscanf(in,"%*ld %ld %ld",&num_docs,&num_Subjects);
   // printf("Cabinets =%ld, Docs =%ld , sub=%ld",num_Cabinets,num_docs,num_Subjects);
  }
  
  if(argc==2){
    fscanf(in,"%ld %ld %ld",&num_Cabinets,&num_docs,&num_Subjects);
  //  printf("Cabinets =%ld, Docs =%ld , sub=%ld\n",num_Cabinets,num_docs,num_Subjects);
  }
  
  cabinet *cab_arr=(cabinet *)calloc(num_Cabinets,sizeof(cabinet));
  doc *doc_arr = (doc *)malloc(sizeof(doc)*num_docs);
  
  /*Initialization of documents with their scores and their cabinet no. **/
  
  int i=0,j=0,k;
  
  for(i=0;i<num_docs;i++){
    (doc_arr + i)->sub_score = (double *)malloc(sizeof(double)*num_Subjects);
    fscanf(in,"%ld",&((doc_arr+i)->id));
    k=i%num_Cabinets;
    (doc_arr+i)->cabinetNum = k;
    for(j=0;j<num_Subjects;j++){
      fscanf(in,"%lf",(doc_arr+i)->sub_score+j);
     // printf("For doc_id=%ld , cabinet=%ld sub_id=%d ,score=%lf \n",(doc_arr+i)->id ,(doc_arr+i)->cabinetNum ,j,*((doc_arr+i)->sub_score+j));    
    }
      
  }
  
  /*Allocating memory in cabinets */
 #pragma omp parallel for private(i)
  for(i=0;i<num_Cabinets;i++){
    printf("No of threads is %d\n",omp_get_num_threads());
    int no_thds= NO_THDS;
    cabinet *cbnt;
    cbnt =  (cab_arr+i);
    cbnt->docs_sub_avg = (double *)calloc(num_Subjects,sizeof(double));
    cbnt->docs_sub_sum = (double *)calloc(num_Subjects,sizeof(double));
    cbnt->doc_ptr_arr =(doc_ptr **)calloc(no_thds,sizeof(doc_ptr *));
    cbnt->last =(doc_ptr **)malloc(sizeof(doc_ptr *)*no_thds);
    cbnt->doc_cnt =(long int*)calloc(no_thds,sizeof(long int));
  }
  
  /* Assigning documents to respective cabinets*/
  #pragma omp parallel for private(i,j)
  for(i=0;i<num_docs;i++){
    //int no_thds= omp_get_num_threads();
    int tid = omp_get_thread_num();
    //printf("No of threads is %d\n",omp_get_num_threads());
    doc_ptr *dptr;
    cabinet *cbnt;long int l;
    int doc_cab_id = (doc_arr+i)->cabinetNum;
     cbnt =  (cab_arr+doc_cab_id);
     l= (*(cbnt->doc_cnt+tid));
     //printf("Length =%ld\n",l);
     if(l==0){
      
    *(cbnt->doc_ptr_arr+tid) =(doc_ptr *)malloc(sizeof(doc_ptr)*1);
    dptr = *(cbnt->doc_ptr_arr+tid);
     (dptr + 0)->dct = doc_arr+i;
     *(cbnt->last+tid) =  (dptr + 0);
     *(cbnt->doc_cnt+tid) = 1;
    
    }
    else{
      (*(cbnt->last+tid))->next = (doc_ptr *)malloc(sizeof(doc_ptr)*1);
      dptr = (*(cbnt->last+tid))->next;
      dptr->dct=doc_arr+i;
      (*(cbnt->last+tid))= dptr;
     *(cbnt->doc_cnt+tid) = *(cbnt->doc_cnt+tid) +1;
    }
    
    
  }
  /* printing the details in the cabinet..can be commented.Parallel on outer loop */
 /* for(i=0;i<num_Cabinets;i++){
   
    doc_ptr **doc_ptr_arr=(cab_arr+i)->doc_ptr_arr; 
    doc_ptr *tdptr;
    for(j=0;j<NO_THDS;j++){
      tdptr = *(doc_ptr_arr+j);
     // printf("No of documents in Cabinet=%d adn threadId=%d are %ld\n",i,j, *((cab_arr+i)->doc_cnt+j));
      while(tdptr!=NULL){
	//printf("From cabninet=%d, thread_id=%d , doc_id=%ld\n",i,j,tdptr->dct->id);
	tdptr =tdptr->next;
      }
    }
    
  }*/
  
  /*  Minimizing the distances of documents to the cabinets */
  int flag=1;
  for(;flag==1;){
    flag=0;
    /*Compute the averages for each cabinet and free the contents of the cabinets as they are going to 
     be filled with new document pointers*/
    int i=0;cabinet *cbnt;
//#pragma omp parallel for private(i,cbnt)
    for(i=0;i<num_Cabinets;i++){
    cbnt = cab_arr+i;
    doc_ptr **doc_ptr_arr=(cab_arr+i)->doc_ptr_arr; 
    doc_ptr *tdptr,*td;
    doc *dct;int total_doc_cnt=0;
    int b=0;
    for(b=0;b<num_Subjects;b++){
      *(cbnt->docs_sub_sum+b) = 0.0f;
      *(cbnt->docs_sub_avg+b) = 0.0f;
    }
    for(j=0;j<NO_THDS;j++){
      int k=0;
      total_doc_cnt +=*(cbnt->doc_cnt+j);
      tdptr = *(doc_ptr_arr+j);
      while(tdptr!=NULL){
	dct = tdptr->dct;
	for(k=0;k<num_Subjects;k++){
	  *(cbnt->docs_sub_sum+k) +=*(dct->sub_score+k);
	}
	//printf("From cabninet=%d, thread_id=%d , doc_id=%ld\n",i,j,tdptr->dct->id);
	td=tdptr;
	tdptr =tdptr->next;
	free(td);
      }
      *(doc_ptr_arr+j) =NULL;
      *(cbnt->doc_cnt+j) = 0;
    }
    //printf("Total number of documents in cabinet=%d are %d\n",i,total_doc_cnt);
    int l;
    /* Finding averages of subjects for the cabinet */
     //printf("Cabinet_Id=%d ",i);
    for(l=0;l<num_Subjects;l++){
      *(cbnt->docs_sub_avg+l) = *(cbnt->docs_sub_sum+l)/total_doc_cnt;
      //printf(" sub_id=%d avg_score=%lf \n",l,*(cbnt->docs_sub_avg+l));
    }
    //printf("\n");
  }
 // printf("Execution is correct till here \n");
  /* For each document finding the cabinet where it will have minimal score */
#pragma omp parallel for private(i)
  for(i=0;i<num_docs;i++){
    int k,j;double min_score=DBL_MAX;double score,sc;cabinet *min_cab; long int min_cab_num;
    for(k=0;k<num_Cabinets;k++){
      score=0.0f;
      for(j=0;j<num_Subjects;j++){
	sc = (*((doc_arr+i)->sub_score+j) - *((cab_arr+k)->docs_sub_avg+j));
	score = score+(sc*sc);
      }
      if(sqrt(score) < min_score){
	min_score=sqrt(score);
	min_cab=cab_arr+k;
	min_cab_num=k;
	}
    }
    //printf("For Document=%d the cabinetNo=%ld\n",i,min_cab_num);
    if((doc_arr+i)->cabinetNum == min_cab_num){
      
      doc_ptr *dp;
      int tid = omp_get_thread_num();
     // int tid =0;
      if(*(min_cab->doc_ptr_arr+tid) == NULL){
	*(min_cab->doc_ptr_arr+tid)=(doc_ptr *)malloc(sizeof(doc_ptr)*1);
	 dp = *(min_cab->doc_ptr_arr+tid);
	 dp->dct=doc_arr+i;
	// dp->dct->cabinetNum=min_cab_num;
	 *(min_cab->last+tid) = dp;
	  (*(min_cab->last+tid))->next=NULL;
	 *(min_cab->doc_cnt+tid) +=1;
	}
      else{
	(*(min_cab->last+tid))->next = (doc_ptr *)malloc(sizeof(doc_ptr));
	(*(min_cab->last+tid))->next->dct=doc_arr+i;
	(*(min_cab->last+tid))->next->dct->cabinetNum=min_cab_num;
	 *(min_cab->doc_cnt+tid) +=1;
	  *(min_cab->last+tid) = (*(min_cab->last+tid))->next;
	 (*(min_cab->last+tid))->next=NULL;
      }
    }
    else{
      /*change the flag to 1
       * asssign this document pointer in its respective cabinet
       */
      (doc_arr+i)->cabinetNum = min_cab_num;
      flag=1;doc_ptr *dp;
      int tid = omp_get_thread_num();
      if(*(min_cab->doc_ptr_arr+tid) == NULL){
	*(min_cab->doc_ptr_arr+tid)=(doc_ptr *)malloc(sizeof(doc_ptr)*1);
	 dp = *(min_cab->doc_ptr_arr+tid);
	 dp->dct=doc_arr+i;
	 dp->dct->cabinetNum=min_cab_num;
	 *(min_cab->last+tid) = dp;
	  (*(min_cab->last+tid))->next=NULL;
	 *(min_cab->doc_cnt+tid) +=1;
	}
      else{
	(*(min_cab->last+tid))->next = (doc_ptr *)malloc(sizeof(doc_ptr));
	(*(min_cab->last+tid))->next->dct=doc_arr+i;
	(*(min_cab->last+tid))->next->dct->cabinetNum=min_cab_num;
	 *(min_cab->doc_cnt+tid) +=1;
	  *(min_cab->last+tid) = (*(min_cab->last+tid))->next;
	 (*(min_cab->last+tid))->next=NULL;
      }
    }
  }
  
  }
  printf("PRINTING THE RESULT\n");
  int m=0;
  FILE *fg=fopen("/home/pushparaj/pdc/xfile.out","w");
  if(fg==NULL){
    printf("output file is not valid\n");
  }
  for(m=0;m<num_docs;m++){
    int rs=fprintf(fg,"%d %ld\n",m,((doc_arr+m)->cabinetNum));
   // printf("%d %ld\n",m,((doc_arr+m)->cabinetNum));
  }
 // printf("error value for file %d\n",ferror(fg));
  double end_time = omp_get_wtime();
  
  printf("Time for execution is %lf \n",end_time-start_time);
  fclose(in);fclose(fg);
  //printf("No of threads is %d\n",omp_get_num_threads());
  
  
  /* Computing the averages in each cabinet and  Minimizing the distances of documents with respect to cabibets */
 /* for(;flag==1;){
    
    
    
  }*/
  
  /*int no_threads= omp_get_num_threads(); 
  int thread_id = omp_get_thread_num();

  printf("No_of threads = %d and thread id is =%d\n",no_threads,thread_id);*/
  
  }
  
  
  
  
  
  
  
  

 