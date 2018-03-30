#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <sys/time.h>
#define  threadnum 5;
#define  totalit 1000000000;
int isint(char moji){
    return moji >= '0' && moji <= '9';
}

void time_out(FILE* fp, int* newtime,int* oldtime, int nowtime,
        int* newtimeu, int*oldtimeu, int nowtimeu, int threads, int type){
    *newtime = nowtime;
    *newtimeu = nowtimeu;
    fprintf(fp, "%f %d %d\n",  *newtime-(double)*oldtime+(*newtimeu-*oldtimeu)*1e-6, threads, type);
    *oldtime = *newtime;
    *oldtimeu = *newtimeu;
}

void get_breaks(int* breaks, int threads, int imax){
    int i;
    for(i=0; i!=threads; i++)
        breaks[i] = (i*imax)/threads;
    breaks[threads] = imax;
}

int main(int argc, char* argv[]){
    FILE *fp;
    fp = fopen("pi_results.data","a");
    int type =0, threads = threadnum;
    int myid, numprocs, sum;
    double pi;
    int i = 0;
    int imax = 0, temp_int;
    struct timeval s;
    int* breaks;
    int new_time, new_timeu;
    gettimeofday(&s, NULL);
    int old_time=s.tv_sec, old_timeu=s.tv_usec;
    if(argc != 1){
        threads = atoi(argv[2]);
        imax = atoi(argv[1]);
    }else imax = 1e5;
    breaks =  (int*)malloc(sizeof(int)*(threads+1));

    get_breaks(breaks, threads, imax);

    gettimeofday(&s, NULL);
    time_out( fp,&new_time, &old_time, s.tv_sec, &new_timeu, &old_timeu, s.tv_usec, threads, 0);
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid==0){
        gettimeofday(&s, NULL);
        time_out( fp,&new_time, &old_time, s.tv_sec, &new_timeu, &old_timeu, s.tv_usec,threads, 1);
    }
    double x,y;
    int count=0;
    srand(myid);
    int j = 0;
    for(j=breaks[myid]; j!=breaks[myid+1]; j++){
        x = rand()/((double)RAND_MAX+1);
        y = rand()/((double)RAND_MAX+1);
        if(x*x + y*y < 1)
            count++;
    }
    if(myid==0){
        gettimeofday(&s, NULL);
        time_out( fp,&new_time, &old_time, s.tv_sec, &new_timeu, &old_timeu, s.tv_usec, threads, 2);
    }
    MPI_Reduce(&count, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(myid==0){
        gettimeofday(&s, NULL);
        time_out( fp,&new_time, &old_time, s.tv_sec, &new_timeu, &old_timeu, s.tv_usec, threads, 3);
    }
    if(myid==0){
        pi=4*(double)sum/imax;
        printf("i=%d pi=%f\n",imax, pi);

    }
    MPI_Finalize(); 
    if(myid==0){
        gettimeofday(&s, NULL);
        time_out( fp,&new_time, &old_time, s.tv_sec, &new_timeu, &old_timeu, s.tv_usec,threads, 4);
    }
    fclose(fp); 
    return 0;
}
