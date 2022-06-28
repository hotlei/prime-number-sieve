/*  Copyright (C) Lei Yan */
#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <stdlib.h>

int main(argc,argv)
int argc;    char *argv[];
{
    int         myid, numprocs, resultlen;
    char        name[MPI_MAX_PROCESSOR_NAME];
    double      startwtime, endwtime;
    int         j,m,k,p,n,n64;
    register    t1,t2,i;
    static int  key[32],mask[32];
    static int  ikey[32],imask[32],tran[32];
    int         *arr, *seed;
    char        *little;
    int         num,root,root2;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(name, &resultlen); 

    startwtime = MPI_Wtime();
    num = 16000000;
    num = (num + numprocs * 64 - 1) / (numprocs * 64);
    num = num * numprocs * 64;
    n = num/numprocs;
    n64 = n / 64;
    root = sqrt(num) +1;
    root2 = sqrt(root) + 1;

    arr = (int*)malloc(n64 * sizeof(int));
    seed = (int*)malloc(root * sizeof(int) / 4);

    for(i=0;i<n64;i++)
        arr[i] = 0xffffffff;
    if(myid == 0) arr[0] = 0xfffffffe;
    key[0]=1;
    for(i=1;i<32;i++) key[i] = key[i-1] << 1;
    for(i=0;i<32;i++) ikey[i] = ~key[i];

        little = (char*)malloc(root * sizeof(char));
        for(i=0;i<root;i++){
            little[i] = 'p';
        }
        for(i=2;i<root2;i++){
            if(little[i] == 'p'){
                for(j=i*i;j<root;j=j+i){
                    little[j] = '*';
                }
            }
        }
        j = 0;
        for(i=3;i<root;i++){
            if(little[i] == 'p'){
                seed[j] = i;
                j++;
            }
        }
        free(little);

    if(myid == 0){
        endwtime = MPI_Wtime();
        printf("prepare time proc %d = %f\n", myid,endwtime-startwtime); 
        printf("num %d root %d root2 %d found %d prime number\n",
                 num,root,root2,j);
    }

    for(m=0;m<j;m++){
        p = seed[m];
        i = myid * n / p / 2;
        if (myid == 0) 
            i = (i+p)*p;
        else
            i = i * p * 2 + p;
        i = i - myid * n;
        if(i<0) i = i+2*p;
        if(p>30)
            for(;i<n;i=i+2*p){
                arr[i>>6] = arr[i>>6] & ikey[(i/2)&31];
            }
        else{
            mask[0] = 0;
            for(k=0;k<64;k=k+2*p){
                mask[0] = mask[0] + key[k/2];
            }
            imask[0] = ~mask[0];
            for(k=3;k<64;k=k+2){
                mask[k/2] = mask[k/2-1] << 1;
                if(k%p == 0)
                    mask[k/2] = mask[k/2] +1;
            }
            for(k=1;k<32;k++)
                imask[k] = ~mask[k];
           
            for(k=1;k<64;k=k+2){        
                t2 = (63 - k + 2*p)/p/2;
                t2 = t2 * p * 2;
                tran[k/2] = ((k+t2)/2) % 32;
            }
            if(myid == 0){
                t1 = p + 64;
                for(;i<t1;i=i+2*p)
                    arr[i>>6] = arr[i>>6] & ikey[(i/2)&31];
            }
            t1 = i / 64;
            t2 = (i/2) % 32;
            for( ; t1<n64 ; t1++){
                arr[t1] = arr[t1] & imask[t2];
                t2 = tran[t2];
            }
        }
    }

/*    MPI_Barrier(MPI_COMM_WORLD);*/

    endwtime = MPI_Wtime();
    printf("wall clock time proc %d = %f\n", myid,endwtime-startwtime); 

    if(myid == numprocs - 1){
        printf("result:\n");
        k = 0;
        for(i=n-101;i<n;i=i+2){
            if((arr[i/64] & key[(i/2)%32]) != 0){
                k++;
                printf("%d  \n",i+n*myid);
            }
        }
        printf("%x  %x  %x  %x  \n",
                    arr[n64-1],arr[n64-2],arr[n64-3],arr[n64-4]);
    }
    free(seed);
    free(arr);

    MPI_Finalize();
}
