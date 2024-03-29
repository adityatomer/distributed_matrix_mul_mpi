#include<iostream>
#include<cstdio>
#include<vector>
#include<mpi.h>
#include<cmath>
#include<chrono>
#include<stdio.h>
#include<chrono>
#include<iostream>
#include<math.h>
#include<time.h>
#include<array>
#include <ctime>
#include <stdlib.h>
#define mod 10
// using namespace std;
using namespace std;
// typedef int64_t __int64;
typedef std::vector<std::vector<int> > Matrix;

Matrix getSmallerMatrix(Matrix m, int row_st,int row_end, int col_st, int col_end){
	Matrix newm;
	for(int i=0;i<row_end-row_st;++i){
		newm.push_back(std::vector<int>());
		for(int j=0;j<col_end-col_st;++j){
			newm[i].push_back(m[row_st+i][col_st+j]);
		}
	}
	return newm;
}

//0 for Head
//1 for Tail
time_t g_seed;
void seedRandomNumber(){
	srand(time(NULL));
	g_seed=rand();
}

inline uint64_t getRandomNumber() { 
	g_seed = (214013*g_seed+2531011); 
	return ((g_seed>>16)&0x7FFF)%mod; 
} 


Matrix getMatrixOfSizeR(int n, bool isRandom=true){
	Matrix A;
	for(int i=0;i<n;++i){
		A.push_back(std::vector<int>());
		for(int j=0;j<n;++j){
			if(isRandom){
				A[i].push_back(getRandomNumber());
			}else{
				A[i].push_back(j);
			}
			
		}
	}
	return A;
}

void print(Matrix A){
	std::cout<<"print called with A size: "<<A.size()<<std::endl;
	int n=A.size();
	for(int i=0;i<n;++i){
		for(int j=0;j<n;++j){
			std::cout<<A[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<"************\n";
}
void print(int rank, std::string stg, int **A, int n){
        std::cout<<"\nprint called with (array) size: "<<n<<std::endl;
        for(int i=0;i<n;++i){
                for(int j=0;j<n;++j){
                        std::cout<<stg<<"_"<<rank<<"_"<<A[i][j]<<" ";
                }
                std::cout<<std::endl;
        }
        std::cout<<"************\n";
}
void print1DVector(std::vector<int>testVector, std::string stg){
	for(int i=0;i<testVector.size();++i){
        	std::cout<<stg<<" "<<testVector[i]<<"\n";
        }
	std::cout<<"\n";
}

void allocate2dVector(){

}

int malloc2DInt(int ***array, int row, int col) { 
  
     int *pointer = (int *)malloc(row*col*sizeof(int));
     if (!pointer) return -1;
 

     (*array) = (int **)malloc(row*sizeof(int*));
     if (!(*array)) {
        free(pointer);
        return -1;
     }
 
     for (int i=0; i<row; i++)
        (*array)[i] = &(pointer[i*col]);
 
     return 0;
}

void performMatMul(int **C,int **A, int **B, int myrank, int n){
	for (int i = 0; i < n; i++)
    	{
        	for (int k = 0; k < n; k++)
        	{
            	
            		for (int j = 0; j < n; j++)
                		C[i][j] += A[i][k]*B[k][j];
        	}
    	}	
}

int free2DIntArr(int ***array) {
    free(&((*array)[0][0]));

    free(*array);

    return 0;
}

void mm_rotate_A_rotate_B(int **c, int **a, int **b, int n, int myrank, int processor, int blocksize, MPI_Datatype smallMatType,MPI_Comm comm){
	MPI_Status status;
	int tag=123456;
	int sqrtP=sqrt(processor);

	int row=myrank/sqrtP;
	int col=myrank%sqrtP;

	int dest_A_row=row;
	int dest_A_col=(sqrtP+col-row)%sqrtP;
	int dest_A_rank=sqrtP*dest_A_row+dest_A_col;

	int src_A_row=row;
	int src_A_col=(-sqrtP+col+row)%sqrtP;
	int src_A_rank=((myrank+row) % sqrtP) + (row*sqrtP);

	int dest_B_row=(sqrtP+col-row)%sqrtP;
	int dest_B_col=col;
	int dest_B_rank=sqrtP*((sqrtP+row-col)%sqrtP)+col;
	int src_B_rank= sqrtP*((row+col)%sqrtP)+col;

	MPI_Sendrecv_replace(&(a[0][0]),blocksize*blocksize, MPI_INT,dest_A_rank,tag,src_A_rank,tag,comm,&status);
	MPI_Sendrecv_replace(&(b[0][0]),blocksize*blocksize, MPI_INT,dest_B_rank,tag,src_B_rank,tag,comm,&status);

	for(int l=1;l<=sqrtP;++l){
		performMatMul(c,a,b, myrank,blocksize);
		if(l < sqrtP){
			MPI_Sendrecv_replace(&(a[0][0]),blocksize*blocksize, MPI_INT,(row*sqrtP+(sqrtP+col-1)%sqrtP),tag, (row*sqrtP+(sqrtP+col+1)%sqrtP) ,tag,comm,&status);
			MPI_Sendrecv_replace(&(b[0][0]),blocksize*blocksize, MPI_INT,(col+sqrtP*((sqrtP+row-1)%sqrtP)),tag, (col+sqrtP*((sqrtP+row+1)%sqrtP)) ,tag,comm,&status);
		}
	}
}

int main(int argc, char *argv[]){
	seedRandomNumber();
	int n=atoi(argv[1]);
	int world_size,myrank;
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);	
	int blockcount = sqrt(world_size);
	int blocksize  = n/blockcount;
	int sendcount[world_size];
	int displaycount[world_size];
	int **A;
	int **B;
	int **C;
	int **smallMat_A;
	int **smallMat_B;
	int **smallMat_C;
	malloc2DInt(&smallMat_A, blocksize, blocksize);
	malloc2DInt(&smallMat_B, blocksize, blocksize);
	malloc2DInt(&smallMat_C, blocksize, blocksize);
	__int64 now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	/*__cilkrts_set_param("nworkers", argv[1]);*/
	if(myrank==0){
		malloc2DInt(&A, n, n);
		malloc2DInt(&B, n, n);
		malloc2DInt(&C, n, n);
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){
				A[i][j] = (i*j+1) % 10;
				B[i][j] = (i*j+1) % 10;
				C[i][j] = 0;
			}
		}
		std::cout<<"A =>\n";
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){
				std::cout<<A[i][j]<<" ";
			}
			std::cout<<"\n";
		}
		/*
		std::cout<<"B =>\n";
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){
				std::cout<<B[i][j]<<" ";
			}
			std::cout<<"\n";
		}
		*/
		std::cout<<"world_size: "<<world_size<<" "<<" blockcount: "<<blockcount<<" blocksize: "<<blocksize<<"\n";
	}
	
	//smallMat.resize(blocksize, std::vector<int>(blocksize,0));
	MPI_Datatype type, smallMatType;
	int sizes[2]={n,n};
	int subSizes[2]={blocksize,blocksize};
	int starts[2]={0,0};

	MPI_Type_create_subarray(2, sizes, subSizes, starts, MPI_ORDER_C, MPI_INT, &type);
	MPI_Type_create_resized(type, 0, blocksize*sizeof(int), &smallMatType);
	MPI_Type_commit(&smallMatType);
	int *globalptr_A=NULL;
	int *globalptr_B=NULL;
	int *globalptr_C=NULL;
    	if(myrank==0){
		globalptr_A=&(A[0][0]);
		globalptr_B=&(B[0][0]);
		globalptr_C=&(C[0][0]);
		int cnt=0;
		for(int i=0; i<blockcount; ++i){
			for(int j=0; j<blockcount; ++j){
				sendcount[cnt]=1;
				displaycount[cnt++]= i * blockcount * blocksize + j;
			}
		}
	}
	MPI_Scatterv(globalptr_A, sendcount, displaycount, smallMatType, &(smallMat_A[0][0]),
                 blocksize*blocksize, MPI_INT,0, MPI_COMM_WORLD);
	MPI_Scatterv(globalptr_B, sendcount, displaycount, smallMatType, &(smallMat_B[0][0]),
                 blocksize*blocksize, MPI_INT,0, MPI_COMM_WORLD);
	MPI_Scatterv(globalptr_C, sendcount, displaycount, smallMatType, &(smallMat_C[0][0]),
                 blocksize*blocksize, MPI_INT,0, MPI_COMM_WORLD);
	for (int p=0; p<world_size; p++) {
		MPI_Barrier(MPI_COMM_WORLD);
	}
	// performMatMul(smallMat_C,smallMat_A, smallMat_B, myrank,blocksize);
	mm_rotate_A_rotate_B(smallMat_C, smallMat_A, smallMat_B, blocksize, myrank, world_size, blocksize, smallMatType,MPI_COMM_WORLD);

	MPI_Gatherv(&(smallMat_C[0][0]), blocksize*blocksize,  MPI_INT, globalptr_C, sendcount, 
		displaycount, smallMatType,0, MPI_COMM_WORLD);

	__int64 now1 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	float timeTaken=float (now1- now);
	

	free2DIntArr(&smallMat_A);
	free2DIntArr(&smallMat_B);
	free2DIntArr(&smallMat_C);
	if(myrank==0){
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){	
				std::cout<<C[i][j]<<" ";
			}
			std::cout<<"\n";
		}
		free2DIntArr(&A);
		free2DIntArr(&B);
		free2DIntArr(&C);
	}		
	MPI_Type_free(&smallMatType);
	MPI_Finalize();	
	std::cout<<"time taken: "<<timeTaken<<"\n";

	return 0;
}
