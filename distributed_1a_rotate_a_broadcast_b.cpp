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
#define mod 100
// typedef int64_t __int64;
using namespace std;
typedef std::vector<std::vector<int> > Matrix;

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
int ** getSmallerMatrix(int **m, int n,int row_st,int row_end, int col_st, int col_end){
	int **newm; 
	malloc2DInt(&newm, n, n);
	for(int i=0;i<row_end-row_st;++i){
		for(int j=0;j<col_end-col_st;++j){
			newm[i][j] = (m[row_st+i][col_st+j]);
		}
	}
	return newm;
}
void mergeMatrix(int **mat, int **smallMat, int small_n_row, int small_n_col, int small_n_size){
	// **mergeMatr(int **mat, int **smallMat, int n, int small_n_row, int small_n_col, int small_n_size)
	for(int i=0; i<small_n_size; ++i){
		for(int j=0; j<small_n_size; ++j){
			mat[small_n_row+i][small_n_col+j] = smallMat[i][j];
		}
	}
}
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
	cout<<"print called with A size: "<<A.size()<<endl;
	int n=A.size();
	for(int i=0;i<n;++i){
		for(int j=0;j<n;++j){
			cout<<A[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<"************\n";
}
void print(int rank, std::string stg, int **A, int n){
        cout<<"\nprint called with (array) size: "<<n<<endl;
        for(int i=0;i<n;++i){
                for(int j=0;j<n;++j){
                        cout<<stg<<"_"<<rank<<"_"<<A[i][j]<<" ";
                }
                cout<<endl;
        }
        cout<<"************\n";
}
void print1DVector(std::vector<int>testVector, std::string stg){
	for(int i=0;i<testVector.size();++i){
        	cout<<stg<<" "<<testVector[i]<<"\n";
        }
	cout<<"\n";
}

void allocate2dVector(){

}



void performMatMul(int **C,int **A, int **B, int n){
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

void print(int **mat, int n, int rank){
	cout<<"printprint rank: "<<rank<<" "<<"\n";
	for(int i=0;i<n;++i){
		for(int j=0;j<n;++j){
			cout<<mat[i][j]<<" ";
		}cout<<"\n";
	}
}

void mm_rotate_A_broadcast_B(int **c, int **a, int **b, int myrank, int world_size, int blocksize, MPI_Comm COL_COMM_WORLD	){
	MPI_Status status;
	int tag=123456;
	int sqrtP=sqrt(world_size);

	int row=myrank/sqrtP;
	int col=myrank%sqrtP;
	int *local_buffer_pntr;
	int **local_allocated_buffer;
	malloc2DInt(&local_allocated_buffer, blocksize, blocksize);
	
	for(int l=1;l<=sqrtP;++l){
		int k=(col+l-1)%sqrtP;
		if(k==row){
			local_allocated_buffer=b;
		}

		MPI_Bcast(&(local_allocated_buffer[0][0]), blocksize*blocksize, MPI_INT, k, COL_COMM_WORLD);
		performMatMul(c,a,local_allocated_buffer,blocksize);
		if(l < sqrtP){
			MPI_Sendrecv_replace(&(a[0][0]),blocksize*blocksize, MPI_INT,(row*sqrtP+(sqrtP+col-1)%sqrtP),tag, (row*sqrtP+(sqrtP+col+1)%sqrtP) ,tag,MPI_COMM_WORLD,&status);
		}
	}
}

int main(int argc, char *argv[]){
	seedRandomNumber();
	int n=pow(2,atoi(argv[1]));
	int flag = atoi(argv[2]);
	int world_size,myrank;
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);	
	int blockcount = sqrt(world_size);
	int blocksize  = n/blockcount;
	int sendcount[world_size];
	int displaycount[world_size];
	MPI_Comm COL_COMM_WORLD;
	MPI_Comm_split(MPI_COMM_WORLD, myrank % blockcount, myrank, &COL_COMM_WORLD);
	int **A;
	int **B;
	int **C;
	int tag_A=1;
	int tag_B=2;
	int tag_C=3;
	int sqrtP=sqrt(world_size);
	__int64 now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	int **small_A;
	int **small_B;
	int **small_C;
	if(myrank==0){
		malloc2DInt(&A, n, n);
		malloc2DInt(&B, n, n);
		malloc2DInt(&C, n, n);

		int rank_local=1;
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){
				A[i][j] = getRandomNumber();
				if(getRandomNumber()%2 ==0){
				A[i][j] *= -1;
				}

				B[i][j] =getRandomNumber();
				if(getRandomNumber()%2 ==0){
				B[i][j] *= -1;
				}
				C[i][j] = 0;
			}
		}
		cout<<"A completed =>\n";
		if(flag==1){
			for (int i=0; i<n; i++) {
				for (int j=0; j<n; j++){
					cout<<A[i][j]<<" ";
				}
				cout<<"\n";
			}
		}
		cout<<"****************************\n\n";

		cout<<"B completed =>\n";
		if(flag==1){
			print(B, n, myrank);
		}
		int rankTemp=0;
		MPI_Request request_A[world_size];
		MPI_Request request_B[world_size];
		MPI_Request request_C[world_size];
		
		for(int i=0;i<sqrtP;++i){
			for(int j=0;j<sqrtP;++j){
				int row=i*blocksize;
				int col=j*blocksize;
				int **smallMat_A;
				int **smallMat_B;
				int **smallMat_C;
				// cout<<"row start: "<<row<<" "<<"row end: "<<row+blocksize<<" col: "<<col<<" col end: "<<col+blocksize<<"\n";
				smallMat_A=getSmallerMatrix( A,  blocksize, row, row+blocksize,  col,  col+blocksize);
				smallMat_B=getSmallerMatrix( B,  blocksize, row, row+blocksize,  col,  col+blocksize);
				smallMat_C=getSmallerMatrix( C,  blocksize, row, row+blocksize,  col,  col+blocksize);

				if(rankTemp !=0){	
					MPI_Send(&(smallMat_A[0][0]), blocksize*blocksize, MPI_INT, rankTemp, tag_A, MPI_COMM_WORLD);
					MPI_Send(&(smallMat_B[0][0]), blocksize*blocksize, MPI_INT, rankTemp, tag_B, MPI_COMM_WORLD);
					MPI_Send(&(smallMat_C[0][0]), blocksize*blocksize, MPI_INT, rankTemp, tag_C, MPI_COMM_WORLD);
				}
				rankTemp++;
			}
		}
		small_A=getSmallerMatrix( A,  blocksize, 0, 0+blocksize,  0,  0+blocksize);
		small_B=getSmallerMatrix( B,  blocksize, 0, 0+blocksize,  0,  0+blocksize);
		small_C=getSmallerMatrix( C,  blocksize, 0, 0+blocksize,  0,  0+blocksize);
		free2DIntArr(&A);
		free2DIntArr(&B);
	}

	MPI_Request request_recv_matrices[3];
	int **smallMat_A;
	int **smallMat_B;
	int **smallMat_C;
	malloc2DInt(&smallMat_A, blocksize, blocksize);
	malloc2DInt(&smallMat_B, blocksize, blocksize);
	malloc2DInt(&smallMat_C, blocksize, blocksize);
	if(myrank==0){
	smallMat_C =small_C;
	smallMat_A = small_A;
	smallMat_B= small_B;
	}
	if(myrank!=0){
		
		MPI_Status status1[3];
		MPI_Recv(&(smallMat_A[0][0]), blocksize*blocksize, MPI_INT, 0, tag_A, MPI_COMM_WORLD,&status1[0]);
		MPI_Recv(&(smallMat_B[0][0]), blocksize*blocksize, MPI_INT, 0, tag_B, MPI_COMM_WORLD,&status1[1]);
		MPI_Recv(&(smallMat_C[0][0]), blocksize*blocksize, MPI_INT, 0, tag_C, MPI_COMM_WORLD,&status1[2]);
		

		// print(smallMat_A, blocksize, myrank);
		// cout<<"*************HHHHHHHHH\n";
		// print(smallMat_B, blocksize, myrank);
		// cout<<"*************\n";
		// print(smallMat_C, blocksize, myrank);
	}
	
	mm_rotate_A_broadcast_B(smallMat_C, smallMat_A, smallMat_B, myrank, world_size, blocksize,COL_COMM_WORLD);
	/***********SENDING MATRIX AFTER CALCULATING*************************/
	// print(smallMat_C, blocksize, myrank);
	int row=myrank/sqrtP;
	int col=myrank%sqrtP;
	MPI_Request request_send;

	if(myrank !=0){
	MPI_Send(&(smallMat_C[0][0]), blocksize*blocksize, MPI_INT, 0, tag_C, MPI_COMM_WORLD);
	}	

	if(myrank==0){
		mergeMatrix(C, smallMat_C, 0, 0, blocksize);
		int **smallMat_C_received;
		MPI_Request request_recv[world_size];
		MPI_Status status[world_size];
		int cnt=1;
		for(int i=0;i<sqrtP;++i){
			for(int j=0;j<sqrtP;++j){
				if(i==0 && j==0){
					continue;
				}
				int row=i*blocksize;
				int col=j*blocksize;
				malloc2DInt(&smallMat_C_received, blocksize, blocksize);
				MPI_Recv(&(smallMat_C_received[0][0]), blocksize*blocksize, MPI_INT, cnt, tag_C, MPI_COMM_WORLD,&status[cnt]);
				// MPI_Wait(&request_recv[cnt], &status);
				mergeMatrix(C, smallMat_C_received, row, col, blocksize);

				cnt++;
				
			}
		}
		__int64 now1 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		float timeTaken=float (now1- now);
		cout<<"****************************\n\n";
		cout<<"Matrix multiplication done:: C =>\n";
		if(flag==1){
			print(C, n, myrank);
		}
		std::cout<<"time taken: "<<timeTaken<<"\n";
		if(sqrtP!=1){
		free2DIntArr(&smallMat_C_received);
		}
	}

	free2DIntArr(&smallMat_A);
	free2DIntArr(&smallMat_B);
	free2DIntArr(&smallMat_C);
	
	MPI_Finalize();	

	return 0;
}
