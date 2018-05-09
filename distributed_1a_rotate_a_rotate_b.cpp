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

void print(int **mat, int n, int rank){
	cout<<"printprint rank: "<<rank<<" "<<"\n";
	for(int i=0;i<n;++i){
		for(int j=0;j<n;++j){
			cout<<mat[i][j]<<" ";
		}cout<<"\n";
	}
}

void mm_rotate_A_rotate_B(int **c, int **a, int **b, int myrank, int processor, int blocksize,MPI_Comm comm){
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
	int **A;
	int **B;
	int **C;
	int tag_A=1;
	int tag_B=2;
	int tag_C=3;
	int sqrtP=sqrt(world_size);
	__int64 now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
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
				cout<<"row start: "<<row<<" "<<"row end: "<<row+blocksize<<" col: "<<col<<" col end: "<<col+blocksize<<"\n";
				smallMat_A=getSmallerMatrix( A,  blocksize, row, row+blocksize,  col,  col+blocksize);
				// print(smallMat_A,n, myrank);
				smallMat_B=getSmallerMatrix( B,  blocksize, row, row+blocksize,  col,  col+blocksize);
				smallMat_C=getSmallerMatrix( C,  blocksize, row, row+blocksize,  col,  col+blocksize);

				MPI_Isend(&(smallMat_A[0][0]), blocksize*blocksize, MPI_INT, rankTemp, tag_A, MPI_COMM_WORLD, &request_A[i*sqrtP+j]);
				MPI_Isend(&(smallMat_B[0][0]), blocksize*blocksize, MPI_INT, rankTemp, tag_B, MPI_COMM_WORLD, &request_B[i*sqrtP+j]);
				MPI_Isend(&(smallMat_C[0][0]), blocksize*blocksize, MPI_INT, rankTemp, tag_C, MPI_COMM_WORLD, &request_C[i*sqrtP+j]);
				rankTemp++;
			}
		}
		// return 0;

		// MPI_Wait(&request_A[0],&status);
		// MPI_Wait(&request_B[0],&status);
		// MPI_Wait(&request_C[0],&status);
		free2DIntArr(&A);
		free2DIntArr(&B);
		// free2DIntArr(&C);
	}
	// return 0;

	MPI_Request request_recv_matrices[3];
	
	int **smallMat_A;
	int **smallMat_B;
	int **smallMat_C;
	malloc2DInt(&smallMat_A, blocksize, blocksize);
	malloc2DInt(&smallMat_B, blocksize, blocksize);
	malloc2DInt(&smallMat_C, blocksize, blocksize);

	MPI_Irecv(&(smallMat_A[0][0]), blocksize*blocksize, MPI_INT, 0, tag_A, MPI_COMM_WORLD,&request_recv_matrices[0]);
	MPI_Irecv(&(smallMat_B[0][0]), blocksize*blocksize, MPI_INT, 0, tag_B, MPI_COMM_WORLD,&request_recv_matrices[1]);
	MPI_Irecv(&(smallMat_C[0][0]), blocksize*blocksize, MPI_INT, 0, tag_C, MPI_COMM_WORLD,&request_recv_matrices[2]);
	
	MPI_Wait(&request_recv_matrices[0], &status);
	MPI_Wait(&request_recv_matrices[1], &status);
	MPI_Wait(&request_recv_matrices[2], &status);

	// if(myrank == 0){
	// cout<<"the code owdwdvcoenfovin";
	// print(smallMat_A, blocksize, myrank);
	// }
	// cout<<"*************HHHHHHHHH\n";
	// print(smallMat_B, blocksize, myrank);
	// cout<<"*************\n";
	// print(smallMat_C, blocksize, myrank);
	// return 0;
	// performMatMul(smallMat_C,smallMat_A, smallMat_B, myrank,blocksize);
	mm_rotate_A_rotate_B(smallMat_C, smallMat_A, smallMat_B, myrank, world_size, blocksize,MPI_COMM_WORLD);
	// if(myrank == 0){
	// print(smallMat_C, blocksize, myrank);
	// }
	/***********SENDING MATRIX AFTER CALCULATING*************************/
	int row=myrank/sqrtP;
	int col=myrank%sqrtP;
	MPI_Request request_send;

	MPI_Isend(&(smallMat_C[0][0]), blocksize*blocksize, MPI_INT, 0, tag_C, MPI_COMM_WORLD, &request_send);
	if(myrank==0){
		

		// MPI_Request request_send_root[world_size];

		// MPI_Wait(&request_send,&status);
		int **smallMat_C_received;
		MPI_Request request_recv[world_size];
		int cnt=0;
		for(int i=0;i<sqrtP;++i){
			for(int j=0;j<sqrtP;++j){
				int row=i*blocksize;
				int col=j*blocksize;
				malloc2DInt(&smallMat_C_received, blocksize, blocksize);
				MPI_Irecv(&(smallMat_C_received[0][0]), blocksize*blocksize, MPI_INT, cnt, tag_C, MPI_COMM_WORLD,&request_recv[cnt]);
				MPI_Wait(&request_recv[cnt], &status);
				
				// print(smallMat_C,blocksize,myrank);
				mergeMatrix(C, smallMat_C_received, row, col, blocksize);

				cnt++;
				
			}
		}
		__int64 now1 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		float timeTaken=float (now1- now);	
		// cnt=0;
		// for(int i=0;i<sqrtP;++i){
		// 	for(int j=0;j<sqrtP;++j){
		// 		MPI_Wait(&request_recv[cnt], &status);
		// 		print(smallMat_C_received, blocksize, myrank);
		// 		mergeMatrix(C, smallMat_C_received, row, col, blocksize);
		// 		cnt++;
		// 	}
		// }
		cout<<"****************************\n\n";
		cout<<"Matrix multiplication done:: C =>\n";
		if(flag==1){
			print(C, n, myrank);
		}
		std::cout<<"time taken: "<<timeTaken<<"\n";
		// cout<<"Printing Matrix C\n";
		// for(int i=0;i<n;++i){
		// 	for(int j=0;j<n;++j){
		// 		cout<<C[i][j]<<" ";
		// 	}cout<<"\n";
		// }
		free2DIntArr(&smallMat_C_received);
	}

	free2DIntArr(&smallMat_A);
	free2DIntArr(&smallMat_B);
	free2DIntArr(&smallMat_C);
	// free2DIntArr(&A);
	// free2DIntArr(&B);
	// free2DIntArr(&C);
	// free2DIntArr(&smallMat_B);
	// free2DIntArr(&smallMat_C);
	// if(myrank==0){
	// 	for (int i=0; i<n; i++) {
	// 		for (int j=0; j<n; j++){	
	// 			cout<<C[i][j]<<" ";
	// 		}
	// 		cout<<"\n";
	// 	}
	// 	free2DIntArr(&A);
	// 	free2DIntArr(&B);
	// 	free2DIntArr(&C);
	// }		
	
	MPI_Finalize();	

	return 0;
}
