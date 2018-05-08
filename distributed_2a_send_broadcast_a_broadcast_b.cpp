#include<iostream>
#include<cstdio>
#include<vector>
#include<mpi.h>
#include<cmath>
#define mod 10
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



void getMul_IKJ(int **z,int **x, int **y,int z_row,int z_col, int x_row,int x_col,int y_row,int y_col, int n){
for(int i=0;i<n;++i){
		for(int k=0;k<n;++k){
			for(int j=0;j<n;++j){
				z[z_row+i][z_col+j]+=x[x_row+i][x_col+k]*y[y_row+k][y_col+j];
			}
		}
	}
}

void ParRecMM(int **z, int **x, int **y,int z_row,int z_col, int x_row,int x_col,int y_row,int y_col, int n,int m){
	if(n==m){
		return getMul_IKJ(z,x,y,z_row, z_col, x_row, x_col, y_row, y_col,n);
	}
	ParRecMM(z, x, y, z_row, z_col, x_row, x_col, y_row, y_col, n/2,m);
	ParRecMM(z, x, y, z_row, z_col+n/2, x_row, x_col, y_row, y_col+n/2, n/2,m);
   	ParRecMM(z, x, y, z_row+n/2, z_col, x_row+n/2, x_col, y_row, y_col, n/2,m);
	ParRecMM(z, x, y, z_row+n/2, z_col+n/2, x_row+n/2, x_col, y_row, y_col+n/2, n/2,m);
	// cilk_sync;
	ParRecMM(z, x, y, z_row, z_col, x_row, x_col+n/2, y_row+n/2, y_col, n/2,m);
	ParRecMM(z, x, y, z_row, z_col+n/2, x_row, x_col+n/2, y_row+n/2, y_col+n/2, n/2,m);
	ParRecMM(z , x, y, z_row+n/2, z_col, x_row+n/2, x_col+n/2, y_row+n/2, y_col, n/2,m);
	ParRecMM(z, x, y, z_row+n/2, z_col+n/2, x_row+n/2, x_col+n/2, y_row+n/2, y_col+n/2, n/2,m);
	// cilk_sync;
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

void mm_broadcast_A_broadcast_B(int **c, int **a, int **b, int myrank, int processor, int blocksize,MPI_Comm comm,MPI_Comm COL_COMM_WORLD, MPI_Comm ROW_COMM_WORLD){
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

	// for(int l=1;l<=sqrtP;++l){
	// 	performMatMul(c,a,b, myrank,blocksize);
	// 	if(l < sqrtP){
	// 		MPI_Sendrecv_replace(&(a[0][0]),blocksize*blocksize, MPI_INT,(row*sqrtP+(sqrtP+col-1)%sqrtP),tag, (row*sqrtP+(sqrtP+col+1)%sqrtP) ,tag,comm,&status);
	// 		MPI_Sendrecv_replace(&(b[0][0]),blocksize*blocksize, MPI_INT,(col+sqrtP*((sqrtP+row-1)%sqrtP)),tag, (col+sqrtP*((sqrtP+row+1)%sqrtP)) ,tag,comm,&status);
	// 	}
	// }
	int **local_allocated_buffer;
	malloc2DInt(&local_allocated_buffer, blocksize, blocksize);
	int **local_allocated_buffer_a;
	malloc2DInt(&local_allocated_buffer_a, blocksize, blocksize);
	for(int l=1;l<=sqrtP;++l){
		int k=(l-1);
		if(k==col){
			// copy2DMatrix(a,local_allocated_buffer_a,blocksize,blocksize);
			local_allocated_buffer_a =a ;
		}
		MPI_Bcast(&(local_allocated_buffer_a[0][0]), blocksize*blocksize, MPI_INT, k, ROW_COMM_WORLD);
		if(k==row){
			// local_buffer_pntr=&(b[0][0]);
			// copy2DMatrix(b,local_allocated_buffer,blocksize,blocksize);
			local_allocated_buffer=b;
		}

		MPI_Bcast(&(local_allocated_buffer[0][0]), blocksize*blocksize, MPI_INT, k, COL_COMM_WORLD);

		// performMatMul(c,local_allocated_buffer_a,local_allocated_buffer,myrank,blocksize);
		ParRecMM(c,local_allocated_buffer_a,local_allocated_buffer,0,0,0,0,0,0,blocksize,4);
		// if(l < sqrtP){
		// 	MPI_Sendrecv_replace(&(a[0][0]),blocksize*blocksize, MPI_INT,(row*sqrtP+(sqrtP+col-1)%sqrtP),tag, (row*sqrtP+(sqrtP+col+1)%sqrtP) ,tag,MPI_COMM_WORLD,&status);
		// }
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
	MPI_Comm COL_COMM_WORLD;
	MPI_Comm_split(MPI_COMM_WORLD, myrank % blockcount, myrank, &COL_COMM_WORLD);
	MPI_Comm ROW_COMM_WORLD;
	MPI_Comm_split(MPI_COMM_WORLD, myrank/blockcount, myrank, &ROW_COMM_WORLD);
	int **A;
	int **B;
	int **C;
	int tag_A=1;
	int tag_B=2;
	int tag_C=3;
	int sqrtP=sqrt(world_size);
	if(myrank==0){
		malloc2DInt(&A, n, n);
		malloc2DInt(&B, n, n);
		malloc2DInt(&C, n, n);

		int rank_local=1;
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){
				A[i][j] = (i*j+1) % 10;
				B[i][j] = (i*j+1) % 10;
				C[i][j] = 0;
			}
		}
		cout<<"A =>\n";
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){
				cout<<A[i][j]<<" ";
			}
			cout<<"\n";
		}
		cout<<"****************************\n\n";
		cout<<"B=>\n";
		print(B, n, myrank);
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
		// free2DIntArr(&A);
		// free2DIntArr(&B);
	}

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

	// print(smallMat_A, blocksize, myrank);
	// cout<<"*************HHHHHHHHH\n";
	// print(smallMat_B, blocksize, myrank);
	// cout<<"*************\n";
	// print(smallMat_C, blocksize, myrank);
	
	// performMatMul(smallMat_C,smallMat_A, smallMat_B, myrank,blocksize);
	mm_broadcast_A_broadcast_B(smallMat_C, smallMat_A, smallMat_B, myrank, world_size, blocksize,MPI_COMM_WORLD,COL_COMM_WORLD,ROW_COMM_WORLD);
	/***********SENDING MATRIX AFTER CALCULATING*************************/

	print(smallMat_C,blocksize,myrank);
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
				mergeMatrix(C, smallMat_C_received, row, col, blocksize);

				cnt++;
				
			}
		}

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
		cout<<"C =>\n";
		print(C, n, myrank);
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
