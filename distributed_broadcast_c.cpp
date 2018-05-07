#include<iostream>
#include<cstdio>
#include<vector>
#include<mpi.h>
#include<cmath>
#define mod 10
using namespace std;
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

void copy2DMatrix(int **mat, int **local, int row, int col){
	for(int i=0;i<row;++i){
		for(int j=0;j<col;++j){
			local[i][j]=mat[i][j];
		}
	}
}

void mm_rotate_A_broadcast_B(int **c, int **a, int **b, int myrank, int world_size, int blocksize, MPI_Comm COL_COMM_WORLD, MPI_Comm ROW_COMM_WORLD){

	MPI_Status status;
	int tag=123456;
	int sqrtP=sqrt(world_size);

	int row=myrank/sqrtP;
	int col=myrank%sqrtP;
	int **local_allocated_buffer;
	int **local_allocated_buffer_a;
	malloc2DInt(&local_allocated_buffer, blocksize, blocksize);
	malloc2DInt(&local_allocated_buffer_a, blocksize, blocksize);
	for(int l=1;l<=sqrtP;++l){
		int k=(l-1);
		if(k==col){
			copy2DMatrix(a,local_allocated_buffer_a,blocksize,blocksize);
		}
		if(k==row){
			// local_buffer_pntr=&(b[0][0]);
			copy2DMatrix(b,local_allocated_buffer,blocksize,blocksize);
		}
		MPI_Bcast(&(local_allocated_buffer_a[0][0]), blocksize*blocksize, MPI_INT, k, ROW_COMM_WORLD);

		MPI_Bcast(&(local_allocated_buffer[0][0]), blocksize*blocksize, MPI_INT, k, COL_COMM_WORLD);

		performMatMul(c,local_allocated_buffer_a,local_allocated_buffer,blocksize);
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
	int sqrtP=sqrt(world_size);
	MPI_Comm COL_COMM_WORLD;
	MPI_Comm ROW_COMM_WORLD;
	MPI_Comm_split(MPI_COMM_WORLD, myrank % sqrtP, myrank, &COL_COMM_WORLD);
	MPI_Comm_split(MPI_COMM_WORLD, sqrtP + myrank/sqrtP, myrank, &ROW_COMM_WORLD);
	int blockcount = sqrt(world_size);
	int blocksize  = n/blockcount;
	int sendcount[world_size];
	int displaycount[world_size];
	int **A;
	int **B;
	int **C;
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
		cout<<"A =>\n";
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){
				cout<<A[i][j]<<" ";
			}
			cout<<"\n";
		}
		/*
		cout<<"B =>\n";
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){
				cout<<B[i][j]<<" ";
			}
			cout<<"\n";
		}
		*/
		cout<<"world_size: "<<world_size<<" "<<" blockcount: "<<blockcount<<" blocksize: "<<blocksize<<endl;
	}
	int **smallMat_A;
	int **smallMat_B;
	int **smallMat_C;
	malloc2DInt(&smallMat_A, blocksize, blocksize);
	malloc2DInt(&smallMat_B, blocksize, blocksize);
	malloc2DInt(&smallMat_C, blocksize, blocksize);
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
	mm_rotate_A_broadcast_B(smallMat_C,smallMat_A, smallMat_B, myrank, world_size, blocksize,COL_COMM_WORLD,ROW_COMM_WORLD);

	MPI_Gatherv(&(smallMat_C[0][0]), blocksize*blocksize,  MPI_INT, globalptr_C, sendcount, 
		displaycount, smallMatType,0, MPI_COMM_WORLD);
	free2DIntArr(&smallMat_A);
	free2DIntArr(&smallMat_B);
	free2DIntArr(&smallMat_C);
	if(myrank==0){
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){	
				cout<<C[i][j]<<" ";
			}
			cout<<"\n";
		}
		free2DIntArr(&A);
		free2DIntArr(&B);
		free2DIntArr(&C);
	}		
	MPI_Type_free(&smallMatType);
	MPI_Finalize();	

	return 0;
}

