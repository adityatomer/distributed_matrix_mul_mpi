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

int ** getSmallerMatrix(int **m, int row_st,int row_end, int col_st, int col_end){
	int **newm; 
	malloc2DInt(&newm, n, n);
	for(int i=0;i<row_end-row_st;++i){
		for(int j=0;j<col_end-col_st;++j){
			newm[i][j] = (m[row_st+i][col_st+j]);
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
	int tag=123;
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
		int sqrtP=sqrt(world_size);
		int rankTemp=0;
		for(int i=0;i<sqrtP;++i){
			for(int j=0;j<sqrtP;++j){
				rankTemp++;
				int row=rankTemp/sqrtP;
				int col=rankTemp%sqrtP;
				int **smallMat_A;
				int **smallMat_B;
				int **smallMat_C;

				smallMat_A=getSmallerMatrix( A,  row, row+blocksize,  col,  col+blocksize);
				smallMat_B=getSmallerMatrix( B,  row, row+blocksize,  col,  col+blocksize);
				smallMat_C=getSmallerMatrix( C,  row, row+blocksize,  col,  col+blocksize);

				/*
				int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request)
				*/
				MPI_Request request;
				MPI_Isend(&(smallMat_A[0][0]), blocksize*blocksize, MPI_INT, rankTemp, tag,MPI_COMM_WORLD,&request);
				MPI_Isend(&(smallMat_B[0][0]), blocksize*blocksize, MPI_INT, rankTemp, tag,MPI_COMM_WORLD,&request);
				MPI_Isend(&(smallMat_C[0][0]), blocksize*blocksize, MPI_INT, rankTemp, tag,MPI_COMM_WORLD,&request);
			}
		}
	}

	MPI_Request request_1;
	MPI_Status status;
	int **smallMat_A;
	int **smallMat_B;
	int **smallMat_C;
	/*
	int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
	              int tag, MPI_Comm comm, MPI_Request *request)
	*/
	
	MPI_Irecv(&(smallMat_A[0][0]), blocksize*blocksize, MPI_INT, 0, tag,MPI_COMM_WORLD,&request);
	MPI_Irecv(&(smallMat_B[0][0]), blocksize*blocksize, MPI_INT, 0, tag,MPI_COMM_WORLD,&request);
	MPI_Irecv(&(smallMat_C[0][0]), blocksize*blocksize, MPI_INT, 0, tag,MPI_COMM_WORLD,&request);
	MPI_Wait(&request_1, &status);


	// performMatMul(smallMat_C,smallMat_A, smallMat_B, myrank,blocksize);
	mm_rotate_A_rotate_B(smallMat_C, smallMat_A, smallMat_B, blocksize, myrank, world_size, blocksize, smallMatType,MPI_COMM_WORLD);

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
