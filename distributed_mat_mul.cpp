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


void mm_rotate_A_rotate_B(Matrix c, Matrix A, Matrix B, int n, int p){

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
        cout<<"print called with A(array) size: "<<n<<endl;
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
int malloc2dint(int ***array, int n, int m) {

    /* allocate the n*m contiguous items */
    int *p = (int *)malloc(n*m*sizeof(int));
    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*array) = (int **)malloc(n*sizeof(int*));
    if (!(*array)) {
       free(p);
       return -1;
    }

    /* set up the pointers into the contiguous memory */
    for (int i=0; i<n; i++)
       (*array)[i] = &(p[i*m]);

    return 0;
}

void performMatMul(int **C,int **A, int **B, int myrank, int n){
	print(myrank, "mat_A",A,n);
	print(myrank, "mat_B",B,n);
	

	for(int i=0;i<n;++i){
		for(int k=0;k<n;++k){
			for(int j=0;j<n;++j){
				C[i][j]+=A[i][k]*B[k][j];
			}
		}
	}
	print(myrank, "mat_C",C, n);
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
	if(myrank==0){
		malloc2dint(&A, n, n);
		malloc2dint(&B, n, n);
		malloc2dint(&C, n, n);
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){
				A[i][j] = i*j+1;
				B[i][j] = i*j+1;
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
		cout<<"B =>\n";
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++){
				cout<<B[i][j]<<" ";
			}
			cout<<"\n";
		}
		cout<<"world_size: "<<world_size<<" "<<" blockcount: "<<blockcount<<" blocksize: "<<blocksize<<endl;
	}

	int **smallMat_A;
	int **smallMat_B;
	int **smallMat_C;
	malloc2dint(&smallMat_A, blocksize, blocksize);
	malloc2dint(&smallMat_B, blocksize, blocksize);
	malloc2dint(&smallMat_C, blocksize, blocksize);
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
		for(int i=0; i<cnt; ++i){
			cout<<"displaycount: "<<displaycount[i]<<"\n";
		}
	}


	MPI_Scatterv(globalptr_A, sendcount, displaycount, smallMatType, &(smallMat_A[0][0]),
                 world_size, MPI_INT,0, MPI_COMM_WORLD);
	MPI_Scatterv(globalptr_B, sendcount, displaycount, smallMatType, &(smallMat_B[0][0]),
                 world_size, MPI_INT,0, MPI_COMM_WORLD);
	MPI_Scatterv(globalptr_C, sendcount, displaycount, smallMatType, &(smallMat_C[0][0]),
                 world_size, MPI_INT,0, MPI_COMM_WORLD);
	

	for (int p=0; p<world_size; p++) {
		MPI_Barrier(MPI_COMM_WORLD);
	}

	performMatMul(smallMat_C,smallMat_A, smallMat_B, myrank,blocksize);
	


	MPI_Type_free(&smallMatType);
	cout<<"-----------------\n";
	// print(myrank, "mat_A",smallMat_A,blocksize);
	// print(myrank, "mat_B",smallMat_B,blocksize);
	// print(myrank, "mat_C",smallMat_C,blocksize);

	

	return 0;
}
