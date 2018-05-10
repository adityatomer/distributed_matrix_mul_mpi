

Go to Path =====>. "/work/05567/varun31/stampede2/Group19_VarunAgarwal_111491232#_AdityaTomer_111491409#/code/distributed_matrix_mul_mpi"

 Run make 

For question 1a : mpirun -n 4 ./rot_A_rot_B 10 0 -ppn 1
				  
				  mpirun -n 4 ./rot_A_broadcast_B 10 0 -ppn 1

				  mpirun -n 4 ./broadcast_A_broadcast_B 10 0 -ppn 1
				  

For question 1b: for((i=10;i<=14;i=i+1)); do mpirun -n 16 ./rot_A_rot_B $i 0 -ppn 1 ; done
					
				 for((i=10;i<=14;i=i+1)); do mpirun -n 4 ./rot_A_rot_B $i 0 -ppn 1 ; done
				 
				 for((i=10;i<=14;i=i+1)); do mpirun -n 1 ./rot_A_rot_B $i 0 -ppn 1 ; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 16 ./rot_A_broadcast_B $i 0 -ppn 1 ; done	

				 for((i=10;i<=14;i=i+1)); do mpirun -n 4 ./rot_A_broadcast_B $i 0 -ppn 1 ; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 1 ./rot_A_broadcast_B $i 0 -ppn 1 ; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 16 ./broadcast_A_broadcast_B $i 0 -ppn 1 ; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 4 ./broadcast_A_broadcast_B $i 0 -ppn 1 ; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 1 ./broadcast_A_broadcast_B $i 0 -ppn 1 ; done

				 

For question 1c: for((i=10;i<=14;i=i+1)); do mpirun -n 16 ./rot_A_rot_B $i 0 -ppn 68 ; done
					
				 for((i=10;i<=14;i=i+1)); do mpirun -n 4 ./rot_A_rot_B $i 0 -ppn 68 ; done
				 
				 for((i=10;i<=14;i=i+1)); do mpirun -n 1 ./rot_A_rot_B $i 0 -ppn 68 ; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 16 ./rot_A_broadcast_B $i 0 -ppn 68 ; done	

				 for((i=10;i<=14;i=i+1)); do mpirun -n 4 ./rot_A_broadcast_B $i 0 -ppn 68 ; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 1 ./rot_A_broadcast_B $i 0 -ppn 68 ; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 16 ./broadcast_A_broadcast_B $i 0 -ppn 68 ; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 4 ./broadcast_A_broadcast_B $i 0 -ppn 68; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 1 ./broadcast_A_broadcast_B $i 0 -ppn 68; done


For question 1d: mpirun -n 4 ./_broadcast_A_broadcast_B 10 0 -ppn 1


For question 1e: for((i=10;i<=14;i=i+1)); do mpirun -n 16 ./_broadcast_A_broadcast_B $i 0 -ppn 1 ; done
					
				 for((i=10;i<=14;i=i+1)); do mpirun -n 4 ./_broadcast_A_broadcast_B $i 0 -ppn 1 ; done
				 
				 for((i=10;i<=14;i=i+1)); do mpirun -n 1 ./_broadcast_A_broadcast_B $i 0 -ppn 1 ; done


				 

For question 1f: for((i=10;i<=14;i=i+1)); do mpirun -n 16 ./_broadcast_A_broadcast_B $i 0 -ppn 68 ; done
					
				 for((i=10;i<=14;i=i+1)); do mpirun -n 4 ./_broadcast_A_broadcast_B $i 0 -ppn 68 ; done
				 
				 for((i=10;i<=14;i=i+1)); do mpirun -n 1 ./_broadcast_A_broadcast_B $i 0 -ppn 68 ; done

				 

For question 2a : mpirun -n 4 ./2a_scatter 10 0 -ppn 1

				  mpirun -n 4 ./2a_send 10 0 -ppn 1

For question 2b: for((i=10;i<=14;i=i+1)); do mpirun -n 16 ./2a_scatter $i 0 -ppn 1 ; done
					
				 for((i=10;i<=14;i=i+1)); do mpirun -n 4 ./2a_scatter $i 0 -ppn 1 ; done
				 
				 for((i=10;i<=14;i=i+1)); do mpirun -n 1 ./2a_scatter $i 0 -ppn 1 ; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 16 ./2a_send $i 0 -ppn 1 ; done	

				 for((i=10;i<=14;i=i+1)); do mpirun -n 4 ./2a_send $i 0 -ppn 1 ; done

				 for((i=10;i<=14;i=i+1)); do mpirun -n 1 ./2a_send $i 0 -ppn 1 ; done

				