CC=mpicc
CFLAGS=-std=c++11	

all:
		$(CC) $(CFLAGS) -o rot_A_rot_B distributed_1a_rotate_a_rotate_b.cpp
		$(CC) $(CFLAGS) -o rot_A_broadcast_B distributed_1a_rotate_a_broadcast_b.cpp
		$(CC) $(CFLAGS) -o broadcast_A_broadcast_B distributed_1a_broadcast_a_broadcast_b.cpp
		$(CC) $(CFLAGS) -o _rot_A_rot_B distributed_1d_rotate_a_rotate_b.cpp
		$(CC) $(CFLAGS) -o _rot_A_broadcast_B distributed_1d_rotate_a_broadcast_b.cpp
		$(CC) $(CFLAGS) -o _broadcast_A_broadcast_B distributed_1d_broadcast_a_broadcast_b.cpp
		$(CC) $(CFLAGS) -o 2a_scatter distributed_2a_scatter_broadcast_a_broadcast_b.cpp
		$(CC) $(CFLAGS) -o 2a_send distributed_2a_send_broadcast_a_broadcast_b.cpp

clean:
		rm -rf rot_A_rot_B
		rm -rf rot_A_broadcast_B
		rm -rf broadcast_A_broadcast_B
		rm -rf _rot_A_rot_B
		rm -rf _rot_A_broadcast_B
		rm -rf _broadcast_A_broadcast_B
		rm -rf 2a_scatter
		rm -rf 2a_send
