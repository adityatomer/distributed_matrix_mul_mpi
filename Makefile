CC=mpicc
CFLAGS=-std=c++11

all:
		$(CC) $(CFLAGS) -o rot_A_rot_B distributed_1a_rotate_a_rotate_b.cpp
		$(CC) $(CFLAGS) -o rot_A_broadcast_B distributed_1a_rotate_a_broadcast_b.cpp
		$(CC) $(CFLAGS) -o broadcast_A_broadcast_B distributed_1a_broadcast_a_broadcast_b.cpp

clean:
		rm -rf rot_A_rot_B
		rm -rf rot_A_broadcast_B
		rm -rf broadcast_A_broadcast_B