include /workspace/example.mk
#include ~/common.mk

CC=mpic++

LDIR =

OBJ1 = non_homogeneous_activity.o
OBJ2 = non_homogeneous_activity_control_rd.o
OBJ3 = non_homogeneous_activity_hanangata.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)


non_homogeneous_activity: $(OBJ1)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

non_homogeneous_activity_control_rd: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

non_homogeneous_activity_hanangata: $(OBJ3)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)


all: non_homogeneous_activity non_homogeneous_activity_control_rd non_homogeneous_activity_hanangata

run: all
	mpirun -np 4 /non_homogeneous_activity

.PHONY: clean all run

clean:
	rm -f *.o *~ non_homogeneous_activity

