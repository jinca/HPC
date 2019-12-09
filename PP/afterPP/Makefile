SRC=MD.c control.c
OBJ=$(SRC:.c=.o)
CC=icc 
CFLAGS=-O3 -ipo -restrict


all: MD

MD: $(OBJ)
	$(CC) $(CFLAGS)  -o $@  $(OBJ) -lm


output.dat: MD input.dat
	./MD


clean:
	rm -f MD $(OBJ) 

$(OBJ) : coord.h Makefile


