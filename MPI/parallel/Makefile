MF= Makefile

# for cirrus
CC=	mpicc

CFLAGS=	-cc=icc	-std=c99	-O3
LFLAGS=	-lm

EXE=	imagenewp768x768_1000

INC= \
	pgmio.h

SRC= \
	imagenewp768x768_1000_src.c \
	pgmio.c


#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(OBJ):	$(INC)

$(EXE):	$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

clean:
	rm -f $(OBJ) $(EXE) core
