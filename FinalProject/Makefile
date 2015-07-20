CC = gcc
CFLAGS = -g -Wall -pedantic -Werror -std=c99
.c.o:  ; $(CC) -c $(CFLAGS) $<

OBJ = 	helper.o\
      	init.o\
      	boundary_val.o\
      	uvp.o\
      	main.o\
      	visual.o\
		sor.o\
		surface.o\


all:  $(OBJ)
	$(CC) $(CFLAGS) -o sim $(OBJ)  -lm

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	rm $(OBJ)

helper.o      : helper.h
init.o        : helper.h init.h surface.h
boundary_val.o: helper.h boundary_val.h
uvp.o         : helper.h uvp.h
visual.o      : helper.h surface.h
sor.o	      : sor.h
surface.o     : surface.h 

main.o        : helper.h init.h boundary_val.h uvp.h visual.h sor.h surface.h
