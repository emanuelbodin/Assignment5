CFLAGS= -g -pg -Wall -Ofast

galsim: galsim.o
	gcc -o galsim galsim.o -lm

galsim.o: galsim.c 
	gcc $(CFLAGS) -c galsim.c

clean:
	rm -f ./galsim *.o
