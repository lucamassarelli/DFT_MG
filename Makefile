CC = gcc #select compiler

install: DFT_MG.c Hartree.c rt_nonfinite.c rtGetNaN.c rtGetInf.c
	$(CC) -c Hartree.c -o Hartree.o
	$(CC) -c rt_nonfinite.c -o rt_nonfinite.o
	$(CC) -c rtGetInf.c -o rtGetInf.o
	$(CC) -c rtGetNaN.c -o rtGetNaN.o
	ar rcs libHartree.a Hartree.o
	ar rcs librt_nonfinite.a rt_nonfinite.o
	ar rcs librtGetInf.a rtGetInf.o
	ar rcs librtGetNaN.a rtGetNaN.o

	$(CC) DFT_MG.c -o DFT_MG.o -L $(CURDIR) -lHartree -lrt_nonfinite -lrtGetInf -lrtGetNaN -lm

#clean installation
clean:
ifeq ($(OS),Windows_NT) 
	del *.a
	del Hartree.o
	del rt_nonfinite.o
	del rtGetInf.o
	del rtGetNaN.o
else
	rm Hartree.o
	rm rt_nonfinite.o
	rm rtGetInf.o
	rm rtGetNaN.o
	rm *.a
endif
