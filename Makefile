compile:
	cd src; \
	gcc -fpic -c -lpthread mfp.c; \
	gcc -shared -o mfp.so mfp.o

	