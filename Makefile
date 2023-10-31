compile:
	gcc -fpic -c mfp.c 
	gcc -shared -o mfp.so mfp.o