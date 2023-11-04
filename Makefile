compile:
	cd dp; \
	gcc -fpic -c mfp.c; \
	gcc -shared -o mfp.so mfp.o; \
	gcc -fpic -c mfp_dynamic_computing.c; \
	gcc -shared -o mfp_dynamic_computing.so mfp_dynamic_computing.o; \
	gcc -fpic -c mfp_unrestricted_length.c; \
	gcc -shared -o mfp_unrestricted_length.so mfp_unrestricted_length.o

	