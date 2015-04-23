all:
	gcc dp.c -fno-builtin -O3 -lm -o dp
	./dp_test
