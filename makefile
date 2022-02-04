TARGET = SCFT_SPHERE
SRC_DIR = ./src
SRC = $(wildcard $(SRC_DIR)/*.c)
CC = gcc
CFLAGS = -std=c99 -O3 -Wall  -I ./include
LFLAGS = -lshtns_omp -lfftw3_omp -lfftw3 -lm -lgsl -lgslcblas -lstdc++ -fopenmp

$(TARGET):$(SRC)
		$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)