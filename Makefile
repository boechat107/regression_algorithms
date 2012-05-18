CXX = gcc
CFLAGS = -Wall -std=c99 -g

#all: svm-train svm-predict svm-scale
all: aakr-predict

#lib: svm.o
#	$(CXX) -shared -dynamiclib -Wl,-soname,libsvm.so.$(SHVER) svm.o -o libsvm.so.$(SHVER)
#

aakr-predict: aakr-predict.c aakr.o
	$(CXX) $(CFLAGS) aakr-predict.c aakr.o -o aakr-predict -lm

aakr.o: aakr.c aakr.h
	$(CXX) $(CFLAGS) -c aakr.c 

clean:
	rm -f aakr.o aakr-predict

mem:
	valgrind --leak-check=yes ./aakr-predict train.data test.data

#svm-predict: svm-predict.c svm.o
#	$(CXX) $(CFLAGS) svm-predict.c svm.o -o svm-predict -lm
#svm-train: svm-train.c svm.o
#	$(CXX) $(CFLAGS) svm-train.c svm.o -o svm-train -lm
#svm-scale: svm-scale.c
#	$(CXX) $(CFLAGS) svm-scale.c -o svm-scale
#svm.o: svm.cpp svm.h
#	$(CXX) $(CFLAGS) -c svm.cpp
#clean:
#	rm -f *~ svm.o svm-train svm-predict svm-scale libsvm.so.$(SHVER)
