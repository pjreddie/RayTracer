CC=g++
LDFLAGS=-O2 -lm -lpthread -L/usr/X11R6/lib -lm -lpthread -lX11

RayTracer: Ray.o
	$(CC) $(LDFLAGS) Ray.o -o RayTracer
Ray.o: Ray.cpp CImg/CImg.h
	$(CC) -c Ray.cpp
