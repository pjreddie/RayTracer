CC=g++
LDFLAGS=-O2 -lm -lpthread -L/usr/X11R6/lib -lm -lpthread -lX11

RayTracer: RayTracer.o
	$(CC) $(LDFLAGS) RayTracer.o -o RayTracer
RayTracer.o: RayTracer.cpp CImg/CImg.h
	$(CC) -c RayTracer.cpp
