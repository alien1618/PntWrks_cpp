#!/bin/bash
echo "************************"
echo "PntWrks"
echo "************************"
echo "Generating makefile..."
mkdir logs
mkdir out
mkdir pics
mkdir obj
rm -r out/*
rm -r pics/*
rm obj/main.o
rm run
rm a.out
rm makefile.mak
g++ src/makefile.cpp
./a.out
echo "Compiling source code..."
make -f makefile.mak
echo "Running executable..."
./run
echo "Cleaning..."
rm run
rm makefile.mak
rm a.out
rm obj/main.o
echo "Script complete..."
