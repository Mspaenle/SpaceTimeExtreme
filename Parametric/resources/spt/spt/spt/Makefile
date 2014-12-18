#!/bin/bash

CC=gcc

$CC  -I/Library/Frameworks/R.framework/Resources/include  -DNDEBUG  -I/usr/local/include -I/usr/local/include/freetype2 -I/opt/X11/include    -fPIC  -Wall -mtune=core2 -g -O2  -c main.c -o main.o -fopenmp
$CC -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/usr/local/lib -o main.so main.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation


gcc -I/Library/Frameworks/R.framework/Resources/include -o main.exe main.c timer.c spt.c 
