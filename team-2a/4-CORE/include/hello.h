#ifndef _HELLO
#define _HELLO

#include <stdio.h>
#include <math.h>
#include <cstring>
#include "core_cuda.h"

__global__ void helloFromGPU(void);

#endif // _HELLO