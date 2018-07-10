#ifndef _MY_MEX_MEM_H_
#define _MY_MEX_MEM_H_

#include "mex.h"

/* memory allocations - redirect to MATLAB memory menager */
void* operator new(size_t size);
void* operator new[](size_t size);
void operator delete(void* ptr);
void operator delete[](void* ptr);

#endif // _MY_MEX_MEM_H_
