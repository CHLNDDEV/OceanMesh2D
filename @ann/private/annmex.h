#ifndef _ANNMEX_H_
#define _ANNMEX_H_

#include "mex.h"


#ifdef USE64BITS
    #include <tmwtypes.h>
    typedef uint64_T pointer_t;
    #define mxPOINTER_CLASS mxUINT64_CLASS
#else
    typedef unsigned long int pointer_t;
    #define mxPOINTER_CLASS mxUINT32_CLASS
#endif // USE64BITS
typedef unsigned int index_t;


#include "ANN.h"

 
#define CLASS_INTEGRITY_CHECK   0x12345678
#define BUCKET_SIZE             1
#define KD_SPLIT_RULE           ANN_KD_SUGGEST

/* use matlab's memory manager */
#include "my_mex_mem.h" 

// this enum affects definitions in private/modes.m 
typedef enum {
    OPEN = 1,
    CLOSE = 2,
    KSEARCH = 3,
    PRISEARCH=4,
    FRSEARCH=5,
} ann_mex_mode_t;

/* wrapper class */
class ann_mex_t {
    private:
        unsigned long int m_clintck;    // a flag verifying memoreys  integrity
        ANNkd_tree*       m_pTree;      
        index_t           m_idim;
        index_t           m_inpoints;
        ANNpointArray     pa;           // allocate point array - kd tree does not do it.
        
    public:
        ann_mex_t(const mxArray* mxPts);    // construct ann from a matrix
               
        ~ann_mex_t();                       // destruct -
        
        void annkSearch(					// approx k near neighbor search
            const mxArray*	mxQ,			// query point(s)
            int				k,				// number of near neighbors to return
            mxArray**		mxIdx,			// nearest neighbor array (modified)
            mxArray**       mxDst,			// dist to near neighbors (modified)
            double			eps=0.0);		// error bound
        
        void annkPriSearch( 				// priority k near neighbor search
            const mxArray*	mxQ,			// query point(s)
            int				k,				// number of near neighbors to return
            mxArray**		mxIdx,			// nearest neighbor array (modified)
            mxArray**       mxDst,			// dist to near neighbors (modified)
            double			eps=0.0);		// error bound
        
        int annkFRSearch(					// approx fixed-radius kNN search
            const mxArray*	mxQ,			// query point (only one at each time)
            ANNdist         sqRad,			// squared radius of query ball
            int				k,				// number of near neighbors to return
            mxArray**		mxIdx,			// nearest neighbor array (modified)
            mxArray**       mxDst,			// dist to near neighbors (modified)
            double			eps=0.0);		// error bound
    
        
        bool IsGood() { return m_clintck == CLASS_INTEGRITY_CHECK; };
        
        int  CordBitSiz() { return sizeof(ANNcoord); } ;    // get inside info on the data type ann uses
        bool isCordInt()  { ANNcoord tmp; tmp = 0.5; return tmp != 0.5; };
};
    
    
template<class T>
void GetScalar(const mxArray* x, T& scalar)
{
    if ( mxGetNumberOfElements(x) != 1 )
        mexErrMsgIdAndTxt("annmex:GetScalar","input is not a scalar!");
    void *p = mxGetData(x);
    switch (mxGetClassID(x)) {
        case mxCHAR_CLASS:
            scalar = *(char*)p;
            break;
        case mxDOUBLE_CLASS:
            scalar = *(double*)p;
            break;
        case mxSINGLE_CLASS:
            scalar = *(float*)p;
            break;
        case mxINT8_CLASS:
            scalar = *(char*)p;
            break;
        case mxUINT8_CLASS:
            scalar = *(unsigned char*)p;
            break;
        case mxINT16_CLASS:
            scalar = *(short*)p;
            break;
        case mxUINT16_CLASS:
            scalar = *(unsigned short*)p;
            break;
        case mxINT32_CLASS:
            scalar = *(int*)p;
            break;
        case mxUINT32_CLASS:
            scalar = *(unsigned int*)p;
            break;
#ifdef USE64BITS
        case mxINT64_CLASS:
            scalar = *(int64_T*)p;
            break;
        case mxUINT64_CLASS:
            scalar = *(uint64_T*)p;
            break;
#endif /* 64 bits machines */
            default:
                mexErrMsgIdAndTxt("annmex:GetScalar","unsupported data type");
    }
}


#endif // _ANNMEX_H_
