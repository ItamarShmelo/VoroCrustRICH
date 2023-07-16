#ifndef _HILBERT_TYPES_H
#define _HILBERT_TYPES_H

typedef double coord_t;
typedef unsigned long int hilbert_index_t;

typedef struct _3DPoint
{
    coord_t x;
    coord_t y;
    coord_t z;
} _3DPoint;

typedef struct _2DPoint
{
    coord_t x;
    coord_t y;
} _2DPoint;

#endif // _HILBERT_TYPES_H