#ifndef TESSELATE_H
#define TESSELATE_H

#include "triangle.h"

typedef struct{
    REAL *pointlist;                                               /* In / out */
    REAL *pointattributelist;                                      /* In / out */
    int *pointmarkerlist;                                          /* In / out */
    int numberofpoints;                                            /* In / out */
    int numberofpointattributes;                                   /* In / out */

    int *trianglelist;                                             /* In / out */
    REAL *triangleattributelist;                                   /* In / out */
    REAL *trianglearealist;                                         /* In only */
    int *neighborlist;                                             /* Out only */
    int numberoftriangles;                                         /* In / out */
    int numberofcorners;                                           /* In / out */
    int numberoftriangleattributes;                                /* In / out */

    int *segmentlist;                                              /* In / out */
    int *segmentmarkerlist;                                        /* In / out */
    int numberofsegments;                                          /* In / out */

    REAL *holelist;                        /* In / pointer to array copied out */
    int numberofholes;                                      /* In / copied out */

    REAL *regionlist;                      /* In / pointer to array copied out */
    int numberofregions;                                    /* In / copied out */

    int *edgelist;                                                 /* Out only */
    int *edgemarkerlist;            /* Not used with Voronoi diagram; out only */
    REAL *normlist;                /* Used only with Voronoi diagram; out only */
    int numberofedges;                                             /* Out only */
} my_triangulation;

my_triangulation tesselate_unit_square(int, char*);
my_triangulation tesselate_unit_simplex(REAL, char*);
my_triangulation tesselate_unit_simplex_uniform_edges(int, char*);
my_triangulation tesselate_triangle_uniform_edges(REAL, REAL, REAL, REAL, REAL, REAL, int, char*);

#endif // TESSELATE_H
