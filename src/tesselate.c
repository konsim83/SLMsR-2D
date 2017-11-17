/*****************************************************************************/
/*                                                                           */
/*  (tricall.c)                                                              */
/*                                                                           */
/*  Example program that demonstrates how to call Triangle.                  */
/*                                                                           */
/*  Accompanies Triangle Version 1.6                                         */
/*  July 19, 1996                                                            */
/*                                                                           */
/*  This file is placed in the public domain (but the file that it calls     */
/*  is still copyrighted!) by                                                */
/*  Jonathan Richard Shewchuk                                                */
/*  2360 Woolsey #H                                                          */
/*  Berkeley, California  94705-1927                                         */
/*  jrs@cs.berkeley.edu                                                      */
/*                                                                           */
/*****************************************************************************/

/* If SINGLE is defined when triangle.o is compiled, it should also be       */
/*   defined here.  If not, it should not be defined here.                   */

/* #define SINGLE */

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

int const Debug = 0;
int const  print_flag = 0;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allocate.h"
#include "triangle.h"
#include "tesselate.h"



/*****************************************************************************/
/*                                                                           */
/*  report()   Print the input or output.                                    */
/*                                                                           */
/*****************************************************************************/
void report(io, markers, reporttriangles, reportneighbors, reportsegments,
            reportedges, reportnorms)
struct triangulateio *io;
int markers;
int reporttriangles;
int reportneighbors;
int reportsegments;
int reportedges;
int reportnorms;
{
    int i, j;
    
    for (i = 0; i < io->numberofpoints; i++)
    {
	printf("Point %4d:", i);
	for (j = 0; j < 2; j++)
	{
	    printf("  %.3f", io->pointlist[i * 2 + j]);
	}
	if (io->numberofpointattributes > 0)
	{
	    printf("   attributes");
	}
	for (j = 0; j < io->numberofpointattributes; j++)
	{
	    printf("  %.3f",
		   io->pointattributelist[i * io->numberofpointattributes + j]);
	}
	if (markers)
	{
	    printf("   marker %3d\n", io->pointmarkerlist[i]);
	}
	else
	{
	    printf("\n");
	}
    }
    printf("\n");
  
    if (reporttriangles || reportneighbors)
    {
	for (i = 0; i < io->numberoftriangles; i++)
	{
	    if (reporttriangles)
	    {
		printf("Triangle %4d    has points:   ", i);
		for (j = 0; j < io->numberofcorners; j++)
		{
		    printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
		}
		if (io->numberoftriangleattributes > 0)
		{
		    printf("   attributes");
		}
		for (j = 0; j < io->numberoftriangleattributes; j++)
		{
		    printf("  %.3f", io->triangleattributelist[i *
							       io->numberoftriangleattributes + j]);
		}
		printf("\n");
	    }
	    if (reportneighbors)
	    {
		printf("Triangle %4d    has neighbors:", i);
		for (j = 0; j < 3; j++)
		{
		    printf("  %4d", io->neighborlist[i * 3 + j]);
		}
		printf("\n");
	    }
	}
	printf("\n");
    }

    if (reportsegments)
    {
	for (i = 0; i < io->numberofsegments; i++)
	{
	    printf("Segment %4d points:", i);
	    for (j = 0; j < 2; j++)
	    {
		printf("  %4d", io->segmentlist[i * 2 + j]);
	    }
	    if (markers)
	    {
		printf("   marker %d\n", io->segmentmarkerlist[i]);
	    }
	    else
	    {
		printf("\n");
	    }
	}
	printf("\n");
    }

    if (reportedges)
    {
	for (i = 0; i < io->numberofedges; i++)
	{
	    printf("Edge %4d points:", i);
	    for (j = 0; j < 2; j++)
	    {
		printf("  %4d", io->edgelist[i * 2 + j]);
	    }
	    if (reportnorms && (io->edgelist[i * 2 + 1] == -1))
	    {
		for (j = 0; j < 2; j++)
		{
		    printf("  %.3f", io->normlist[i * 2 + j]);
		}
	    }
	    if (markers)
	    {
		printf("   marker %d\n", io->edgemarkerlist[i]);
	    }
	    else
	    {
		printf("\n");
	    }
	}
	printf("\n");
    }
}



/*****************************************************************************/
/*                                                                           */
/*  ()   Create a mesh of the unit square.                                       */
/*                                                                           */
/*****************************************************************************/
my_triangulation tesselate_unit_square(int n_edges_per_segment, char* switches)
{
    struct triangulateio in, out;

    my_triangulation output;
  
    int i;  
    REAL h = 1.0 / n_edges_per_segment;


    in.numberofpoints = 4 * n_edges_per_segment;
    in.numberofsegments = in.numberofpoints;
    in.numberofpointattributes = 0;
    in.numberofholes = 0;
    in.numberofregions = 0; 
  
    /* Define input points. */
    in.pointlist = real_alloc(in.numberofpoints * 2, "in.pointlist");
    // Edge 0
    for(i = 0; i < n_edges_per_segment; i++)
    {
	in.pointlist[2*i] = i*h;
	in.pointlist[2*i+1] = 0;
    }
    // Edge 1
    for(i = 0; i < n_edges_per_segment; i++)
    {
	in.pointlist[2*i+2*n_edges_per_segment] = 1;
	in.pointlist[2*i+1+2*n_edges_per_segment] = h*i;
    }
    // Edge 2
    for(i = 0; i < n_edges_per_segment; i++)
    {
	in.pointlist[2*i+4*n_edges_per_segment] = 1 - i*h;
	in.pointlist[2*i+1+4*n_edges_per_segment] = 1;
    }
    // Edge 3
    for(i = 0; i < n_edges_per_segment; i++)
    {
	in.pointlist[2*i+6*n_edges_per_segment] = 0;
	in.pointlist[2*i+1+6*n_edges_per_segment] = 1 - i*h;
    }


    in.pointattributelist = real_alloc(0, "in.pointattributelist");
    in.pointmarkerlist = int_alloc(0, "in.pointmarkerlist");

    
    in.segmentlist = int_alloc(in.numberofpoints * 2, "in.segmentlist");
    in.segmentmarkerlist = int_alloc(0, "in.segmentmarkerlist");
    for(i = 0; i < in.numberofsegments; i++)
    {
	in.segmentlist[2*i] = (i+1) % in.numberofpoints;
	in.segmentlist[2*i+1] = i;
    }



    if (print_flag != 0)
    {
	printf("\n\n---------------   Points   ---------------\n");
	for(i = 0; i < in.numberofpoints; i++)
	{
	    printf("Point index   %4d   with coordinates   ", i);
	    printf("(%.3f, %.3f).\n", in.pointlist[2*i], in.pointlist[2*i+1]);
	}
	printf("---------------   End points   ---------------\n\n");

	printf("\n\n---------------   Segments   ---------------\n");
	for(i = 0; i < in.numberofsegments; i++)
	{
	    printf("Segement index   %4d   has points with indices   (%4d, %4d).\n", i, in.segmentlist[2*i], in.segmentlist[2*i+1]);
	}
	printf("---------------   End segments   ---------------\n\n");
    }
  
    
    out.pointlist = real_alloc(0, "out.pointlist");
    out.pointattributelist = real_alloc(0, "out.pointattributelist");
    out.pointmarkerlist = int_alloc(0,"out.pointmarkerlist") ;
    out.trianglelist = int_alloc(0,"out.trianglelist");
    out.triangleattributelist = real_alloc(0, "out.triangleattributelist");
    out.neighborlist = int_alloc(0,"out.neighborlist");
    out.segmentlist = int_alloc(0,"out.segmentlist");
    out.segmentmarkerlist = int_alloc(0,"out.segmentmarkerlist");
    out.edgelist = int_alloc(0,"out.edgelist");
    out.edgemarkerlist = int_alloc(0,"out.edgemarkerlist");
    

    /* Triangulate the points. */
    printf("\n\n\n+++++++++++++++++++   Meshing of input   +++++++++++++++++++\n\n");
    triangulate(switches, &in, &out, (struct triangulateio *) NULL);
    printf("+++++++++++++++++++   End meshing   +++++++++++++++++++\n\n\n");


    if (print_flag != 0)
    {
	printf("\n\n\n===================   Initial triangulation   ===================\n");
	report(&out, 1, 1, 1, 0, 1, 0);
	printf("Number of Points:   %4d\n", out.numberofpoints);
	printf("Number of Edges:   %4d\n", out.numberofedges);
	printf("Number of Triangles:   %4d\n", out.numberoftriangles);
	printf("===================   End initial triangulation   ===================\n\n\n");
    }
  
   
    /* Release the memory. */
    in.pointlist = real_free(in.pointlist, "in.pointlist");
    in.pointattributelist = real_free(in.pointattributelist, "in.pointattributelist");
    in.pointmarkerlist = int_free(in.pointmarkerlist, "in.pointmarkerlist");
    in.segmentlist = int_free(in.segmentlist, "in.segmentlist");
    in.segmentmarkerlist = int_free(in.segmentmarkerlist, "in.segmentmarkerlist");
    
    
    output.pointlist  = out.pointlist;
    output.pointattributelist = out.pointattributelist;
    output.pointmarkerlist = out.pointmarkerlist;
    output.numberofpoints = out.numberofpoints;
    output.numberofpointattributes = out.numberofpointattributes;

    output.trianglelist = out.trianglelist;
    output.triangleattributelist = out.triangleattributelist;
    output.trianglearealist = out.trianglearealist;
    output.neighborlist = out.neighborlist;
    output.numberoftriangles = out.numberoftriangles;
    output.numberofcorners = out.numberofcorners;
    output.numberoftriangleattributes = out.numberoftriangleattributes;

    output.segmentlist = out.segmentlist;
    output.segmentmarkerlist = out.segmentmarkerlist;
    output.numberofsegments = out.numberofsegments;

    output.holelist = out.holelist;
    output.numberofholes = out.numberofholes;

    output.regionlist = out.regionlist;
    output.numberofregions = out.numberofregions;

    output.edgelist = out.edgelist;
    output.edgemarkerlist = out.edgemarkerlist;
    output.normlist = out.normlist;
    output.numberofedges = out.numberofedges;
    
    
    return output;
}



/*****************************************************************************/
/*                                                                           */
/*  ()   Create a mesh of the unit simplex.                                       */
/*                                                                           */
/*****************************************************************************/
my_triangulation tesselate_unit_simplex(REAL dummy, char* switches)
{
    struct triangulateio in, out;
    my_triangulation output;
    int i;  

    in.numberofpoints = 3;
    in.numberofsegments = in.numberofpoints;
    in.numberofpointattributes = 0;
    in.numberofholes = 0;
    in.numberofregions = 0; 
  
    /* Define input points. */
    in.pointlist = real_alloc(in.numberofpoints * 2, "in.pointlist");
    in.pointlist[0] = 0;
    in.pointlist[1] = 0;
    in.pointlist[2] = 1;
    in.pointlist[3] = 0;
    in.pointlist[4] = 0;
    in.pointlist[5] = 1;



    in.pointattributelist = real_alloc(0, "in.pointattributelist");
    in.pointmarkerlist = int_alloc(0, "in.pointmarkerlist");

    
    in.segmentlist = int_alloc(in.numberofpoints, "in.segmentlist");
    in.segmentmarkerlist = int_alloc(0, "in.segmentmarkerlist");

    in.segmentlist[0] = 0;
    in.segmentlist[1] = 1;
    in.segmentlist[2] = 1;
    in.segmentlist[3] = 2;
    in.segmentlist[4] = 2;
    in.segmentlist[5] = 0;


    if (print_flag != 0)
    {
	printf("\n\n---------------   Points   ---------------\n");
	for(i = 0; i < in.numberofpoints; i++)
	{
	    printf("Point index   %4d   with coordinates   ", i);
	    printf("(%.3f, %.3f).\n", in.pointlist[2*i], in.pointlist[2*i+1]);
	}
	printf("---------------   End points   ---------------\n\n");

	printf("\n\n---------------   Segments   ---------------\n");
	for(i = 0; i < in.numberofsegments; i++)
	{
	    printf("Segement index   %4d   has points with indices   (%4d, %4d).\n", i, in.segmentlist[2*i], in.segmentlist[2*i+1]);
	}
	printf("---------------   End segments   ---------------\n\n");
    }
  
    
    out.pointlist = real_alloc(0, "out.pointlist");
    out.pointattributelist = real_alloc(0, "out.pointattributelist");
    out.pointmarkerlist = int_alloc(0,"out.pointmarkerlist") ;
    out.trianglelist = int_alloc(0,"out.trianglelist");
    out.triangleattributelist = real_alloc(0, "out.triangleattributelist");
    out.neighborlist = int_alloc(0,"out.neighborlist");
    out.segmentlist = int_alloc(0,"out.segmentlist");
    out.segmentmarkerlist = int_alloc(0,"out.segmentmarkerlist");
    out.edgelist = int_alloc(0,"out.edgelist");
    out.edgemarkerlist = int_alloc(0,"out.edgemarkerlist");
    

    /* Triangulate the points. */
    printf("\n\n\n+++++++++++++++++++   Meshing of input   +++++++++++++++++++\n\n");
    triangulate(switches, &in, &out, (struct triangulateio *) NULL);
    printf("+++++++++++++++++++   End meshing   +++++++++++++++++++\n\n\n");


    if (print_flag != 0)
    {
	printf("\n\n\n===================   Initial triangulation   ===================\n");
	report(&out, 1, 1, 1, 0, 1, 0);
	printf("Number of Points:   %4d\n", out.numberofpoints);
	printf("Number of Edges:   %4d\n", out.numberofedges);
	printf("Number of Triangles:   %4d\n", out.numberoftriangles);
	printf("===================   End initial triangulation   ===================\n\n\n");
    }
  
   
    /* Release the memory. */
    in.pointlist = real_free(in.pointlist, "in.pointlist");
    in.pointattributelist = real_free(in.pointattributelist, "in.pointattributelist");
    in.pointmarkerlist = int_free(in.pointmarkerlist, "in.pointmarkerlist");
    in.segmentlist = int_free(in.segmentlist, "in.segmentlist");
    in.segmentmarkerlist = int_free(in.segmentmarkerlist, "in.segmentmarkerlist");
    
    
    output.pointlist  = out.pointlist;
    output.pointattributelist = out.pointattributelist;
    output.pointmarkerlist = out.pointmarkerlist;
    output.numberofpoints = out.numberofpoints;
    output.numberofpointattributes = out.numberofpointattributes;

    output.trianglelist = out.trianglelist;
    output.triangleattributelist = out.triangleattributelist;
    output.trianglearealist = out.trianglearealist;
    output.neighborlist = out.neighborlist;
    output.numberoftriangles = out.numberoftriangles;
    output.numberofcorners = out.numberofcorners;
    output.numberoftriangleattributes = out.numberoftriangleattributes;

    output.segmentlist = out.segmentlist;
    output.segmentmarkerlist = out.segmentmarkerlist;
    output.numberofsegments = out.numberofsegments;

    output.holelist = out.holelist;
    output.numberofholes = out.numberofholes;

    output.regionlist = out.regionlist;
    output.numberofregions = out.numberofregions;

    output.edgelist = out.edgelist;
    output.edgemarkerlist = out.edgemarkerlist;
    output.normlist = out.normlist;
    output.numberofedges = out.numberofedges;
    
    
    return output;
}


/*****************************************************************************/


my_triangulation tesselate_unit_simplex_uniform_edges(int n_edges_per_segment, char* switches)
{
    struct triangulateio in, out;
    my_triangulation output;
    int i;  

    REAL h = 1.0 / n_edges_per_segment;

    in.numberofpoints = 3 * n_edges_per_segment;
    in.numberofsegments = in.numberofpoints;
    in.numberofpointattributes = 0;
    in.numberofholes = 0;
    in.numberofregions = 0; 
  
    /* Define input points. */
    in.pointlist = real_alloc(in.numberofpoints * 2, "in.pointlist");
    // Edge 0
    for(i = 0; i < n_edges_per_segment; i++)
    {
	in.pointlist[2*i] = i*h;
	in.pointlist[2*i+1] = 0;
    }
    // Edge 1
    for(i = 0; i < n_edges_per_segment; i++)
    {
	in.pointlist[2*i+2*n_edges_per_segment] = 1-i*h;
	in.pointlist[2*i+1+2*n_edges_per_segment] = h*i;
    }
    // Edge 2
    for(i = 0; i < n_edges_per_segment; i++)
    {
	in.pointlist[2*i+4*n_edges_per_segment] = 0;
	in.pointlist[2*i+1+4*n_edges_per_segment] = 1 - i*h;
    }
    

    in.pointattributelist = real_alloc(0, "in.pointattributelist");
    in.pointmarkerlist = int_alloc(0, "in.pointmarkerlist");

    
    in.segmentlist = int_alloc(in.numberofpoints * 2, "in.segmentlist");
    in.segmentmarkerlist = int_alloc(0, "in.segmentmarkerlist");
    for(i = 0; i < in.numberofsegments; i++)
    {
	in.segmentlist[2*i] = (i+1) % in.numberofpoints;
	in.segmentlist[2*i+1] = i;
    }


    if (print_flag != 0)
    {
	printf("\n\n---------------   Points   ---------------\n");
	for(i = 0; i < in.numberofpoints; i++)
	{
	    printf("Point index   %4d   with coordinates   ", i);
	    printf("(%.3f, %.3f).\n", in.pointlist[2*i], in.pointlist[2*i+1]);
	}
	printf("---------------   End points   ---------------\n\n");

	printf("\n\n---------------   Segments   ---------------\n");
	for(i = 0; i < in.numberofsegments; i++)
	{
	    printf("Segement index   %4d   has points with indices   (%4d, %4d).\n", i, in.segmentlist[2*i], in.segmentlist[2*i+1]);
	}
	printf("---------------   End segments   ---------------\n\n");
    }
  
    
    out.pointlist = real_alloc(0, "out.pointlist");
    out.pointattributelist = real_alloc(0, "out.pointattributelist");
    out.pointmarkerlist = int_alloc(0,"out.pointmarkerlist") ;
    out.trianglelist = int_alloc(0,"out.trianglelist");
    out.triangleattributelist = real_alloc(0, "out.triangleattributelist");
    out.neighborlist = int_alloc(0,"out.neighborlist");
    out.segmentlist = int_alloc(0,"out.segmentlist");
    out.segmentmarkerlist = int_alloc(0,"out.segmentmarkerlist");
    out.edgelist = int_alloc(0,"out.edgelist");
    out.edgemarkerlist = int_alloc(0,"out.edgemarkerlist");
    

    /* Triangulate the points. */
    printf("\n\n\n+++++++++++++++++++   Meshing of input   +++++++++++++++++++\n\n");
    triangulate(switches, &in, &out, (struct triangulateio *) NULL);
    printf("+++++++++++++++++++   End meshing   +++++++++++++++++++\n\n\n");


    if (print_flag != 0)
    {
	printf("\n\n\n===================   Initial triangulation   ===================\n");
	report(&out, 1, 1, 1, 0, 1, 0);
	printf("Number of Points:   %4d\n", out.numberofpoints);
	printf("Number of Edges:   %4d\n", out.numberofedges);
	printf("Number of Triangles:   %4d\n", out.numberoftriangles);
	printf("===================   End initial triangulation   ===================\n\n\n");
    }
  
   
    /* Release the memory. */
    in.pointlist = real_free(in.pointlist, "in.pointlist");
    in.pointattributelist = real_free(in.pointattributelist, "in.pointattributelist");
    in.pointmarkerlist = int_free(in.pointmarkerlist, "in.pointmarkerlist");
    in.segmentlist = int_free(in.segmentlist, "in.segmentlist");
    in.segmentmarkerlist = int_free(in.segmentmarkerlist, "in.segmentmarkerlist");
    
    
    output.pointlist  = out.pointlist;
    output.pointattributelist = out.pointattributelist;
    output.pointmarkerlist = out.pointmarkerlist;
    output.numberofpoints = out.numberofpoints;
    output.numberofpointattributes = out.numberofpointattributes;

    output.trianglelist = out.trianglelist;
    output.triangleattributelist = out.triangleattributelist;
    output.trianglearealist = out.trianglearealist;
    output.neighborlist = out.neighborlist;
    output.numberoftriangles = out.numberoftriangles;
    output.numberofcorners = out.numberofcorners;
    output.numberoftriangleattributes = out.numberoftriangleattributes;

    output.segmentlist = out.segmentlist;
    output.segmentmarkerlist = out.segmentmarkerlist;
    output.numberofsegments = out.numberofsegments;

    output.holelist = out.holelist;
    output.numberofholes = out.numberofholes;

    output.regionlist = out.regionlist;
    output.numberofregions = out.numberofregions;

    output.edgelist = out.edgelist;
    output.edgemarkerlist = out.edgemarkerlist;
    output.normlist = out.normlist;
    output.numberofedges = out.numberofedges;
    
    
    return output;
}


/*****************************************************************************/


my_triangulation tesselate_triangle_uniform_edges(REAL a_x, REAL a_y,
						  REAL b_x, REAL b_y,
						  REAL c_x, REAL c_y,
						  int n_edges_per_segment, char* switches)
{
    struct triangulateio in, out;
    my_triangulation output;
    int i;  

    REAL h = 1.0 / n_edges_per_segment;

    in.numberofpoints = 3 * n_edges_per_segment;
    in.numberofsegments = in.numberofpoints;
    in.numberofpointattributes = 0;
    in.numberofholes = 0;
    in.numberofregions = 0; 
  
    /* Define input points. */
    in.pointlist = real_alloc(in.numberofpoints * 2, "in.pointlist");
    // Edge 0
    for(i = 0; i < n_edges_per_segment; i++)
    {
	in.pointlist[2*i] = a_x*(1 - i*h) + b_x*i*h;
	in.pointlist[2*i+1] = a_y*(1 - i*h) + b_y*i*h;
    }
    // Edge 1
    for(i = 0; i < n_edges_per_segment; i++)
    {
	in.pointlist[2*i+2*n_edges_per_segment] = b_x*(1 - i*h) + c_x*i*h;
	in.pointlist[2*i+1+2*n_edges_per_segment] = b_y*(1 - i*h) + c_y*i*h;
    }
    // Edge 2
    for(i = 0; i < n_edges_per_segment; i++)
    {
	in.pointlist[2*i+4*n_edges_per_segment] = c_x*(1 - i*h) + a_x*i*h;
	in.pointlist[2*i+1+4*n_edges_per_segment] = c_y*(1 - i*h) + a_y*i*h;
    }

    in.pointattributelist = real_alloc(0, "in.pointattributelist");
    in.pointmarkerlist = int_alloc(0, "in.pointmarkerlist");

    
    in.segmentlist = int_alloc(in.numberofpoints * 2, "in.segmentlist");
    in.segmentmarkerlist = int_alloc(0, "in.segmentmarkerlist");
    for(i = 0; i < in.numberofsegments; i++)
    {
	in.segmentlist[2*i] = (i+1) % in.numberofpoints;
	in.segmentlist[2*i+1] = i;
    }


    if (print_flag != 0)
    {
	printf("\n\n---------------   Points   ---------------\n");
	for(i = 0; i < in.numberofpoints; i++)
	{
	    printf("Point index   %4d   with coordinates   ", i);
	    printf("(%.3f, %.3f).\n", in.pointlist[2*i], in.pointlist[2*i+1]);
	}
	printf("---------------   End points   ---------------\n\n");

	printf("\n\n---------------   Segments   ---------------\n");
	for(i = 0; i < in.numberofsegments; i++)
	{
	    printf("Segement index   %4d   has points with indices   (%4d, %4d).\n", i, in.segmentlist[2*i], in.segmentlist[2*i+1]);
	}
	printf("---------------   End segments   ---------------\n\n");
    }
  
    
    out.pointlist = real_alloc(0, "out.pointlist");
    out.pointattributelist = real_alloc(0, "out.pointattributelist");
    out.pointmarkerlist = int_alloc(0,"out.pointmarkerlist") ;
    out.trianglelist = int_alloc(0,"out.trianglelist");
    out.triangleattributelist = real_alloc(0, "out.triangleattributelist");
    out.neighborlist = int_alloc(0,"out.neighborlist");
    out.segmentlist = int_alloc(0,"out.segmentlist");
    out.segmentmarkerlist = int_alloc(0,"out.segmentmarkerlist");
    out.edgelist = int_alloc(0,"out.edgelist");
    out.edgemarkerlist = int_alloc(0,"out.edgemarkerlist");

    /* Triangulate the points. */
    printf("\n\n\n+++++++++++++++++++   Meshing of input   +++++++++++++++++++\n\n");
    triangulate(switches, &in, &out, (struct triangulateio *) NULL);
    printf("+++++++++++++++++++   End meshing   +++++++++++++++++++\n\n\n");


    if (print_flag != 0)
    {
	printf("\n\n\n===================   Initial triangulation   ===================\n");
	report(&out, 1, 1, 1, 0, 1, 0);
	printf("Number of Points:   %4d\n", out.numberofpoints);
	printf("Number of Edges:   %4d\n", out.numberofedges);
	printf("Number of Triangles:   %4d\n", out.numberoftriangles);
	printf("===================   End initial triangulation   ===================\n\n\n");
    }
  
   
    /* Release the memory. */
    in.pointlist = real_free(in.pointlist, "in.pointlist");
    in.pointattributelist = real_free(in.pointattributelist, "in.pointattributelist");
    in.pointmarkerlist = int_free(in.pointmarkerlist, "in.pointmarkerlist");
    in.segmentlist = int_free(in.segmentlist, "in.segmentlist");
    in.segmentmarkerlist = int_free(in.segmentmarkerlist, "in.segmentmarkerlist");
    
    
    output.pointlist  = out.pointlist;
    output.pointattributelist = out.pointattributelist;
    output.pointmarkerlist = out.pointmarkerlist;
    output.numberofpoints = out.numberofpoints;
    output.numberofpointattributes = out.numberofpointattributes;

    output.trianglelist = out.trianglelist;
    output.triangleattributelist = out.triangleattributelist;
    output.trianglearealist = out.trianglearealist;
    output.neighborlist = out.neighborlist;
    output.numberoftriangles = out.numberoftriangles;
    output.numberofcorners = out.numberofcorners;
    output.numberoftriangleattributes = out.numberoftriangleattributes;

    output.segmentlist = out.segmentlist;
    output.segmentmarkerlist = out.segmentmarkerlist;
    output.numberofsegments = out.numberofsegments;

    output.holelist = out.holelist;
    output.numberofholes = out.numberofholes;

    output.regionlist = out.regionlist;
    output.numberofregions = out.numberofregions;

    output.edgelist = out.edgelist;
    output.edgemarkerlist = out.edgemarkerlist;
    output.normlist = out.normlist;
    output.numberofedges = out.numberofedges;
    
    
    return output;
}
