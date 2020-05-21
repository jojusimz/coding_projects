
# DELAUNAY TRIGULATION PROJECT 

**Triangulation** is the decomposition of 2D geometrical planar object P into a set of triangles. It involves finding a set of triangles with pairwise non-intersecting interiors whose union geometrically approximates P.

**Delaunay triangulation**: Given a set input discrete points in a plane the Delaunay triangulation discretizes such that no point is inside the circumcircle of any triangle. Delaunay triangulations maximize the minimum angle of all the angles of the triangles in the triangulation and thus avoids the occurence of sliver triangles.

https://www.ics.uci.edu/~eppstein/pubs/BerEpp-CEG-95.pdf

**Images of triangulation**

![image](https://user-images.githubusercontent.com/60849864/82563624-4cc6e600-9b6f-11ea-8601-298b7188bbf3.png)

# Documentation

See the documentation file for a more detailed description of the project.

# Project Overview

The aim of this project was to implement a code which verifies if a set of triangulation points from a " .TRI file " and " .NODE file" fulfilled the delaunay triangulation criteria. 

The code was to be enabled for the following:
 
 1. File input and output i.e. can read an input tri.file -> manipulate -> output a tri.file.

 2. Check if a triangulation is Delaunay.
 
 3.	Performing the integral of a function  f(x,y) over a triangulation domain.
 
 4.	Queries of the form “ given x,y coordinates locate the triangle containing it "
 
 5.	Figure 4.1 below illustrates the data structure employed to encapsulate the triangulation. 
 
 ![image](https://user-images.githubusercontent.com/60849864/82563539-2c972700-9b6f-11ea-83b7-c0e5cb84e0c2.png)
 
 
 # Features employed in this project
 
 * OOP:
 
 1. Template functions and classes.
 
 2. Operator overloading.
 
 3. Functors
 
* OpenMP OPTIMIZATION

 # COMPILATION
 
 The code was developed to be compiled in a UNIX envrionment using the instructions:
 
          g++ -O3-o main.exe project_triangulation.cpp

          main.exe
 
