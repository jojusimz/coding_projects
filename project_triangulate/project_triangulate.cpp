#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <iterator>
#include <functional>
#include <cmath>
#include <omp.h>

//#include "interval.h"

using namespace std;
template < typename F> class f_x_y;
template < class s> class is_point_in_tri;

using namespace std;
template<class T>
class vertex
{

public:
//constructors
    vertex ( int vid, const T& xx, const T& yy, const T& zz=0);
    vertex ( int vid, const vector<T> &coords);
    vertex (const vertex<T> &vv);  // copy constructor
//destructor
   // ~vertex();
//file reading and writing
    void read_node_file (const string& vertexfile);
    void write_node_file( const string& vertexfile);

// utility member functions
    void v_set_id(const int& id);
    void v_set_x(const T& xx);
    void v_set_y(const T& yy);
    void v_set_z(const T& zz);
    int v_get_id()const;
    int v_set_id()const;
    T v_get_x()const;
    T v_get_y()const;
    T v_get_z()const;

// operator overloading
    template<typename s>
    friend ostream& operator <<(std::ostream& out,vertex< s >& vv);            // output streaming
    void operator=(const vertex<T>& vv);                                       // assignment operator = overloaded to allow VERTEX = VERTEX operation
    vertex<T> operator()(const vertex<T>& vv);
private:
    int v_id;
    vector<T> coordinates;                                                      // z coordinate
    int v_dimension;
};
//.........................................................................................................................................
template<class T>
class triangle
{
public:
    triangle (int id, vertex<T> *v0, vertex<T> *v1, vertex<T> *v2): t_id(id),vertex_0(v0),vertex_1(v1),vertex_2(v2) {};  // takes pointer to vertex.
    triangle(const triangle<T> &tt);                                                                                     //copy constructor  - NOT DEFINED SINCE TIME CONSTRAINT ON PROJECT DID NOT PERMIT IT

    //destructor
   // ~triangle();
    //utility member functions get and set
    int t_get_id()const;
    vertex<T> t_get_v0()const;
    T t_get_v0_x()const;
    T t_get_v0_y()const;
    T t_get_v0_z()const;
    vertex<T> t_get_v1()const;
    T t_get_v1_x()const;
    T t_get_v1_y()const;
    T t_get_v1_z()const;
    vertex<T> t_get_v2()const;
    T t_get_v2_x()const;
    T t_get_v2_y()const;
    T t_get_v2_z()const;
    void t_set(int& id, vertex<T> *v0, vertex<T> *v1, vertex<T> *v2);

    //calculate area
    T calc_area_t ()const;
    T calc_area_t (const T &x0, const T &y0, const T &x1, const T &y1, const T &x2, const T &y2)const;

    // calculate integral - two methods provided: constant approximation and linear interpolation
    template <typename s > T calc_f_integral1 ( const f_x_y<s> &ff_x_y  )const;
    template <typename s > T calc_f_integral2 ( const f_x_y<s> &ff_x_y  )const;           //

    //compute circumcentre e.t.c
    vector<T> circumcentre_x_y()const;
    bool t_is_in_circumcircle_(const T &x, const T& y)const;
    bool t_radius_compare(const T& x,const T& y,const T& cx,const T& cy,const T& r)const;
    //bool is_point_in_tri(const T &x, const T &y)const;

    // member functions to find triangle neighbour
    bool is_a_neighbour(const triangle<T> &tri_)const;
    void set_neighbour( triangle<T> *tri_);
    vector<triangle <T> *> t_get_neighbours()const;  // returns pointer to triangle neighbours

    //operator overload for output streaming
    template<typename s>
    friend ostream& operator <<(ostream& out,triangle< s >& tri);
    triangle<T> operator()(const triangle<T>& tt)const;


private:
    int t_id;
    vertex<T>* vertex_0;
    vertex<T>* vertex_1;
    vertex<T>* vertex_2;
    vector  < triangle <T> *> neighbours;                  //points to triangles neighbours
};
//...................................................................................................................................
template<class T>
class t_mesh
{
public:
    //constructor
    t_mesh( ) {};                        // default constructor
    t_mesh(const t_mesh<T> & mesh);      // copy constructor     - NOT DEFINED SINCE TIME CONSTRAINT ON PROJECT DID NOT PERMIT

    //destructor
    //~t_mesh();

    void create_delaunay( const T & xx, const T & yy);
    // NOT DEFINED SINCE TIME CONSTRAINT OF PROJECT DID NOT PERMIT THIS
    //a declaration of the interface to a member function that re-triangulates
    // i.e. given a set of points xx and yy, find the triangle that contains the points
    // and re-triangulate using bowyer watson to make delaunay.
    // this function was not implemented since it was not a requirement for the exercise.

    //file reading and writing
    void read_triangulation_file (const string& trianglefile);                                //reads from file
    void write_triangulation_file (const string& trianglefile);                                 //writes to file
    void print_neighbour();                                                                  // prints triangle + its neighbours id
    void find_neighbours();                                                                  // finds the neighbours
    bool is_delaunay();                                                                       //checks for delaunay - the method chosen involved checking the neighbours of each triangle
    // it should be noted that approximation errors interfere with the method
    // incorporating exact maths / interval maths component could help manage
    // this problem.
    //
    bool find_t(const int & t_id)const;                                                     // given a triangles id and finds if it is in a triangulation - SEEMS Redundant
    int  find_t(const T & xx, const T&  yy);                                                 // given a set of points find the triangle it falls inside, returns triangle id.

    template<typename s> T mesh_calc_f_integral1( const f_x_y<s> &ff_x_y )const;
    template<typename s> T mesh_calc_f_integral2( const f_x_y<s> &ff_x_y )const;

    //operator overload for output streaming
    template<typename s>
    friend ostream& operator <<(ostream& out,t_mesh< s >& tt);
    //assignment operator overloading
    void operator=(const t_mesh<T>& mesh);
    t_mesh<T> operator()(const t_mesh<T> &mesh);
private:
    vector< triangle <T> >  t_i;
    vector< vertex <T> >    v_i;
    int t_dimension;
    int t_attributes_size;
    int neighbours_have_been_set;
};



// definitions.......................................................................................
template<class T>
vertex<T>::vertex(int vid, const T& xx, const T& yy, const T& zz)
    : v_id(vid)
{
    coordinates.push_back(xx);
    coordinates.push_back(yy);
    coordinates.push_back(zz);

}

template<class T>
vertex<T>::vertex (const vertex<T> &vv)  // copy constructor
{
    v_id = vv.v_id;
    coordinates = vv.coordinates;

}

template<class T>
vertex<T>::vertex ( int vid, const vector<T> &coords)
{
    v_id = vid;
    coordinates.push_back(coords[0]);
    coordinates.push_back(coords[1]);
    if (coords.size() == 3) coordinates.push_back(coords[2]);

}


/*template<class T>
vertex<T>::~vertex()
{

}*/
template<class T>
vertex<T> vertex<T>::operator()(const vertex<T>& vv)                                       // assignment operator = overloaded to allow VERTEX = VERTEX operation
{
    vertex<T> temp (vv);
    return vv;

}

template<class T>
void vertex<T>::operator=(const vertex<T>& vv)
{
    v_id= vv.v_id;
    coordinates = vv.coordinates;

}

template<class T>
inline void vertex<T>::v_set_x(const T& xx)
{
    coordinates[0]=xx;
}

template<class T>
inline void vertex<T>::v_set_y(const T& yy)
{
    if(v_dimension>2) coordinates[1]=yy;
}

template<class T>
inline void vertex<T>::v_set_z(const T& zz)
{
    if(v_dimension==3) coordinates[2]=zz;
}

template<class T>
inline void vertex<T>::v_set_id(const int& id)
{
    this->v_id = id;
}

template<class T>
inline int vertex<T>:: v_get_id()const
{
    return this->v_id;
}

template<class T>
inline T vertex<T>:: v_get_x()const
{
    return this->coordinates[0];
}

template<class T>
inline T vertex<T>:: v_get_y()const
{
    return this->coordinates[1];
}

template<class T>
inline T vertex<T>:: v_get_z()const
{
    try
    {

        if (v_dimension > 3)
            throw 1;
    }
    catch (int e)
    {
        cout<< "Error in accessing..... no z coordinates" << e <<endl;
        return -99999999;                                                     // design rule : this would be replaced by some message value that indicates a problem
    }
    return this->coordinates[2];

}
template <typename s>
ostream& operator <<(ostream& out,vertex< s >& vv)                        // output stream for vertex prints: id x y z
{
    out<<vv.v_id<<" "<< vv.coordinates[0] <<" "<< vv.coordinates[1] <<" ";

    if( vv.coordinates.size()==3 ) out << vv.coordinates[2]<<endl;

    else out << endl;

    return out;
}
//............................................................................................................................................
//copy constructor
template<class T>
triangle<T>::triangle(const triangle<T> &tt)
{
    t_id = tt.t_id;
    vertex_0 =tt.vertex_0;
    vertex_1 =tt.vertex_1;
    vertex_2 =tt.vertex_2;
}
/*
//triangle destructor
template< class T>
triangle<T>::~triangle()
{
//delete [] neighbours;
//delete [] t_attributes;
}
*/
template<class T>
triangle<T> triangle<T>::operator()(const triangle<T>& tt)const
{
    return triangle<T> (tt) ;

}

//triangle get and set functions
template<class T>
void triangle<T>::t_set(int &id, vertex<T> *v0, vertex<T> *v1, vertex<T> *v2)
{
    this->t_id = id;
    this->vertex_0=v0;
    this->vertex_1=v1;
    this->vertex_2=v2;
}
template<class T>
inline T triangle<T>:: t_get_v0_x() const              //return x coordinate of Vertex V0 of triangle
{
    return this->vertex_0->v_get_x();
}

template<class T>
inline T triangle<T>:: t_get_v0_y() const              //return y coordinate of Vertex V0 of triangle
{
    return this->vertex_0->v_get_y();
}

template<class T>
inline T triangle<T>:: t_get_v0_z()const                //return z coordinate of Vertex V0 of triangle
{
    return this->vertex_0->v_get_z();
}

template<class T>
inline T triangle<T>:: t_get_v1_x()const              //return X coordinate of Vertex V1 of triangle
{
    return this->vertex_1->v_get_x();
}

template<class T>
inline T triangle<T>:: t_get_v1_y()const              //return Y coordinate of Vertex V1 of triangle
{
    return this->vertex_1->v_get_y();
}

template<class T>
inline T triangle<T>:: t_get_v1_z()const              //return Z coordinate of Vertex V1 of triangle
{
    return this->vertex_1->v_get_z();
}

template<class T>
inline T triangle<T>:: t_get_v2_x()const              //return X coordinate of Vertex V2 of triangle
{
    return this->vertex_2->v_get_x();
}

template<class T>
inline T triangle<T>:: t_get_v2_y()const              //return Y coordinate of Vertex V2 of triangle
{
    return this->vertex_2->v_get_y();
}

template<class T>
inline T triangle<T>:: t_get_v2_z()const              //return Z coordinate of Vertex V2 of triangle
{
    return this->vertex_2->v_get_z();
}

template<class T>
inline int triangle<T>::t_get_id()const              //return ID of triangle
{
    return t_id;
}

template<class T>
inline vertex<T> triangle<T>::t_get_v0()const              //return vertex V0 of triangle
{
    return *vertex_0;
}

template<class T>
inline vertex<T> triangle<T>::t_get_v1()const              //return vertex V1 of triangle
{
    return *vertex_1;
}

template<class T>
inline vertex<T> triangle<T>::t_get_v2()const              //return vertex V2 of triangle
{
    return *vertex_2;
}

template<class T>
T triangle<T>::calc_area_t ()const                        //return area of triangle
{

    T x0 = this->t_get_v0_x();
    T x1 = this->t_get_v1_x();
    T x2 = this->t_get_v2_x();
    T y0 = this->t_get_v0_y();
    T y1 = this->t_get_v1_y();
    T y2 = this->t_get_v2_y();

    return ( fabs(0.5* (x0 *(y1-y2) +  x1*(y2-y0) + x2*(y0-y1) )  ));
}

template<class T>
T triangle<T>::calc_area_t (const T& x0, const T& y0, const T& x1, const T& y1, const T& x2, const T& y2)const   // returns area of triangle when a set of vertices is provided.
{

    return ( fabs(0.5* (x0 *(y1-y2) +  x1*(y2-y0) + x2*(y0-y1) )  ));
}

template<class T>
vector<T> triangle<T>::circumcentre_x_y( ) const                       // returns the circumcentre of the triangle
{
    vector<T> o_x_y;
    T x0 = this->t_get_v0_x();
    T x1 = this->t_get_v1_x();
    T x2 = this->t_get_v2_x();
    T y0 = this->t_get_v0_y();
    T y1 = this->t_get_v1_y();
    T y2 = this->t_get_v2_y();

    T a =( x0 * x0) + ( y0 * y0 ) - ( x1 * x1 ) - ( y1 * y1 );
    T b =( x1 * x1) + ( y1 * y1 ) - ( x2 * x2 ) - ( y2 * y2 );
    T aa = x0 - x1;
    T bb = y0 - y1;
    T cc = x1 - x2;
    T dd = y1 - y2;


    try
    {
        if((bb*cc-aa*dd)==0)
        {
            throw "CANNOT divide by zero: error computing circumcentre.....check values!!";
        }
        T ox = 0.5*( b *bb - dd * a) / ( bb * cc - aa *dd) ;
        T oy = 0.5*( a *cc - aa * b) / ( bb * cc - aa *dd) ;
        vector<T> o_x_y(2,0);
        o_x_y.push_back(ox);
        o_x_y.push_back(oy);
    }                               // checks for division by 0
    catch ( const char *msg)
    {
        cerr<<msg<<endl;
        return o_x_y;
    }
    return (o_x_y);
}

template < typename T>
template < typename s >
T triangle<T>::calc_f_integral1( const f_x_y<s> &ff_x_y ) const                              // Integral 1: constant value approximation
{
    vector<T> temp (2,0);

    T result (0.0);

    temp = this->circumcentre_x_y();
    T x0 = this->t_get_v0_x();
    T x1 = this->t_get_v1_x();
    T x2 = this->t_get_v2_x();
    T y0 = this->t_get_v0_y();
    T y1 = this->t_get_v1_y();
    T y2 = this->t_get_v2_y();

    T area1 =  calc_area_t(temp[0],temp[1],x1,y1,x2,y2);
    T area2 =  calc_area_t(temp[0],temp[1],x0,y0,x2,y2);
    T area3 =  calc_area_t(temp[0],temp[1],x1,y1,x0,y0);

    result = area1*ff_x_y(x0,y0) + area2*ff_x_y(x1,y1) + area3*ff_x_y(x2,y2);

    return (result);
}

template < typename T>
template < typename s >
T triangle<T>::calc_f_integral2 ( const f_x_y<s> &ff_x_y  ) const             //integral 2: Linear interpolation approximation
{
    vector<T> temp (2,0);

    T result (0.0);

    T x0 = this->t_get_v0_x();
    T x1 = this->t_get_v1_x();
    T x2 = this->t_get_v2_x();
    T y0 = this->t_get_v0_y();
    T y1 = this->t_get_v1_y();
    T y2 = this->t_get_v2_y();

    T area =  calc_area_t(x0,y0,x1,y1,x2,y2);

    result =  area*0.333333*( ff_x_y(x0,y0) + ff_x_y(x1,y1) + ff_x_y(x2,y2) );

    return (result);
}


template<typename s>
ostream& operator <<(std::ostream& out,triangle< s >& tri)                        // output stream for triangle prints : id , v0_id, v1_id : v2_id
{

    out<<tri.t_id<<" "<< tri.vertex_0->v_get_id()<<" "<< tri.vertex_1->v_get_id()<<" "<< tri.vertex_2->v_get_id()<<endl;

}

template<class T>
vector<triangle<T> *> triangle<T>::t_get_neighbours() const                  // get neighbours of triangle - returns vector of pointers to neighbours
{
    return this->neighbours;
}

template<class T>
void triangle<T>::set_neighbour( triangle<T> *tri_)                          // sets triangle neighbours.
{
    this->neighbours.push_back(tri_);
}
template<class T>
bool triangle<T>::t_radius_compare(const T& x,const T& y,const T& cx,const T& cy,const T& r)const             //compares the radius - used to find if point is in circumcircle
{
    return ((sqrt((x-cx)*(x-cx)+(y-cy)*(y-cy))-r)<-0.5);
}

template<class T>
bool triangle<T>:: t_is_in_circumcircle_(const T& x,const T& y)const                    // checks if point is in circumcircle
{

    T x0 = this->vertex_0->v_get_x();
    T x1 = this->vertex_1->v_get_x();
    T x2 = this->vertex_2->v_get_x();
    T y0 = this->vertex_0->v_get_y();
    T y1 = this->vertex_1->v_get_y();
    T y2 = this->vertex_2->v_get_y();

    T m1 = (x0-x1)/(y1-y0);
    T m2 = (x0-x2)/(y2-y0);
    T b1 = ((y0+y1)*0.5) - 0.5*m1*(x0+x1);
    T b2 = ((y0+y2)*0.5) - 0.5*m2*(x0+x2);
    T xx = (b2-b1)/(m1-m2);
    T yy = m1*xx + b1;

    return t_radius_compare(x,y,xx,yy,sqrt((xx-x0)*(xx-x0)+(yy-y0)*(yy-y0)));
}

template<class T>
bool triangle<T>:: is_a_neighbour(const triangle<T> &tri_)const            // returns true if tri_ is a neighbour
{

    bool check0 ( ( this->vertex_0->v_get_id()== tri_.vertex_0->v_get_id()) ||
                  ( this->vertex_0->v_get_id()== tri_.vertex_1->v_get_id()) ||
                  ( this->vertex_0->v_get_id()== tri_.vertex_2->v_get_id()) ) ;

    bool check1 ( ( this->vertex_1->v_get_id()== tri_.vertex_0->v_get_id()) ||
                  ( this->vertex_1->v_get_id()== tri_.vertex_1->v_get_id()) ||
                  ( this->vertex_1->v_get_id()== tri_.vertex_2->v_get_id()) ) ;

    bool check2 ( ( this->vertex_2->v_get_id()== tri_.vertex_0->v_get_id()) ||
                  ( this->vertex_2->v_get_id()== tri_.vertex_1->v_get_id()) ||
                  ( this->vertex_2->v_get_id()== tri_.vertex_2->v_get_id()) ) ;

    return ((check0&&check1) || (check0&&check2) || (check1&&check2));
}

//.................................................................................................................................................
/*
template < typename T >
t_mesh<T>::~t_mesh()
{
    //delete[] t_i;
    //delete[] v_i;
}
*/
template < typename T >
t_mesh<T>::t_mesh(const t_mesh<T> & mesh)
{
    t_i = mesh.t_i;
    v_i = mesh.v_i;
    t_dimension = mesh.t_dimension;
    neighbours_have_been_set = mesh.neighbours_have_been_set;

}

template < typename T >
void::t_mesh<T>::create_delaunay(const T& xx,const T&yy)                    //NOT DEFINED -  for future implementation
{


}

template < typename T >
void t_mesh<T>::operator=(const t_mesh<T>& mesh)
{
    t_i = mesh.t_i;
    v_i = mesh.v_i;
    t_dimension = mesh.t_dimension;
    neighbours_have_been_set = mesh.neighbours_have_been_set;

}

template < typename T >
t_mesh<T> t_mesh<T>::operator()(const t_mesh<T> &mesh)
{
    return t_mesh<T> (mesh);
}


template < typename T >
template < typename s >
T t_mesh<T>::mesh_calc_f_integral1( const f_x_y<s> &ff_x_y ) const                            // Integral 1: constant value approximation
{
    typename vector< triangle<T> >::const_iterator first = this->t_i.begin() ;
    typename vector< triangle<T> >::const_iterator last = this->t_i.end();
    T sum(0.);
    for ( first; first!=last ; first++)
    {

        sum +=first->calc_f_integral1(ff_x_y);
    }
    return sum;

}

template < typename T >
template < typename s >
T t_mesh<T>::mesh_calc_f_integral2( const f_x_y<s> &ff_x_y ) const                              // Integral 1: constant value approximation
{
    typename vector< triangle<T> >::const_iterator first = this->t_i.begin() ;
    typename vector< triangle<T> >::const_iterator last = this->t_i.end();
    T sum(0.);
    for ( first; first!=last ; first++)
    {
        cout<<sum<<endl;
        sum +=first->calc_f_integral2(ff_x_y);
    }
    return sum;

}

template<class T>
void t_mesh<T>::write_triangulation_file (const string& trianglefile)       //creates triangulation file
{
    ofstream out_file;

    out_file.open((trianglefile+"").c_str());
    out_file<<this->v_i.size()<<" "<<this->t_dimension<<" "<<endl;
    int end_ = this->v_i.size();

    for(int ii = 0 ; ii< end_; ++ii)
    {
        out_file<<v_i[ii];
    }

    out_file<<this->t_i.size()<<" "<<this->t_dimension<<endl;

    end_ = this->t_i.size();

    for(int ii = 0 ; ii< end_; ++ii)
    {
        out_file<<t_i[ii];
    }

}


template<class T>
void t_mesh<T>::read_triangulation_file (const string& trianglefile)                      //reads from triangulation file into t_mesh
{
    vector<T> whole_file;
    int no_of_nodes(0), dimension(0), attributes(0);

    ifstream in_file;
    try
    {
        in_file.open((trianglefile+"").c_str());
        exception e;
        if(!in_file) throw e;
    }
    catch(exception e)
    {
        cout<<"file not found!!";
    }

    in_file>>no_of_nodes>>dimension>>attributes;
    copy(istream_iterator <T>(in_file), istream_iterator<T>(),inserter(whole_file,whole_file.begin()));

    this->t_dimension = dimension;

    neighbours_have_been_set = 0;

    vector<T> temp((no_of_nodes)*(dimension+1),0);


    int iterate=dimension+1,i(0);

    int size_no_file = (dimension+1)*(no_of_nodes);                                         // computes the total data in the node_file

    swap_ranges(whole_file.begin(),whole_file.begin()+size_no_file,  temp.begin());


    if( dimension==2 )                                                          //checks dimension and populates appropriately
    {
        for(int ii=0; ii<no_of_nodes; ii++)
        {

            vertex<T> temp_v ((int)temp[i],temp[i+1],temp[i+2] );

            this->v_i.push_back(temp_v);

            i+=iterate;
        }
    }


    else if (dimension==3)                                                          // checks dimension and populates appropriately
    {

        for(int ii=0; ii<no_of_nodes; ii++)
        {

            vertex<T> temp_v ((int)temp[i],temp[i+1],temp[i+2],temp[i+3] );

            v_i.push_back(temp_v);

            //cout<<v_i[ii];

            i+=iterate;
        }
    }

    if( whole_file.size() > size_no_file )                      // ensures triangle vertices is defined in the file and also allows the odd case of wanting to
        // read from the vertices.node into
        // the mesh vertices.
    {
        int t_no = whole_file[size_no_file], t_dimensions = whole_file[size_no_file+1], t_attributes =whole_file[size_no_file+2];

        int size_t_file =t_no*(1+t_dimensions + t_attributes );

        temp.resize(size_t_file,0);

        swap_ranges(whole_file.begin()+size_no_file+3,whole_file.begin()+(size_no_file+size_t_file), temp.begin());

        iterate=t_dimensions+t_attributes+1;                                                // variable used to skip the reading the attributes

        i=0;


        for(int ii=0; ii<t_no; ii++)
        {

            triangle<T> t_temp ( (int)temp[i], &v_i[temp[i+1]],  &v_i[temp[i+2]], &v_i[temp[i+3]] );

            t_i.push_back(t_temp);

            //cout<<t_i[ii]<<endl;

            i+=iterate;
        }
    }
}


template<class T>
int t_mesh<T>::find_t(const T& xx, const T & yy )                   //returns id of first triangle with point
{
    int sizee=t_i.size();

    typename vector< triangle<T> >::iterator it = find_if(this->t_i.begin(),this->t_i.end(),is_point_in_tri<T>(xx,yy));

    if( it != this->t_i.end())
    {
        return it->t_get_id();
    }
    else return -1;                                                      //returns negative int if point isn't in mesh.
}

template <typename s>
ostream& operator <<(std::ostream& out,t_mesh< s >& tt)
{
    vector<int> temp(3,0);
    int i_end = tt.t_i.size();
    for(int i=0; i < i_end; ++i)
    {
        out<<tt.t_i[i]<<endl;
    }
    return out;
}

template<class T>
void t_mesh<T>::print_neighbour()                 // prints id of each triangle and its neighbours
{
    if(neighbours_have_been_set==0)
    {
        cout<<"neighbours have not been set"<<endl;
    }
    else
    {
   //omp_set_dynamic(0);

    int i, var(2);

   // cout<<"enter number of threads" <<endl;

   // cin>>var;

    omp_set_num_threads(var);
    #pragma omp parallel private(i) default(shared)
    {
        #pragma omp for //  schedule(static)
        for (i = 0; i < this->t_i.size(); i++)
        {
            cout<< t_i[i];

            for(int j=0 ; j < t_i[i].t_get_neighbours().size() ; j++)
            {
                cout<<"neighbour :"<<j<<*( this->t_i[i].t_get_neighbours()[j] )<< endl;
            }
        }
    }
    }
}


template<class T>
void t_mesh<T>::find_neighbours()                                                            // finds neighbours of each triangle in a mesh
{
    omp_set_dynamic(0);

    int ii, var(2);

   // cout<<"enter number of threads" <<endl;

    //cin>>var;

    omp_set_num_threads(var);

    #pragma omp parallel private(ii) default(shared)
    {
        #pragma omp for //num_threads(4) //  schedule(static)
        for( ii = 0; ii < t_i.size() ; ii++)                                    //outer for loop searches the triangles and tests the vertices of the neighbours
            {
               for ( int jj=0; jj < t_i.size(); jj++)
        {
            if(t_i[ii].t_get_neighbours().size()>3) break;                                     // stops the inner loop once max neighbours have been found

                if( ( ii!=jj )&& ( this->t_i[ii].is_a_neighbour( this->t_i[jj] )))           // inner for loop searches the mesh to find the triangles that share vertices
                {
                    t_i[ii].set_neighbour(&t_i[jj]);

                }
                                                                                            // stops inner loop once 3 neighbours have been identified.
        }
    }
    }
    neighbours_have_been_set = 1;                                                            // flags true that the neighbours have been set - used in delaunay function

}

template<class T>
bool t_mesh<T>::is_delaunay()
{
    bool check(1),check0(1),check1(1),check2(1);

    if(!neighbours_have_been_set)                                                    //checks if neighbours have been set. if false then sets neighbours and calls functions
    {
        this->find_neighbours();
        this->is_delaunay();
    }
    else
    {


    T x0, x1, x2, y0, y1, y2, xx0, xx1, xx2, yy0, yy1, yy2;
    int ii,total(0);

    omp_set_dynamic(0);

    int var(2);

    //cout<<"enter number of threads" <<endl;

    //cin>>var;

    omp_set_num_threads(var);

    #pragma omp parallel private(ii) default(shared)
    {
        #pragma omp for //  schedule(static)
        for( ii = 0; ii < t_i.size() ; ii++)                                    //outer for loop searches the triangles and tests the vertices of the neighbours
                                                                                    // to check if it falls inside circumcircle.
        {
            for ( int jj=0; jj < t_i[ii].t_get_neighbours().size(); jj++)           // searches vertices of neighbours. maximum loop size = 3
                                                                                    // big O( 3n )
            {

                xx0 = this->t_i[ii].t_get_neighbours()[jj]->t_get_v0_x();
                yy0 = this->t_i[ii].t_get_neighbours()[jj]->t_get_v0_y();

                check0  = ! ( t_i[ii].t_is_in_circumcircle_(xx0, yy0) );

                xx1 = this->t_i[ii].t_get_neighbours()[jj]->t_get_v1_x();
                yy1 = this->t_i[ii].t_get_neighbours()[jj]->t_get_v1_y();

                check1  = ! ( t_i[ii].t_is_in_circumcircle_(xx1, yy1) );

                xx2 = this->t_i[ii].t_get_neighbours()[jj]->t_get_v2_x();
                yy2 = this->t_i[ii].t_get_neighbours()[jj]->t_get_v2_y();

                check2  = ! ( t_i[ii].t_is_in_circumcircle_(xx2, yy2) );

                //total = check0+check1+check2;                                          //  IGNORE - used in debugging!!
                if ((check0+check1+check2)<3 )   check = 0;                              // sum should always be three.
            }
        }
     }
    }
    return check;
}


template<typename T> void input_function(T &a)
{
    int check(0), count(0);

    do
    {
        check=0;

        cout<<endl<<endl<<" Enter 1 ,2 ,3 or 4 to check for speed up of member functions" <<endl;

        cin>>a;

        if( a>5)
        {
        check=1;
        count++;
        cout<<" Error: you have entered the wrong value " <<endl;
        }

    }while((check==1)&&(count<4));
}


int main(int argc, char* argv[])
{
    t_mesh<double> test;

    test.read_triangulation_file("triangulation#2.tri");

    //timer start

    int enter(0);

    input_function(enter);

    if(enter == 1){

        cout<< " .....openmp speedup tests on find_neighbours member function......."<<endl;

        struct timespec Start_timer, End_timer;

        clock_gettime(CLOCK_REALTIME, &Start_timer);        // START CLOCK //timer stop

        test.find_neighbours();

        clock_gettime(CLOCK_REALTIME, &End_timer);       //STOP CLOCK

       //MEASURE TIME ELAPSED
        double time_elapsed = ( End_timer.tv_sec - Start_timer.tv_sec )+ (double)((End_timer.tv_nsec - Start_timer.tv_nsec)*1E-9);

        cout<<"\n\n clock_gettime wall time= "<<time_elapsed<<endl;     //RETRIEVE MEASURED TIME

    }
    else if(enter == 2){

        cout<< " ..... openmp test on is_delaunay member function......."<<endl;

        struct timespec Start_timer, End_timer;

        clock_gettime(CLOCK_REALTIME, &Start_timer);        // START CLOCK //timer stop

        test.is_delaunay();

        clock_gettime(CLOCK_REALTIME, &End_timer);       //STOP CLOCK

       //MEASURE TIME ELAPSED
        double time_elapsed = ( End_timer.tv_sec - Start_timer.tv_sec )+ (double)((End_timer.tv_nsec - Start_timer.tv_nsec)*1E-9);

        cout<<"\n\n clock_gettime wall time= "<<time_elapsed<<endl;     //RETRIEVE MEASURED TIME

    }

    else if(enter == 3){

        cout<< " ..... openmp test on print_neighbour member function......."<<endl;

        struct timespec Start_timer, End_timer;

        clock_gettime(CLOCK_REALTIME, &Start_timer);        // START CLOCK //timer stop

        test.find_neighbours();

        test.print_neighbour();

        clock_gettime(CLOCK_REALTIME, &End_timer);       //STOP CLOCK

       //MEASURE TIME ELAPSED
        double time_elapsed = ( End_timer.tv_sec - Start_timer.tv_sec )+ (double)((End_timer.tv_nsec - Start_timer.tv_nsec)*1E-9);

        cout<<"\n\n clock_gettime wall time= "<<time_elapsed<<endl;     //RETRIEVE MEASURED TIME

    }

    return 0;
}




