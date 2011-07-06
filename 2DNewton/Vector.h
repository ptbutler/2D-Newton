/*
  File: 	 Vector.h
  Author: 	 J.K. Johnstone
  Created:	 11 March 1998
  Last Modified: 22 July 2009
  Purpose: 	 A foundational class, consisting of arrays and vectors.
                 The essential components of the user-defined array are 
		 bounds-checking and dynamic allocation.
		 Almost every other class uses Array, and every piece of geometric 
		 software uses V3f.
		 Functions included are an adhoc collection reflecting the needs
		 of other algorithms (think lazy evaluation): no claims to completeness.
  Dependence:    Independent of other libraries.
  Status:        
  Structure:     Array --> Vec --> (Vec1, Vec2, Vec3, Vec4, Vec5)
*/

#ifndef _VECTOR_
#define _VECTOR_

#include <assert.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

/************************************************************************
       Array = an array of objects.
       Advantages: n and vector are one unit; bounds-checking.
************************************************************************/

template <typename Type>
class Array {

    template <class T> friend ostream& operator<< (ostream &os, Array<T> &a);
    template <class T> friend T        fabs (Array<T> a);

  public:
  
    Array()      { n=0; }
    Array(int N) { v = new Type[n=N]; }
    Array(int N, Type initval) { v=new Type[n=N]; for(int i=0;i<n;i++) v[i]=initval; }
    Array(Type *a, int N);
    Array(Array &a);
    virtual ~Array() { if (n>0) delete [] v; }
    void   create (int N) 
              { if (n>0) delete [] v; v = new Type[n=N]; }  
    void   create (int N, Type initval) 
              { if (n>0) delete [] v; v = new Type[n=N]; for (int i=0; i<n; i++) v[i]=initval; }
    void   create (Type *a, int N) 
              { if (n>0) delete [] v; v = new Type[n=N]; for (int i=0; i<n; i++) v[i] = a[i]; }
    void   create (Array &a, int N);	// copy first N elements of this array
    void   create (Array &a) 
              { if (n>0) delete [] v; v = new Type[n=a.n]; for (int i=0; i<n; i++) v[i]=a[i]; }
    void   allocate (int N) { create(N); } // necessary as unambiguous base version of create
    void   allocateJustLike (Array &a, int nLevel);
    void   grow (int n);
    void   shrink (int truesize);

    int    getn() { return n; }
    Array& operator=(Array&);
    int    operator==(Array&);
    int    operator!=(Array&);
    Array  operator-(Array &a);
    void   operator+=(Array &a);
    Type&  operator[](int index) { assert(index>=0 && index<n); return v[index]; }
    void   clear()      { for (int i=0; i<n; i++) v[i] = 0; } // only for number arrays
    void   set (Type a) { for (int i=0; i<n; i++) v[i] = a; } // only for number arrays

    void   append (Array<Type> &a, Array<Type> &b);
    void   append (Type &a);
    int    binarySearch (int i, int j, Type x);	
    void   bubbleSort (Array<int> &sortIndex);
    void   bubbleSortReverse (Array<int> &sortIndex);
    void   deleteDuplicate();
    void   deleteApproxDuplicate (float eps);
    void   deleteApproxEq (float val, float eps);
    void   filterOut (Array<Type> &w, Array<Type> &filterArr, float eps);
    int	   findInterval (Type x);
    void   insert (Type x, int index);
    void   insert (Array<Type> &a, int index);
    void   insertionSort (Array<int> &sortIndex);
    Type   mean();
    void   mean (Type &m);
    Type   median();
    void   merge (Array<Type> &a, Array<Type> &b);
    void   merge (Array<Type> &a, Array<int> &aindex,
		  Array<Type> &b, Array<int> &bindex,
		  Array<int> &sortIndex);
    void   merge (Array<Type> &a);
    void   mergeSort (Array<int> &sortIndex);
    int    nSubset (int nPtPerSubset, int overlap=0);
    void   removeDuplicate (int nDuplicate, Array<int> &duplicate);
    void   reverseIt();
    void   rotateToFront (int index);
    Type   select (int i, int j, int k);
    void   subsetSize (int nPtPerSubset, int overlap, Array<int> &subsetsize);
    void   swap (int i, int j);
    void   Union (Type x, int &size);
    Type   vecmin (int &index);
    Type   vecmax (int &index);

  protected:
  
    int    n;
    Type  *v;
    int    countDuplicate (Array<int> &duplicate);
    int    countApproxDuplicate (Array<int> &duplicate, float eps);
    int    countApproxEq (float val, Array<int> &duplicate, float eps);
    int    partition (int i, int j, Type pivot);
};

/************************************************************************
       Vec = a point in n-space, not just any array of values.
       Inherits from Array. 
************************************************************************/

template <class Type>
class Vec : public Array<Type> {
  using Array<Type>::n;
                 // using Array<Type>::v;  // don't use, confuses Quaternion's inheritance of v

  public:
  
    Vec() : Array<Type>() {}
    Vec(int dim) : Array<Type> (dim) {}
    Vec(Type *vec, int dim) : Array<Type> (vec,dim) {}
    Vec(Vec &vec) : Array<Type> (vec.v, vec.n) {}
    int   getdimension() { return n; }
    void  operator+= (Vec& rhs) { for (int i=0;i<n; i++) this->v[i] += rhs[i];}
    void  operator-= (Vec& rhs) { for (int i=0;i<n; i++) this->v[i] -= rhs[i];}
    void  operator*= (Type rhs) { for (int i=0;i<n; i++) this->v[i] *= rhs; }
    void  operator/= (Type rhs) { for (int i=0;i<n; i++) this->v[i] /= rhs; }
    Vec   operator+ (Vec& rhs) { Vec result(*this); result+=rhs; return result; } // NOT TESTED
    Vec   operator- (Vec& rhs) { Vec result(*this); result-=rhs; return result; } // NOT TESTED
    Vec&  operator= (Vec& rhs) { for (int i = 0; i < n; i++) this->v[i] = rhs[i]; return *this;}
    float operator* (Vec &vec) { return (dot(vec)); }
    int   operator< (Vec& rhs) {return(this->v[0] < rhs[0]);} // allow vec comp on 1st coord

    float angle();
    float angle (Vec&);
    float angleUnit (Vec&);
    float angle (Vec&, Vec&);
    int   approxeq (Vec& a, float eps);
    float dist (Vec&);
    float dot (Vec&);
    float length();
    void  matrixMult (float **matrix, Vec &product);
    void  matrixMult (float **matrix);
    void  midPt (Vec &a, Vec &b);
    void  normalize();
    void  swap (Vec &vec) { Vec temp=vec; vec=*this; *this=temp; }
    void  swap (int i, int j);
};

/************************************************************************
************************************************************************/

template <class Type>
class Vec1 : public Vec<Type> {

  public:

    Vec1() : Vec<Type>(1) {}
    Vec1(Type *vec) : Vec<Type>(vec,1) {}
    Vec1(Vec1 &vec) : Vec<Type>(vec) {}
    Vec1(Type a)    : Vec<Type>(1) { this->v[0]=a; }
};

/************************************************************************
************************************************************************/

template <class Type>
class Vec2 : public Vec<Type> {

  public:

    Vec2() : Vec<Type>(2) {}
    Vec2(Type *vec) : Vec<Type>(vec,2) {}
    Vec2(Vec2 &vec) : Vec<Type>(vec) {}
    Vec2(Type a, Type b)  : Vec<Type>(2) { create(a,b); }
    void  create (Type &a, Type &b) 	 { this->v[0]=a; this->v[1]=b; }
    void  create (Vec2 &vec) 		 { create (vec[0], vec[1]); }

    float angleFromCenter (Vec2 &center);
    float angleFromX();
    int   contains (float t);
    int   intersectInterval (Vec2 &b, Vec2 &intersection);
    int   leftTurn (Vec2 &P, Vec2 &Q);
    void  rotate (Type c, Type s);
};

/************************************************************************
************************************************************************/

template <class Type>
class Vec3 : public Vec<Type> {

  public:

    Vec3() : Vec<Type>(3) {}
    Vec3(Type *vec) : Vec<Type>(vec,3) {}
    Vec3(Vec3 &vec) : Vec<Type>(vec) {}
    Vec3(Type a, Type b, Type c) : Vec<Type>(3) { create(a,b,c); }
    void  create(Type a, Type b, Type c) { this->v[0]=a; this->v[1]=b; this->v[2]=c; }
    void  create(Vec3 &vec)              { create (vec[0], vec[1], vec[2]); }
    Vec3  operator+ (Vec3& rhs) { Vec3 result(*this); result+=rhs; return result; }// NOT TESTED
    Vec3  operator- (Vec3& rhs) { Vec3 result(*this); result-=rhs; return result; }// NOT TESTED

    void  cross (Vec3 &a, Vec3 &b);
    float distFromPlane (Vec3 &ptOnPlane, Vec3 &unitNormal, int signedDist=0);
    void  rotToCoord (int coordAxis, float **rot);
    void  unitNorm (Vec3 &a, Vec3 &b, Vec3 &c);
    void  xrot (float nDeg);
    void  yrot (float nDeg);
    void  zrot (float nDeg);
};

/************************************************************************
************************************************************************/

template <class Type>
class Vec4 : public Vec<Type> {

  public:

    Vec4() : Vec<Type>(4) {}
    Vec4(Type *vec) : Vec<Type>(vec,4) {}
    Vec4(Vec4 &vec) : Vec<Type>(vec) {}
    Vec4(Type a, Type b, Type c, Type d) : Vec<Type>(4) { create(a,b,c,d); }
    void create(Type a, Type b, Type c, Type d) 
      { this->v[0]=a; this->v[1]=b; this->v[2]=c; this->v[3]=d; }
    void  create(Vec4 &vec) { create (vec[0], vec[1], vec[2], vec[3]); }
};

/************************************************************************
************************************************************************/

template <class Type>
class Vec5 : public Vec<Type> {

  public:

    Vec5() : Vec<Type>(5) {}
    Vec5(Type *vec) : Vec<Type>(vec,5) {}
    Vec5(Vec5 &vec) : Vec<Type>(vec) {}
    Vec5(Type a, Type b, Type c, Type d, Type e) : Vec<Type>(5) { create(a,b,c,d,e); }
    void create(Type a, Type b, Type c, Type d, Type e)
      { this->v[0]=a; this->v[1]=b; this->v[2]=c; this->v[3]=d; this->v[4] = e; }
};

/************************************************************************
************************************************************************/

// two functions needed by rotToCoord
// (we want to avoid dependence on Miscellany, so we don't put there)

/*extern*/ float det (int n, float **a); 
/*extern*/ float subdet (int n, float **a, int row, int col);

/************************************************************************
************************************************************************/

typedef Array<int>    	    IntArr;
typedef Array<float>  	    FloatArr;
typedef Array<double> 	    DoubleArr;
typedef Array<IntArr>  	    IntArrArr;
typedef Array<FloatArr>     FloatArrArr;
typedef Array<IntArrArr>    IntArrArrArr;
typedef Array<FloatArrArr>  FloatArrArrArr;
typedef Array<IntArrArrArr> IntArrArrArrArr;

typedef Vec<float>    	    FloatVec;
typedef Vec1<float>         V1f;
typedef Vec2<int>           V2i;
typedef Vec2<float>         V2f;
typedef Vec2<double>        V2d;
typedef Vec3<int>           V3i;
typedef Vec3<float>         V3f;
typedef Vec3<double>        V3d;
typedef Vec4<int>           V4i;
typedef Vec4<float>         V4f;
typedef Vec4<double>        V4d;
typedef Vec5<int>           V5i;
typedef Vec5<float>         V5f;
typedef Vec5<double>        V5d;

typedef Array<FloatVec>     FloatVecArr;
typedef Array<V1f>          V1fArr;
typedef Array<V2i>    	    V2iArr;
typedef Array<V2f>    	    V2fArr;
typedef Array<V2d>    	    V2dArr;
typedef Array<V3i>    	    V3iArr;
typedef Array<V3f>    	    V3fArr;
typedef Array<V3d>    	    V3dArr;
typedef Array<V4i>    	    V4iArr;
typedef Array<V4f>    	    V4fArr;
typedef Array<V4d>    	    V4dArr;
typedef Array<V2iArr>	    V2iArrArr;
typedef Array<V2fArr>  	    V2fArrArr;
typedef Array<V3iArr>       V3iArrArr;
typedef Array<V3fArr>  	    V3fArrArr;
typedef Array<V4iArr>       V4iArrArr;
typedef Array<V4fArr>       V4fArrArr;
typedef Array<V2iArrArr>    V2iArrArrArr;
typedef Array<V2fArrArr>    V2fArrArrArr;
typedef Array<V3fArrArr>    V3fArrArrArr;
typedef Array<V4iArrArr>    V4iArrArrArr;

//#ifndef WIN32
#include "Vector.cpp"
//#endif

#endif