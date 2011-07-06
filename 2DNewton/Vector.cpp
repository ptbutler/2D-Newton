/*
  File: 	 Vector.cpp
  Author: 	 J.K. Johnstone 
  Created:	 11 March 1998
  Last Modified: 7 July 2009
  Purpose:       A class of arrays and vectors.
  History:	 

     5/9/00:   Added operator!=.
     8/24/00:  Added merge.
     10/26/00: Changed const variables to normal variables 
               in init, constructors, create, operator= and dot.
	       Added create(int) and create (int, Type).
     12/5/00:  Added Vec3:unitNorm.
     12/17/00: Added Vec:operator+, Vec:operator-, Vec:operator*.
	       Motivation (but Blinn is wrong: + is binary, not ternary!): 
	       Jim Blinn's Corner, IEEE CG+A, July/Aug 2000, p. 97
     12/18/00: Moved Vec:dot from private to public.
     12/27/00: Added Vec1 (for use by BezierSurf1).
     1/28/01:  Changed Array:allocate (added if (N>0) check before allocation).
     7/26/01:  Added Vec:midPt and Array:clear.
     8/2/01:   Added Vec2:create and Vec2:circleCenter.
     8/10/01:  Moved volParallelepiped from friend of Vec3 to MiscVector.h.
     8/13/01:  Added Array:append.
     4/16/02:  Added Vec2:leftTurn.
     7/30/02:  Added Array:create(Array &a, int N).
     7/31/02:  Added fabs to distFromPlane to change it from signed distance to distance.
     8/14/02:  Added Array:fabs and Array:operator-.
               Changed countApproxDuplicate.
     8/21/02:  Added deleteApproxEq.
     8/28/02:  Added Array:+=.
     10/9/02:  Added Array:set.
     11/27/02: Added Vec2:intersectInterval.
     12/19/02: Changed default last parameter to 0 in Vec3:distFromPlane.
     3/4/03:   Added Array:shrink.
     5/8/03:   Added reverseIt.
     6/20/03:  Moved VecMatrix to Matrix.h.
               Moved Vec2:intersectInterval (6 parameter version) to LineTri3.cpp,
	       which allows removal of Line.h library.
     6/22/03:  Moved Vec2:circleCenter to LineTri.cpp (that now contains HomogLine.cpp),
               which allows removal of HomogLine.h library.
	       Moved Vec:principalAxes and drawPrincipalAxes to Matrix,
	       which allows removal of Matrix.h library.
	       Added Array:median, and its helper functions select and partition.
	       Added Array:mergeSort, and helper function merge (version with sortIndex).
	       Added Array:binarySearch.
	       Added Vec:operator+ and Vec:operator-.
     6/23/03:  Moved removeDuplicate to protected (it is only a helper function).
               Removed Array:print(): it is subsumed by operator<<.
	       Changed Array:mean to return mean rather than passing back through parameter.
	       Removed Vec:scalarmult, which is superseded by Vec:operator*=.
	       Removed Vec:swap(int i, int j), superseded by Array:swap(int i, int j).
	       Removed Vec2:clear(), which is superseded by Array:clear().
     6/28/03:  Added insertionSort.
     7/25/03:  Added Vec:swap(int i, int j).
     9/17/03:  Added ExpandableArray.
     9/19/03:  Added angleUnit for quicker computation for unit vectors (typically normals).
     11/10/03: Added 2nd Array:append, of just one element.
     7/29/04:  Moved VecMatrix here from Matrix.  It is logically a 3d array.
     9/19/04:  Added nSubset and subsetSize.
     9/21/04:  Added 4-dimensional cross product (Vec4:cross), necessary for
               quaternion rotation in 4-space.
	       Actually, a generalized cross product in n-space is defined.
	       Added 4-dimensional rotation matrix routine (Vec4:rotToCoord).
	       Added Vec:matrixMult (general matrix-vector multiplication).
     5/12/05:  Removed Vec1, Vec5, Vec2Par, Vec3Par, Vec4Par, 
               ExpandableArray, VecMatrix (see archived copies in ~/Cbin/old).
     10/19/05: Added Vector:angle().
     11/4/05:  Added Vec2:contains.
     2/23/06:  Added Array::mean(Type &m)
     10/1/07:  Added Vec3:xrot/yrot/zrot.
     11/23/07: Added Vec3:rotToCoord.
     4/23/08:  Added Vec2:angleFromCenter.
     5/6/08:   Added Array:insert (Type, int) and Array:rotateToFront.
     5/30/08:  Added allocateJustLike.
     12/19/08: Added Union.
     6/11/09:  Moved Vec4 and Vec5 to Cbin/quaternion: overkill here
               (also allows Vector to be independent of Miscellany,
	       since only use was of subdet in rotToCoord).
               Moved VecMatrix to BezierCurve.
	       Removed Array2.
     late 6/09: polishing
     7/6/09:   removed init (instead embedded in constructors)
               replaced allocate(N) by equivalent create(N) 
	       (any calls to allocate in outside software can be replaced by calls to create)
	       changed name of 'reallocate' to 'grow'
     7/22/09:  Need V4f in Vector, so moved Vec4 and Vec5 back from Cbin/quaternion,
               but without Vec4::cross and Vec4::rotToCoord (since they require Miscellany),
	       which are left in quaternion directory
     10/2/09:  Moved det/subdet here from Miscellany, 
	       as helper functions to rotToCoord.
	       Replaced call of deg2rad to explicit code, again
	       to avoid dependence on Miscellany.
*/




#include <assert.h>
#include <iostream>
using namespace std;
#include <math.h>

//#include "basic/Vector.h"
#include "Vector.h"


/****************************************************************************
->os: output stream
->a: array
<-: output stream avec array
****************************************************************************/

template <class Type> ostream& 
operator<< (ostream &os, Array<Type> &a)
{
  if (a.n == 0) { os << "()"; return os; }
  os << "(";
  for (int i=0; i<a.n-1; i++)
    os << a.v[i] << ",";
  os << a.v[a.n-1] << ")";
  return os;
}

/****************************************************************************
  Absolute value of maximal element of an array.

  A strange beast, required in countApproxDuplicate.

->a: input array
<-: absolute value of the array's maximum
****************************************************************************/

template <class Type> Type 
fabs(Array<Type> a)             
{
  Type max = fabs(a[0]);  // this allows fabs(float) or fabs(array) to be called
  for (int i=1; i<a.n; i++)
    if (fabs(a[i]) > max)
        max = fabs(a[i]);
  return max;
}

/****************************************************************************
  Constructor: initialization by built-in vector type.

->a: input array
->N: # of elements to copy over
****************************************************************************/

template <class Type>
Array<Type>::Array (Type *a, int N)
{
  v = new Type[n=N]; 
  for (int i=0; i<n; i++)  v[i] = a[i];
}

/****************************************************************************
  Constructor: initialization by another array object.

->a: input array
****************************************************************************/

template <class Type>
Array<Type>::Array (Array<Type> &a)
{
  v = new Type[n=a.n]; 
  for (int i=0; i<n; i++)  v[i] = a.v[i];
}

/****************************************************************************
  Set array to first N elements of another.
  Unlike constructor, 'create' may be used when object already contains data.

->a: input array
->N: # of elements to copy over
****************************************************************************/

template <class Type> void
Array<Type>::create (Array &a, int N)
{
  if (n>0) delete [] v;
  v = new Type[n=N];
  for (int i=0; i<n; i++) v[i] = a[i];
}

/****************************************************************************
  Allocate this multidimensional array to be the same size (in each dimension)
  as the given multidimensional array.
  The two arrays must have exactly the same number of dimensions, and the same type!  
  For example, this will not work for allocating a FloatArrArrArr 
  to be of the same size as a V3fArrArrArr.

->a: the array to emulate (in memory allocation)
->nLevel: # of levels of allocation
****************************************************************************/

template <class Type> void
Array<Type>::allocateJustLike (Array<Type> &a, int nLevel)
{
  assert (nLevel >= 1);
  create (a.getn());
  if (nLevel > 1)
    for (int i=0; i<a.getn(); i++)
      v[i].allocateJustLike (a[i], nLevel-1);
}

/****************************************************************************
  Add n elements to this array.
  Author: Ross Ptacek.

->n: # of elements to add to the end of the array
****************************************************************************/

template <class Type> void
Array<Type>::grow (int n)
{
  Array<Type> old(*this);
  this->create(this->getn()+n);
  
  for (int i=0; i<old.getn(); i++)
      (*this)[i] = old[i];
}

/****************************************************************************
  Shrink an array (that was declared a maximal size) to its true size.

->truesize: the desired size
****************************************************************************/

template <class Type> void
Array<Type>::shrink (int truesize)
{
  // if (n==0) return;
  int i;
  Type *temp;
  temp = new Type[n=truesize];
  for (i=0; i<truesize; i++) temp[i] = v[i];
  delete [] v;
  v = new Type[n=truesize]; 
  for (i=0; i<truesize; i++) v[i] = temp[i];
  delete [] temp;
}

/****************************************************************************
  Assignment.

->a:
****************************************************************************/

template <class Type> Array<Type>&
Array<Type>::operator= (Array<Type> &a)
{
  if (this == &a) return *this;  
  if (n>0) delete [] v;
  v = new Type[n=a.n]; 
  for (int i=0; i<n; i++)  v[i] = a.v[i];
  return *this;
}

/****************************************************************************
  Equality test.

->a:
<-: equal?
****************************************************************************/

template <class Type> int
Array<Type>::operator== (Array<Type> &a)
{
  if (n != a.n) return(0);
  for (int i=0; i<n; i++) if (v[i] != a[i]) return(0);
  return(1);
}

/****************************************************************************
  Inequality test.

->a:
<-: unequal?
****************************************************************************/

template <class Type> int
Array<Type>::operator!= (Array<Type> &a)
{
  return (!(*this == a));
}

/****************************************************************************
  Element-wise subtraction.

->b: array to be subtracted from this array
<-: difference of these two arrays
****************************************************************************/

template <class Type> Array<Type>
Array<Type>::operator- (Array<Type> &b)
{
  Array<Type> tmp;
  tmp.create(n);
  for (int i=0; i<n; i++)  tmp[i] = v[i] - b[i];
  return tmp;
}

/****************************************************************************
  Splice an array to the end of the present array.

->a: array to be spliced on
****************************************************************************/

template <class Type> void
Array<Type>::operator+= (Array<Type> &a)
{
  int i;
  Array<Type> tmp(n+a.n);  int nTot=0;
  for (i=0; i<n;   i++) tmp[nTot++] = v[i];
  for (i=0; i<a.n; i++) tmp[nTot++] = a.v[i];
  create(tmp);
}

/****************************************************************************
  Set the array to the concatenation of two arrays.

->a,b: the two arrays
****************************************************************************/

template <class Type> void
Array<Type>::append (Array<Type> &a, Array<Type> &b)
{
  int i;
  create(a.n + b.n);
  for (i=0; i<a.n; i++)		v[i] = a[i];
  for (i=0; i<b.n; i++)		v[a.n + i] = b[i];
}

/****************************************************************************
  Append an element to the end of this array.

->a: element to be appended
****************************************************************************/

template <class Type> void
Array<Type>::append (Type &a)
{
  Array<Type> old (*this);
  create (n + 1);
  for (int i=0; i<old.n; i++)	v[i] = old[i];
  v[n-1] = a;
}

/****************************************************************************
  Binary search for x in the (v[i],v[j]) part of this sorted array.
  If x is not in the array, return -1, otherwise return index of x.

        NOT TESTED.

->i,j: interval of array to search
->x: element to be searched for
<-: index of x (or -1 if not present)
****************************************************************************/

template <class Type> int
Array<Type>::binarySearch (int i, int j, Type x)
{
  if (j-i < 4)  // search exhaustively when array is very small
    {
      int index = i; 
      while (index <= j && v[index] != x) index++;
      if (index <= j) return index;
      else return -1;
    }     
  int half = (int) floor ((j-i)/2);
  if      (x ==v[i+half]) return (i+half);
  else if (x < v[i+half]) return binarySearch (i,i+half-1,x);
  else                    return binarySearch (i+half+1,j,x);
}

/****************************************************************************
  Bubblesort (in increasing order).
  Keep track of correspondence between original and new positions:
  new position i = old position sortIndex[i] 

this: array of scalars
<-sortIndex:
****************************************************************************/

template <class Type> void
Array<Type>::bubbleSort (Array<int> &sortIndex)
{
  int i,j;
  sortIndex.create (n);  for (i=0; i<n; i++) sortIndex[i] = i;
  for (j=n-1; j>=1; j--)	// bubble largest elt to position j
    for (i=0; i<j; i++)
      if (v[i+1] < v[i])
       {
        swap (i,i+1);
	sortIndex.swap (i,i+1);
       }
}

/****************************************************************************
  Bubblesort in decreasing order.
  Keep track of correspondence between original and new positions:
  new position i = old position sortIndex[i] 

this: array of scalars
<-sortIndex:
****************************************************************************/

template <class Type> void
Array<Type>::bubbleSortReverse (Array<int> &sortIndex)
{
  int i,j;
  sortIndex.create (n);  for (i=0; i<n; i++) sortIndex[i] = i;
  for (j=n-1; j>=1; j--)	// bubble smallest elt to position j
    for (i=0; i<j; i++)
      if (v[i+1] > v[i])
       {
        swap (i,i+1);
	sortIndex.swap (i,i+1);
       }
}

/****************************************************************************
  Find the duplicates in this *sorted* array,
  in preparation for removing these duplicates.

<-duplicate: marks the duplicates: duplicate[i] = 1 iff v[i] is a duplicate.
<-: # of duplicates
****************************************************************************/

template <class Type> int
Array<Type>::countDuplicate (Array<int> &duplicate)
{
  if (n==0) { duplicate.create(0); return(0); }
  int nDuplicate=0;  duplicate.create(n); duplicate[0]=0;
  for (int i=1; i<n; i++)
    if (v[i] == v[i-1])			// duplicate intersection
       { duplicate[i] = 1; nDuplicate++; }
    else duplicate[i] = 0;
  return (nDuplicate);
}

/****************************************************************************
  Remove the duplicates in this *sorted* array.
****************************************************************************/

template <class Type> void
Array<Type>::deleteDuplicate ()
{
  IntArr duplicate;
  int nDuplicate = countDuplicate (duplicate);
  removeDuplicate (nDuplicate, duplicate);
}

/****************************************************************************
  Find the approximate duplicates in this array (no sorting assumed),
  in preparation for removing these duplicates.

  Only applicable to FloatArr, V2fArr, V3fArr, ....  
  and CommonTangentArr (as of 2/11/04).

<-duplicate: marks the duplicates: duplicate[i] = 1 iff v[i] is a duplicate.
->eps: two entries are considered duplicates if they differ by less than this amt
<-: # of duplicates
****************************************************************************/

template <class Type> int
Array<Type>::countApproxDuplicate (Array<int> &duplicate, float eps)
{
  int i,j;
  if (n==0) { duplicate.create(0); return(0); }
  int nDuplicate=0;  duplicate.create(n); 
  for (i=0; i<n; i++) duplicate[i]=0;
  for (i=0; i<n; i++)
    if (!duplicate[i])
      for (j=i+1; j<n; j++)
       {
	Type diff(v[j]);
	diff -= v[i];
        if (!duplicate[j] && fabs(diff) < eps)
         {
          duplicate[j] = 1;
          nDuplicate++;
         }
       }
  return (nDuplicate);
}

/**********************************************************************  
  // old version for sorted data
  int last=0;
  for (i=1; i<n; i++)
   {
    Array<Type> diff(v[i]); diff -= v[last];
    if (fabs(diff) < eps)	// duplicate intersection
         { duplicate[i] = 1; nDuplicate++; }
    else { duplicate[i] = 0; last = i; }
   }
***********************************************************************/    

/****************************************************************************
  Remove the approximate duplicates in this array (no sorting assumed).
  Only applicable to FloatArr, V2fArr, V3fArr, ....  
  and CommonTangentArr (as of 2/11/04).

->eps: two entries are considered duplicates if they differ by less than this amt
****************************************************************************/

template <class Type> void
Array<Type>::deleteApproxDuplicate (float eps)
{
  IntArr duplicate;
  int nDuplicate = countApproxDuplicate (duplicate, eps);
  removeDuplicate (nDuplicate, duplicate);
}

/****************************************************************************
  Find elements that are approximately equal (within EPS) to VAL.
  Only applicable to FloatArr.

->val: value of interest
<-expendable: array marking the entries approximately equal to 'val'
->eps: amount of allowed deviation from 'val'
<-: # of near-duplicates of val
****************************************************************************/

template <class Type> int
Array<Type>::countApproxEq (float val, Array<int> &expendable, float eps)
{
  if (n==0) { expendable.create(0); return(0); }
  int nExpendable=0;  expendable.create(n); 
  for (int i=0; i<n; i++) 
    if (fabs (v[i] - val) < eps) 	
     {   expendable[i]=1;  nExpendable++; }
    else expendable[i]=0;
  return (nExpendable);
}

/****************************************************************************
  Delete near-duplicates of a value of interest (entries within EPS of VAL).
  Only applicable to FloatArr.

->val: value of interest
->eps: amount of allowed deviation from 'val'
****************************************************************************/

template <class Type> void
Array<Type>::deleteApproxEq (float val, float eps)
{
  IntArr duplicate;
  int nDuplicate = countApproxEq (val, duplicate, eps);
  removeDuplicate (nDuplicate, duplicate);
}

/****************************************************************************
  Filter elements of w out of this array,
  returning the filtered elements in filterArr.

  Restrictions:
  - this array must be FloatArr (since we use fabs).
  - w must be sorted

  Equivalent to computing *this - w (set difference).
  Walk through array, only gathering value if it is not in w.
 
->w: a sorted array of elements to filter out from our array
<-filterArr: the filtered array
->eps: elements within eps are considered equivalent
****************************************************************************/

template <class Type> void
Array<Type>::filterOut (Array<Type> &w, Array<Type> &filterArr, float eps)
{
  int i;
  int nKept=0, wi=0;
  for (i=0; i<n; i++)
   {
    while (wi<w.getn() && w[wi]+eps < v[i])         wi++;
    if    (wi<w.getn() && fabs(w[wi] - v[i]) > eps) nKept++; 
   }
  if (nKept == 0) { filterArr.create(0); return; }
  filterArr.create(nKept);  
  nKept=wi=0;
  for (i=0; i<n; i++)
   {
    while (wi<w.getn() && w[wi]+eps < v[i])	    wi++;
    if    (wi<w.getn() && fabs(w[wi] - v[i]) > eps) filterArr[nKept++] = v[i];
    else if (wi<w.getn()) 
      cout << "Filtering out " << v[i] << " because of " << w[wi] << endl;
   }
}

/****************************************************************************
  In a sorted array, find the consecutive pair [v[i],v[i+1]) that contains x.

this: sorted array
->x:
<-: index of first element of consecutive pair that contains x
****************************************************************************/

template <class Type> int
Array<Type>::findInterval (Type x)
{
  assert (x >= v[0] && x <= v[n-1]);
  if (x==v[n-1]) return (n-1);
  int i=0;
  while (x >= v[i+1]) 	i++;
  return(i);
}

/****************************************************************************
  Insert an element into this array immediately *after* the given index.

->x: element to be inserted
->index: index of the element that should precede the inserted element
****************************************************************************/

template <class Type> void
Array<Type>::insert (Type x, int index)
{
  int i,j;
  Array<Type> tmp(getn() + 1);
  for (i=0,j=0; i<=index;   i++) tmp[j++] = v[i];
  tmp[j++] = x;
  for (i=index+1; i<getn(); i++) tmp[j++] = v[i];
  create (tmp);
}

/****************************************************************************
  Insert an array into this array immediately *after* the given index.

->a: array to be inserted
->index: index of the element that should precede the inserted array
****************************************************************************/

template <class Type> void
Array<Type>::insert (Array<Type> &a, int index)
{
  int i,j;
  Array<Type> tmp(getn() + a.getn());
  for (i=0,j=0; i<=index;   i++) tmp[j++] = v[i];
  for (i=0;     i<a.getn(); i++) tmp[j++] = a[i];
  for (i=index+1; i<getn(); i++) tmp[j++] = v[i];
  create (tmp);
}

/****************************************************************************
  Insertion sort (in increasing order).
  Good for nearly sorted arrays.

        NOT TESTED.									     

this: array of scalars
<-sortIndex: correspondence between original and new positions;
             new position i = old position sortIndex[i]
****************************************************************************/

template <class Type> void
Array<Type>::insertionSort (Array<int> &sortIndex)
{
  int i,j;
  sortIndex.create (n);  for (i=0; i<n; i++) sortIndex[i] = i;
  for (i=1; i<n; i++)
    {  // insert v[i] into its correct position
      Type tmp = v[i];  j=i-1;
      while (j>=0 && tmp < v[j])
	{
	  v[j+1] = v[j]; // push tmp past it
	  sortIndex.swap (j,j+1); // added 5/27/04
	  j--;
	}
      v[j+1] = tmp;
    }
}

/****************************************************************************
  Mean of this numerical array.

<-: mean
****************************************************************************/

template <class Type> Type
Array<Type>::mean ()
{
  assert (n>0);
  Type average = v[0];
  for (int i=1; i<n; i++)  average += v[i];
  average /= n;
  return (average);
}

/****************************************************************************
  Mean of this numerical multidimensional array.

<-m: mean
****************************************************************************/

template <class Type> void
Array<Type>::mean (Type &m)
{
  assert (n>0);
  m = v[0];
  for (int i=1; i<n; i++)  m += v[i];
  m /= n;
}

/****************************************************************************
  Compute median.

	NOT TESTED YET.

<-: median
****************************************************************************/

template <class Type> Type
Array<Type>::median()
{
  if (n%2 == 1) return select (0, n-1, (n-1)/2);
  else          return select (0, n-1, n/2);
}

/****************************************************************************
  Merge two sorted arrays a and b into this sorted array.

->a,b: two sorted arrays
****************************************************************************/

template <class Type> void
Array<Type>::merge (Array<Type> &a, Array<Type> &b)
{
  create(a.n + b.n);
  int pa=0;		// ptr to a
  int pb=0;		// ptr to b
  for (int i=0; i<n; i++)
    if (pa<a.n && (pb == b.n || a[pa] < b[pb])) v[i] = a[pa++]; 
    else 					v[i] = b[pb++];
}

/****************************************************************************
  Merge two sorted arrays a and b into this sorted array,
  while maintaining sort indexes.
  a is the sorted version of the first half of an array
  and b is the sorted version of the second half.
  aindex gives the sort permutation of a, and
  bindex gives the sort permutation of b.
  sortIndex gives the sort permutation of the new sorted array wrt original array.

        NOT TESTED.									   

->a,b: two sorted arrays (during mergesort)
->aindex,bindex: sort indices of these arrays
<-sortIndex: sort index of new merged array
****************************************************************************/

template <class Type> void
Array<Type>::merge (Array<Type> &a, IntArr &aindex,
		    Array<Type> &b, IntArr &bindex,
		    IntArr &sortIndex)
{
  create(a.n + b.n);
  int pa=0;		// ptr to a
  int pb=0;		// ptr to b
  for (int i=0; i<n; i++)
    if (pa<a.n && (pb == b.n || a[pa] < b[pb])) 
         { v[i] = a[pa]; sortIndex[i] = aindex[pa++]; }
    else { v[i] = b[pb]; sortIndex[i] = bindex[pb++] + a.n; }
}

/****************************************************************************
  Merge the sorted array a with this sorted array.

->a: another sorted array
this: *sorted* array
****************************************************************************/

template <class Type> void
Array<Type>::merge (Array<Type> &a)
{
  Array<Type> tmp;	tmp.merge (*this,a);
  create (tmp);
}

/****************************************************************************
  Mergesort (in increasing order).
  Keep track of correspondence between original and new positions:
  new position i = old position sortIndex[i] 

        NOT TESTED.

this: array of scalars
<-sortIndex:
****************************************************************************/

template <class Type> void
Array<Type>::mergeSort (IntArr &sortIndex)
{
  int i,j;
  if (n<10) { bubbleSort (sortIndex); return; }
  sortIndex.create(n);
  for (i=0; i<n; i++) sortIndex[i] = i;
  int n2 = (int) floor(n/2);
  Array<Type> firstHalf; firstHalf.create (*this, n2);
  Array<Type> secondHalf(n-n2); for (i=0,j=n2; j<n; i++,j++) secondHalf[i] = v[j];
  IntArr firstHalfSortIndex, secondHalfSortIndex;
  firstHalf.mergeSort  (firstHalfSortIndex);
  secondHalf.mergeSort (secondHalfSortIndex);
  merge (firstHalf, firstHalfSortIndex, 
	 secondHalf, secondHalfSortIndex, 
	 sortIndex);
}

/****************************************************************************
  Compute the number of subsets of size s in this array,
  with or without overlap.
  Of course, the last subset may not be full.

  The nonoverlap case is obvious.
  In the overlap case, the first subset introduces s elts and the remaining subsets
  introduce s-1 elts.  Thus, after i subsets, we have considered
  i*(s-1) + 1 elements.  If it splits perfectly, n = i*(s-1) + 1, so 
  i = (n-1)/(s-1), but of course there may be some leftover, so we need
  ceil (n-1)/(s-1) subsets.

->s:       # of points per subset (except possibly the last)
->overlap: 1 iff the first element of i+1st subset should
           duplicate the last element of the ith subset
****************************************************************************/

template <class Type> int
Array<Type>::nSubset (int s, int overlap)
{ 
  return (int (overlap ? ceil(float(n-1) / float(s-1))
	               : ceil(float(n)   / float(s))));
}

/****************************************************************************
  Having marked the duplicates in this array (using countDuplicate), 
  remove them.

->nDuplicate: # of duplicates
->duplicate: marks the duplicates: duplicate[i] = 1 iff v[i] is a duplicate
****************************************************************************/

template <class Type> void
Array<Type>::removeDuplicate (int nDuplicate, Array<int> &duplicate)
{
  if (nDuplicate==0) return;
  int newSize = n - nDuplicate;
  Array<Type> arrUnique(newSize);
  for (int i=0, iUnique=0; i<n; i++)	// copy over to smaller array
    if (duplicate[i]==0)
      arrUnique[iUnique++] = v[i];
  create (arrUnique);
}

/****************************************************************************
  Reverse this array.
****************************************************************************/

template <class Type> void
Array<Type>::reverseIt()
{
  for (int i=0; i<floor(n/2); i++) 
    swap (i, n-1-i);
}

/****************************************************************************
  Rotate the array so that the specified index becomes the 1st element.

->index: index of element to be rotated to the front of the array
****************************************************************************/

template <class Type> void
Array<Type>::rotateToFront (int index)
{
  assert (0 <= index && index < n);
  int i,j;
  Array<Type> tmp(n);
  for (i=0, j=index; i<n; i++,j=(j+1)%n)
    tmp[i] = v[j];
  create (tmp);
}

/******************************************************************************
  Partition a range v[i]->v[j] of the array about a pivot.
  Used in quicksort: see Fig 8.13, p. 264 of AHU.
  Partitions v[i],...,v[j] so keys < pivot are at the left
  and keys >= pivot are at the right.

	NOT TESTED YET.

->i,j: range of array to be considered
->pivot:
<-: beginning of group on right
******************************************************************************/

template <class Type> int
Array<Type>::partition (int i, int j, Type pivot)
{
  int l=i,r=j;
  do {
    swap (l, r);
    while (v[l] < pivot) l++;
    while (v[r] >= pivot) r--;
  } while (l <= r);
  return (l);
}

/******************************************************************************
  Select the kth largest element from the subrange v[i],...,v[j].
  See AHU, Data Structures and Algs, p. 286.

	NOT TESTED YET.

->i,j: subrange
->k: 
<-: kth largest element from the subrange v[i],...,v[j]
******************************************************************************/

template <class Type> Type
Array<Type>::select (int i, int j, int k)
{
  int index = i+1;    // base case if v[i],...,v[j] are all equal
  while (index <= j && v[index] == v[i]) index++;
  if (index > j)  return (v[i]);
    
  Type pivot = v[i];
  int m = partition (i,j,pivot);
  if (k <= m-i) return (select (i, m-1, k));
  else          return (select (m, j, k-m+i));
}

/****************************************************************************
  Compute the size of each subset when this array is split into
  subsets of size s, with or without overlap.
  Trivial except for the last two subsets.
  The last two subsets are split evenly, to avoid a tiny last subset.

->s:          # of points per subset (except possibly the last)
->overlap:    1 iff the first element of i+1st subset should
              duplicate the last element of the ith subset
<-subsetSize: size of each subset
****************************************************************************/

template <class Type> void 
Array<Type>::subsetSize (int s, int overlap, IntArr &subsetsize)
{
  int nSub = nSubset(s, overlap);
  subsetsize.create (nSub);
  if (nSub == 1) 
    { subsetsize[0] = n; return; }
  for (int i=0; i<nSub-2; i++)
    subsetsize[i] = s;
  int nLastTwoSubset;
  if (overlap)
      // have now handled (s-1)*(nSub-2) + 1 elements of the set,
      // since first set has s elts and other subsets have s-1 elts
       nLastTwoSubset= n - ((s-1)*(nSub-2)); // # in last 2 subsets
  else nLastTwoSubset= n - s*(nSub-2);
  if (overlap)
   {
    subsetsize[nSub-2] = int (ceil (float(nLastTwoSubset+1) / 2.));
    subsetsize[nSub-1] = int (floor(float(nLastTwoSubset+1) / 2.));
   }
  else
   {
    subsetsize[nSub-2] = int (ceil (float(nLastTwoSubset) / 2.));
    subsetsize[nSub-1] = int (floor(float(nLastTwoSubset) / 2.));
   }
}

/****************************************************************************
  Swap two entries of the array.

->i,j: index of two entries to be swapped
****************************************************************************/

template <class Type> void
Array<Type>::swap (int i, int j)
{
  Type tmp;
  tmp = v[i]; v[i] = v[j];  v[j] = tmp;
}

/********************************************************************************
  Union of this array (interpreted as a set)
  with the set consisting of just one element x.
  The set is stored in the first 'size' elements of the array (although the
  array is larger).
  If x is not already in the set, it is added and the size increases by 1.

->x: element to be added (if not already present)
<->size: size of the set in the array (which will be incremented if x is added)
********************************************************************************/

template <class Type> void
Array<Type>::Union (Type x, int &size)
{
  int present=0;
  for (int i=0; i<size; i++)
    if (x == v[i]) 
      present=1;
  if (!present)
    {
      if (size >= n) 
	{ cout << "Union: array is too small to add an element"<<endl; exit(-1); }
      v[size] = x;
      size++;
    }
}

/****************************************************************************
  Min of array and its index.

<-index: index of minimum element
<-: minimum element
****************************************************************************/

template <class Type> Type
Array<Type>::vecmin (int &index)
{
  index = 0;
  for (int i=1; i<n; i++)
    if (v[i] < v[index])  
      index = i;
  return (v[index]);
}

/****************************************************************************
  Max of array and its index.

<-index: index of max element
<-: max element
****************************************************************************/

template <class Type> Type
Array<Type>::vecmax (int &index)
{
  index = 0;
  for (int i=1; i<n; i++)
    if (v[i] > v[index])  
      index = i;
  return (v[index]);
}

/*********************************************************************************************
 ********************************************************************************************
 ********************************************************************************************
 ********************************************************************************************/

/****************************************************************************
this: unit vector!
<-: angle with positive x-axis, in radians, \in [0,2*PI)
****************************************************************************/

template <class Type> float
Vec<Type>::angle()
{
  if (this->v[1]>=0) 
       return acos(this->v[0]); 
  else return 2*M_PI - acos(this->v[0]);
}

/****************************************************************************
  Angle between two vectors.

  cos(angle between A and B) = A.B / \|A\| \|B\|

->B: other vector
<-: angle in radians \in [0,PI]
****************************************************************************/

template <class Type> float
Vec<Type>::angle (Vec &B)
{
  float a = this->dot(B) / (this->length() * B.length());
  if      (a>1)	 return(0);	// shouldn't be, but there is roundoff error
  else if (a<-1) return(M_PI);
  else 		 return(acos(a));
}

/****************************************************************************
this: unit vector
->B: another unit vector
<-: angle between these two vectors, in radians, \in [0,PI]
****************************************************************************/

template <class Type> float
Vec<Type>::angleUnit (Vec &B)
{
  float a = this->dot(B);
  if      (a>1)	 return(0);	// shouldn't be, but there is roundoff error
  else if (a<-1) return(M_PI);
  else 		 return(acos(a));
}

/****************************************************************************
->P,C: two points
<-: angle of difference vector P-C from positive x-axis, \in [0,PI]
****************************************************************************/

template <class Type> float
Vec<Type>::angle (Vec &P, Vec &C)
{
  V2f V(P);  V -= C;  // P-C
  V.normalize();
  return V.angle();
}

/****************************************************************************
  Does every component differ by at most eps?

->a: vector for comparison
->eps: allowed deviation per dimension
<-: 1 iff a agrees within epsilon
****************************************************************************/

template <class Type> int
Vec<Type>::approxeq (Vec<Type> &a, float eps)
{
  if (this->n != a.n) return(0);
  for (int i=0; i<this->n; i++)
   {
    Type foo = this->v[i];  foo -= a[i];
    if ((foo > 0 && foo>eps) || // don't want to use fabs, committing to Type=float
        (foo < 0 && foo<-eps))	
      return(0);		// this component not within eps
   }
  return(1);
}

/****************************************************************************
  Distance between two points.

->w: 2nd point
<-: distance
****************************************************************************/

template <class Type> float
Vec<Type>::dist (Vec &w)
{
  assert(this->n == w.n);
  Vec<Type> x(this->n);
  for (int i=0; i<this->n; i++) x[i] = w[i] - this->v[i];
  return (x.length());
}

/****************************************************************************
  Dot product.

->w: second vector
<-: inner product
****************************************************************************/

template <class Type> float
Vec<Type>::dot (Vec<Type> &w)
{
  assert(this->n == w.n);
  float dot=0.;
  for (int i=0; i<this->n; i++)
    dot += this->v[i] * w.v[i];
  return(dot);
}

/****************************************************************************
  Length of this vector.

<-: length
****************************************************************************/

template <class Type> float
Vec<Type>::length()
{ 
  float len=0.;
  for (int i=0; i<this->n; i++) len += this->v[i]*this->v[i];
  return(sqrt(len));
}

/****************************************************************************
  Matrix-vector multiplication.

->matrix:  matrix A of same dimension as this vector, stored in row-major C++-style
<-product: Ax, where x is this vector
****************************************************************************/

template <class Type> void
Vec<Type>::matrixMult (float **matrix, Vec<Type> &product)
{
  product.clear();
  for (int i=0; i<this->n; i++) 
    for (int j=0; j<this->n; j++)  // ith row . v
      product[i] += matrix[i][j] * this->v[j];
}

/****************************************************************************
  Matrix-vector multiplication, with product stored back in this vector.

->matrix: matrix A of same dimension as this vector, stored in row-major C++-style
this: changed from x to Ax
****************************************************************************/

template <class Type> void
Vec<Type>::matrixMult (float **matrix)
{
  Vec<Type> product(this->n); product.clear();
  for (int i=0; i<this->n; i++) 
    for (int j=0; j<this->n; j++)  // ith row . v
      product[i] += matrix[i][j] * this->v[j];
  *this = product;
}

/****************************************************************************
  Compute midpoint of two points.

->a,b: two points
this: midpoint
****************************************************************************/

template <class Type> void
Vec<Type>::midPt (Vec<Type> &a, Vec<Type> &b)
{ 
  for (int i=0; i<this->n; i++) this->v[i] = (a.v[i] + b.v[i]) / 2.;
}

/****************************************************************************
this: unit vector (in same direction)
****************************************************************************/

template <class Type> void
Vec<Type>::normalize()
{
  float len;
  len = length();
  if (len != 0) for (int i=0; i<this->n; i++) this->v[i] /= len;
}

/****************************************************************************
  Swap two entries of the array.

->i,j: indices of the two entries
****************************************************************************/

template <class Type> void
Vec<Type>::swap (int i, int j)
{
  Type tmp;
  tmp = this->v[i]; this->v[i] = this->v[j];  this->v[j] = tmp;
}

/*********************************************************************************************
 ********************************************************************************************
 ********************************************************************************************
 ********************************************************************************************/

/****************************************************************************
  Angle (of this point) from a center c.

->this: the point, in 2D
->c: the center
<-: the angle ccw from +ve x-axis to p-c, in degrees (between 0 and 360)
****************************************************************************/

template <class Type> float
Vec2<Type>::angleFromCenter (Vec2 &c)
{
  Vec2<Type> v(*this); v -= c;
  return rad2deg (v.angleFromX());
}

/****************************************************************************
  Angle (of this vector) with x-axis.

  cos(angle) = adjacent / hypotenuse   gives angle between [0,PI].
  Flip if v[1] < 0 (below x-axis).

<-: angle of this vector from x-axis, measured counterclockwise, in radians
    \in [0,2*PI]
****************************************************************************/

template <class Type> float
Vec2<Type>::angleFromX()
{
  float foo = this->v[0] / this->length();
  float angle;
  if (foo>1)       angle = 0; // shouldn't happen, but roundoff
  else if (foo<-1) angle = M_PI;
  else             angle = acos (foo);
  if (this->v[1] >= 0) return (angle);
  else 		 return (2*M_PI - angle);
}

/*****************************************************************************************
  Does t lie in this interval?

this: an interval; if v[1] < v[0], the implied interval is the complement of (v[0],v[1])
->t:
<-: does t lie in this interval?
*****************************************************************************************/

template <class Type> int
Vec2<Type>::contains (float t)
{
  if (this->v[0] < this->v[1])
    return (t >= this->v[0] && t <= this->v[1]);
  else
    return (t <= this->v[0] || t >= this->v[1]);
}

/****************************************************************************
  Intersect two intervals.

->b: second interval
<-intersection: intersection of the two intervals
<-: 0 iff no intersection
****************************************************************************/

template <class Type> int
Vec2<Type>::intersectInterval (Vec2 &b, Vec2 &intersection)
{
  float aMin = min (this->v[0], this->v[1]),
        aMax = max (this->v[0], this->v[1]),
	bMin = min (b.v[0], b.v[1]),
	bMax = max (b.v[0], b.v[1]);
  if      (aMax < bMin || bMax < aMin)	return (0);		// disjoint
  if      (aMin < bMin && bMax < aMax) 	intersection = b;	// contained
  else if (bMin < aMin && aMax < bMax)	intersection = *this;
  else if (aMin < bMin)						// overlap
       { intersection[0] = bMin; intersection[1] = aMax; }
  else { intersection[0] = aMin; intersection[1] = bMax; }
  return 1;
}

/****************************************************************************
->P,Q: two points in 2-space
<-: does PQv form a left turn (where v is this vector)?
****************************************************************************/

template <class Type> int
Vec2<Type>::leftTurn (Vec2 &P, Vec2 &Q)
{
  if (P[0]*Q[1] - P[1]*Q[0] + 
      Q[0]*this->v[1] - Q[1]*this->v[0] + 
      this->v[0]*P[1] - this->v[1]*P[0] > 0)
       return (1);
  else return (0); 
}

/****************************************************************************
  Rotate the 2-vector by theta.

->c: cos(theta)
->s: sin(theta)
****************************************************************************/

template <class Type> void
Vec2<Type>::rotate (Type c, Type s)
{
  Type oldv0 = this->v[0];
  this->v[0] = c*this->v[0]  - s*this->v[1];
  this->v[1] = s*oldv0       + c*this->v[1];
}

/*********************************************************************************************
 ********************************************************************************************
 ********************************************************************************************
 ********************************************************************************************/

/****************************************************************************
  Cross product.

->a,b: two vectors
this: cross product of a and b
****************************************************************************/

template <class Type> void
Vec3<Type>::cross (Vec3 &a, Vec3 &b)
{
  this->v[0] = a[1]*b[2] - a[2]*b[1];
  this->v[1] = a[2]*b[0] - a[0]*b[2];
  this->v[2] = a[0]*b[1] - a[1]*b[0];
}

/****************************************************************************
  Distance from plane.

  Distance of Q from the plane is the solution of
  ((Q+tN) - P) . N = 0, which is t = distance = (P-Q) . N

->P: point of plane
->N: normal of plane
->signedDist: use signed distance?
<-: distance of this point from plane (P,N)
****************************************************************************/

template <class Type> float
Vec3<Type>::distFromPlane (Vec3 &P, Vec3 &N, int signedDist)
{
  Vec3<Type> foo(P);
  foo -= (*this);
  if (signedDist) return (foo * N);
  else 	          return (fabs (foo*N));
}

/***********************************************************************************
  Compute the rotation matrix that rotates a unit vector to a coordinate axis.

  Rotation matrix has ith row = v, and other rows (v x ej) and (v x ej) x v,
  where ej is the coordinate axis furthest from v.
  Better to choose coordinate axis furthest from v, for robustness.
  That is, to avoid a near-singular matrix in cross product.

->this:      unit vector
->coordAxis: the coordinate axis (by index: 0=x-axis, 1=y-axis, 2=z-axis)
<-rot:       3x3 rotation matrix, stored in row-major order C++-style, 
             already allocated
***********************************************************************************/

template <class Type> void
Vec3<Type>::rotToCoord (int coordAxis, float **rot)
{
  assert (fabs(this->length() - 1) < .0001);  // must be unit
  int i,j;
  for (i=0; i<3; i++) for (j=0; j<3; j++) rot[i][j] = 0;
  V3fArr e(3); // coordinate axes
  for (i=0; i<3; i++) for (j=0; j<3; j++) e[i][j] = 0;
  for (i=0; i<3; i++) e[i][i] = 1;
  // find furthest coordinate axis (from both v and -v)
  V3f vneg (-this->v[0], -this->v[1], -this->v[2]);
  V3f distance; 
  for (i=0; i<3; i++) 
    distance[i] = min (this->dist(e[i]), vneg.dist(e[i])); // measure dist from v and -v
  IntArr sortIndex; distance.bubbleSort(sortIndex);
  int furthest = sortIndex[2];

  // build rotation matrix
  V3f row2, row3; // other rows (we already know that v is one of them)
  row2.cross (*this, e[furthest]);
  row3.cross (*this, row2);
  row2.normalize(); row3.normalize();
  for (i=0; i<3; i++)
    {
      rot[coordAxis][i] = this->v[i];
      rot[(coordAxis+1)%3][i] = row2[i];
      rot[(coordAxis+2)%3][i] = row3[i];
    }
  if (fabs(det(3,rot) - 1) > .0001) // choose row order to create rot (det 1), 
                                    // not reflection (det -1)
    for (i=0; i<3; i++) // flip 2 rows
      {
	rot[(coordAxis+1)%3][i] = row3[i];
	rot[(coordAxis+2)%3][i] = row2[i];
      }
  assert (fabs(det(3,rot) - 1) < .0001);
}

/****************************************************************************
  Unit normal of plane through 3 points.

->a,b,c: 3 points
this: unit normal of plane defined by a,b,c
****************************************************************************/

template <class Type> void
Vec3<Type>::unitNorm (Vec3 &a, Vec3 &b, Vec3 &c)
{
  Vec3<Type> ba,ca;
  for (int i=0; i<3; i++) { ba[i] = b[i] - a[i]; ca[i] = c[i] - a[i]; }
  cross (ba,ca);
  this->normalize();
}

/****************************************************************************
  Rotate this vector about the x-axis.

  Mirrors OpenGL's glRotatef (nDeg, 1, 0, 0) command.
  Multiplication by the matrix [1 0 0; 0 c -s; 0 s c].

->nDeg: # of degrees to rotate
****************************************************************************/

template <class Type> void
Vec3<Type>::xrot (float nDeg)
{
  float nRad = (nDeg*2*M_PI)/360;
  float c = cos(nRad), s = sin(nRad);
  Vec3<Type> foo(*this);
  this->v[1] = c*foo[1] - s*foo[2];
  this->v[2] = s*foo[1] + c*foo[2];
}

/****************************************************************************
  Rotate this vector about the y-axis.

  Mirrors OpenGL's glRotatef (nDeg, 0, 1, 0) command.
  Multiplication by the matrix [c 0 s; 0 1 0; -s 0 c].

->nDeg: # of degrees to rotate
****************************************************************************/

template <class Type> void
Vec3<Type>::yrot (float nDeg)
{
  float nRad = (nDeg*2*M_PI)/360;
  float c = cos(nRad), s = sin(nRad);
  Vec3<Type> foo(*this);
  this->v[0] = c*foo[0] + s*foo[2];
  this->v[2] = c*foo[2] - s*foo[0];
}

/****************************************************************************
  Rotate this vector about the z-axis.

  Mirrors OpenGL's glRotatef (nDeg, 0, 0, 1) command.
  Multiplication by the matrix [c -s 0; s c 0; 0 0 1].

->nDeg: # of degrees to rotate
****************************************************************************/

template <class Type> void
Vec3<Type>::zrot (float nDeg)
{
  float nRad = (nDeg*2*M_PI)/360;
  float c = cos(nRad), s = sin(nRad);
  Vec3<Type> foo(*this);
  this->v[0] = c*foo[0] - s*foo[1];
  this->v[1] = s*foo[0] + c*foo[1];
}

/**************************************************************************
  Compute determinant of this nxn matrix.

->n: dimension of matrix
->a: matrix
**************************************************************************/

float det (int n, float **a)
{
  if (n==2) return (a[0][0]*a[1][1] - a[1][0]*a[0][1]);  // base case
  float determ=0.;
  for (int i=0; i<n; i++)
    determ += (i%2==0) ?  a[0][i] * subdet (n, a, 0, i)
    		       : -a[0][i] * subdet (n, a, 0, i);
  return (determ);
}

/**************************************************************************
  Compute determinant of minor, resulting from removal
  of given row and column (where first row/col is 0).

->n: dimension of matrix
->a: square nxn matrix
->row, col: row and column to remove from minor
**************************************************************************/

float subdet (int n, float **a, int row, int col)
{
  int i,j;
  float result;
  float **sub;           // create submatrix
  sub=new float*[n-1]; for(i=0; i<n-1; i++) sub[i]=new float[n-1];
  for (i=0; i<n-1; i++)
    for (j=0; j<n-1; j++)
      if (i<row)
       {
        if (j<col) sub[i][j] = a[i][j];  else sub[i][j] = a[i][j+1];
       }
      else
       {
        if (j<col) sub[i][j] = a[i+1][j]; else sub[i][j] = a[i+1][j+1];
       }
  result = det(n-1,sub);

  for(i=0;i<n-1;i++) delete[] sub[i]; delete[] sub;

  return result;
}

