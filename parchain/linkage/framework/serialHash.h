// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef SEIRAL_HASH_INCLUDED
#define SEIRAL_HASH_INCLUDED
#include "utils.h"
#include "sequence.h"
#include "parallel.h"
using namespace std;

#ifndef A_HASH_LINKAGE_PROB
#define A_HASH_LINKAGE_PROB
#endif

#ifndef A_HASH_LINKAGE_PROBE_THRESH
#define A_HASH_LINKAGE_PROBE_THRESH (m)
#endif

template <class HASH, class intT>
class SerialTable {
 private:
  typedef typename HASH::eType eType;
  typedef typename HASH::kType kType;
  intT m;
  intT mask;
  eType empty;
  HASH hashStruct;
  eType* TA;
  intT* compactL;

  // needs to be in separate routine due to Cilk bugs
  static void clearA(eType* A, intT n, eType v) {
    for (intT i=0; i < n; i++) A[i] = v;
  }


  uintT hashToRange(intT h) {return h & mask;}
  intT firstIndex(kType v) {return hashToRange(hashStruct.hash(v));}
  intT incrementIndex(intT h) {return hashToRange(h+1);}
  bool lessIndex(intT a, intT b) {return 2 * hashToRange(a - b) > m;}

 public:
  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  SerialTable(intT size, HASH hashF) :
    m(1 << utils::log2Up(2*size)), 
    mask(m-1),
    empty(hashF.empty()),
    hashStruct(hashF), 
    TA(newA(eType,m)),
    compactL(NULL) 
      { clearA(TA,m,empty); }

  void setActive(intT mm) {
    m = 1 << utils::log2Up(100+(intT)((float)mm));
    mask = m-1;
  }


  void clear() {
    for (intT i=0; i < m; i++) TA[i] = empty;
  }

  // Deletes the allocated arrays
  void del() {
    free(TA); 
    if (compactL != NULL) free(compactL);
  }

  // prioritized linear probing
  //   a new key will bump an existing key up if it has a higher priority
  //   an equal key will replace an old key if replaceQ(new,old) is true
  // returns 0 if not inserted (i.e. equal and replaceQ false) and 1 otherwise
  bool insert2(eType v) {
    kType vkey = hashStruct.getKey(v);
    intT h = firstIndex(vkey);
    while (1) {
      int cmp;
      eType c = TA[h];
      if (c == empty) {
	TA[h] = v;
	return 1;
      } else {
	cmp = hashStruct.cmp(vkey,hashStruct.getKey(c));
	if (cmp == 1) {  // replace and push current value to another bucket
	  TA[h] = v;
	  v = c;
	  vkey = hashStruct.getKey(v);
	} else if (cmp == 0) { // equal (quit or replace)
	  if (hashStruct.replaceQ(v,TA,h)) {
	    // TA[h] = v;
	    return 1;
	  } else return 0; 
	}
      }
      // move to next bucket
      h = incrementIndex(h);
    }
    return 0; // should never get here
  }

 // prioritized linear probing
  //   a new key will bump an existing key up if it has a higher priority
  //   an equal key will replace an old key if replaceQ(new,old) is true
  // returns 0 if not inserted (i.e. equal and replaceQ false) and 1 otherwise
  // // return <success, reached threshold>
  tuple<bool, bool> insert_thresh(eType v) {
    kType vkey = hashStruct.getKey(v);
    intT h = firstIndex(vkey);
    intT prob_ct = 0;
    while (1) {
      int cmp;
      eType c = TA[h];
      if (c == empty) {
	      TA[h] = v;
	      return make_tuple(true, false);
      } else {
        cmp = hashStruct.cmp(vkey,hashStruct.getKey(c));
        if (cmp == 1) {  // replace and push current value to another bucket
          TA[h] = v;
          v = c;
          vkey = hashStruct.getKey(v);
        } else if (cmp == 0) { // equal (quit or replace)
          if (hashStruct.replaceQ(v,TA,h)) {
            // TA[h] = v;
            return make_tuple(true, false);
          } else return make_tuple(false, false);
        }
      }
      prob_ct++;
      if(prob_ct > A_HASH_LINKAGE_PROBE_THRESH){
        return make_tuple(false, true); 
      }
      // move to next bucket
      h = incrementIndex(h);
    }
    return make_tuple(false, true); // should never get here
  }

  bool insert(eType v) {
    kType vkey = hashStruct.getKey(v);
    intT h = firstIndex(vkey);
    while (1) {
      int cmp;
      eType c = TA[h];
      if (c == empty) {
	TA[h] = v;
	return 1;
      } else {
	cmp = hashStruct.cmp(vkey,hashStruct.getKey(c));
	if (cmp == 1) {  // replace and push current value to another bucket
	  TA[h] = v;
	  v = c;
	  vkey = hashStruct.getKey(v);
	} else if (cmp == 0) { // equal (quit or replace)
	  if (hashStruct.replaceQ(v,c)) {
	    TA[h] = v;
	    return 1;
	  } else return 0; 
	}
      }
      // move to next bucket
      h = incrementIndex(h);
    }
    return 0; // should never get here
  }

  bool deleteVal(kType v) {
    intT i = firstIndex(v);
    int cmp;

    // find element to delete
    eType c = TA[i];
    while (c != empty && (cmp = hashStruct.cmp(v,hashStruct.getKey(c))) < 0) {
      i = incrementIndex(i);
      c = TA[i];
    }
    if (cmp != 0) return false; // does not appear
    intT j = i;

    // shift elements after deleted element down to fill hole
    while (c != empty) {
      do {
	j = incrementIndex(j);
	c = TA[j];
      } while (c != empty && lessIndex(i,firstIndex(hashStruct.getKey(c))));
      TA[i] = c;
      i = j;
    }
    return true;
  }

  // Returns the value if an equal value is found in the table
  // otherwise returns the "empty" element.
  // due to prioritization, can quit early if v is greater than cell
  eType find(kType v) {
    intT h = firstIndex(v);
    eType c = TA[h]; 
    while (1) {
      if (c == empty) return empty;
      int cmp = hashStruct.cmp(v,hashStruct.getKey(c));
      if (cmp >= 0) {
	if (cmp == 1) return empty;
	else return c;
      }
      h = incrementIndex(h);
      c = TA[h];
    }
  }

  eType find_thresh(kType v) {
    intT h = firstIndex(v);
    eType c = TA[h]; 
    intT prob_ct = 0;
    while (1) {
      if (c == empty) return empty;
      int cmp = hashStruct.cmp(v,hashStruct.getKey(c));
      if (cmp >= 0) {
	if (cmp == 1) return empty;
	else return c;
      }
      h = incrementIndex(h);
      c = TA[h];
      prob_ct++;
      if(prob_ct > A_HASH_LINKAGE_PROBE_THRESH){
        return empty;
      }
    }
  }

  // returns the number of entries
  intT count() {
    intT m = 0;
    for (intT i=0; i < m; i++)
      if (TA[i] != empty) m++;
    return m;
  }

  // returns all the current entries compacted into a sequence
  _seq<eType> entries() {
    eType *R = newA(eType,m);
    intT k = 0;
    for (intT i=0; i < m; i++)
      if (TA[i] != empty) R[k++] = TA[i];
    return _seq<eType>(R,k);
  }

  // prints the current entries along with the index they are stored at
  void print() {
    cout << "vals = ";
    for (intT i=0; i < m; i++) 
      if (TA[i] != empty)
	cout << i << ":" << TA[i] << ",";
    cout << endl;
  }
};


#endif // _A_HASH_INCLUDED
