//  ==============
//  QCDSF/ANALYSIS
//  ==============
//! Container for multi-dimensional arrays (derived from QDP++ qdp_multi.h)
//! \file
//! \version $Id: multi.h,v 1.1 2009-12-13 20:24:34 pleiter Exp $

#ifndef RESULT_MULTI_H
#define RESULT_MULTI_H

#include <iostream>
#include <stdlib.h>

using namespace std;

namespace qa
{

  template<class T> class multi1d
  {
    public:
      multi1d() : poke(false)
      { F = 0; n1 = 0; }

      explicit multi1d(int ns1) : poke(false)
      { F = 0; resize(ns1); }

      multi1d(T *f, int ns1) : poke(true)
      { F = f; n1 = ns1; }

      multi1d(const multi1d& s) : poke(false), n1(s.n1), F(0)
      {
        resize(n1);
        for(int i = 0; i < n1; ++i)
          F[i] = s.F[i];
      }

      ~multi1d() {
        if (!poke)
          if (F)
  	  delete[] F;
      }

      void resize(int ns1) { resize(*this, ns1); }

      int size() const { return n1; }
      int size1() const { return n1; }

      multi1d& operator=(const multi1d& s1)
      {
        if (size() != s1.size())		// a simple check avoids resizing always
          resize(s1.size());

        for(int i = 0; i < n1; ++i)
          F[i] = s1.F[i];

        return *this;
      }

      template<class T1>
      multi1d<T>& operator=(const T1& s1)
      {
        if (F == 0) {
  	exit(1);
        }

        for(int i=0; i < n1; ++i)
          F[i] = s1;

        return *this;
      }

      multi1d<T>& operator=(const T* s1)
      {
        if (F == 0) {
  	cerr << "multi1d: left hand side not initialized in =" << endl;
  	exit(1);
        }

        for(int i = 0; i < n1; ++i)
          F[i] = s1[i];

        return *this;
      }

      multi1d<T>& operator+=(const multi1d<T>& s1)
      {
        if (size() != s1.size()) {
  	cerr << "multi1d: Sizes incompatible in +=" << endl;
  	exit(1);
        }

        for(int i = 0; i < n1; ++i)
          F[i] += s1.F[i];

        return *this;
      }

      multi1d<T>& operator+=(const T& s1)
      {
        if (F == 0) {
  	cerr << "multi1d: left hand side not initialized in +=" << endl;
  	exit(1);
        }

        for(int i=0; i < n1; ++i)
          F[i] += s1;

        return *this;
      }

      multi1d<T>& operator-=(const multi1d<T>& s1)
      {
        if (size() != s1.size()) {
  	cerr << "multi1d: Sizes incompatible in -=" << endl;
  	exit(1);
        }

        for(int i = 0; i < n1; ++i)
          F[i] -= s1.F[i];

        return *this;
      }


      multi1d<T>& operator-=(const T& s1)
      {
        if (F == 0) {
  	cerr << "multi1d: left hand side not initialized in -=" << endl;
  	exit(1);
        }

        for(int i = 0; i < n1; ++i)
          F[i] -= s1;

        return *this;
      }


      multi1d<T>& operator*=(const multi1d<T>& s1)
      {
        if (size() != s1.size()) {
  	cerr << "multi1d: Sizes incompatible in *=" << endl;
  	exit(1);
        }

        for(int i = 0; i < n1; ++i)
          F[i] *= s1.F[i];

        return *this;
      }


      multi1d<T>& operator*=(const T& s1)
      {
        if (F == 0) {
  	cerr << "multi1d: left hand side not initialized in *=" << endl;
  	exit(1);
        }

        for(int i=0; i < n1; ++i)
          F[i] *= s1;

        return *this;
      }


      multi1d<T>& operator/=(const multi1d<T>& s1)
      {
        if (size() != s1.size()) {
  	cerr << "multi1d: Sizes incompatible in /=" << endl;
  	exit(1);
        }

        for(int i=0; i < n1; ++i)
          F[i] /= s1.F[i];

        return *this;
      }


      multi1d<T>& operator/=(const T& s1)
      {
        if (F == 0) {
  	cerr << "multi1d: left hand side not initialized in /=" << endl;
  	exit(1);
        }

        for(int i=0; i < n1; ++i)
          F[i] /= s1;
        return *this;
      }


      const T* slice() const { return F; }
      T& operator()(int i) { return F[i]; }
      const T& operator()(int i) const { return F[i]; }
      const T& operator[](int i) const { return F[i]; }
      T& operator[](int i) { return F[i]; }

    private:
      template<typename I>
      void resize(multi1d<I>& disambiguator, int ns1)
      {
        if (F)
          delete [] F;
        n1 = ns1;
        F = new T[n1];
      }
      bool poke;
      int n1;
      T *F;
  };



  template<typename T>
  inline bool operator==(const multi1d<T>& n1, const multi1d<T>& n2)
  {
    if (n1.size() == 0 || n1.size() != n2.size())
      return false;

    for(int i = 0; i < n1.size(); ++i)
      if (n2[i] != n1[i])
        return false;

    return true;
  }


  template<typename T>
  inline bool operator!=(const multi1d<T>& n1, const multi1d<T>& n2)
  {
    return ! (n1 == n2);
  }



  //! a > b
  /*! This definition follows that of string comparison */
  template<typename T>
  inline bool operator>(const multi1d<T>& a, const multi1d<T>& b)
  {
    bool ret = false;
    int  len = (a.size() < b.size()) ? a.size() : b.size();

    for(int i=0; i < len; ++i)
    {
      if (a[i] != b[i])
        return (a[i] > b[i]) ? true : false;
    }

    return (a.size() == b.size()) ? false : (a.size() > b.size()) ? true : false;
  }


  //! a <= b
  /*! This definition follows that of string comparison */
  template<typename T>
  inline bool operator<=(const multi1d<T>& a, const multi1d<T>& b)
  {
    return (a < b) || (a == b);
  }


  //! a >= b
  /*! This definition follows that of string comparison */
  template<typename T>
  inline bool operator>=(const multi1d<T>& a, const multi1d<T>& b)
  {
    return (a > b) || (a == b);
  }


  //--------------------------------------------------------------------------------------------------------------------
  // Basic math support
  //--------------------------------------------------------------------------------------------------------------------

  //! add Arrays
  template< typename T> 
  inline
  multi1d<T> operator+(const multi1d<T>& a, const multi1d<T>& b)
  {
    multi1d<T> c(a); 
    c+=b;
    return c;
  }

  //! subtract Arrays
  template< typename T> 
  inline
  multi1d<T> operator-(const multi1d<T>& a, const multi1d<T>& b)
  {
    multi1d<T> c(a); 
    c-=b;
    return c;
  }

  //! multiply Arrays
  template< typename T> 
  inline
  multi1d<T> operator*(const multi1d<T>& a, const multi1d<T>& b)
  {
    multi1d<T> c(a); 
    c*=b;
    return c;
  }

  //!divide Arrays
  template< typename T> 
  inline
  multi1d<T> operator/(const multi1d<T>& a, const multi1d<T>& b)
  {
    multi1d<T> c(a); 
    c/=b;
    return c;
  }

  //! scalar + Array
  template< typename T> 
  inline
  multi1d<T> operator+(const T& s, const multi1d<T>& a)
  {
    multi1d<T> c(a); 
    c+=s;
    return c;
  }

  //! Array + scalar
  template< typename T> 
  inline
  multi1d<T> operator+(const multi1d<T>& a, const T& s)
  {
    multi1d<T> c(a); 
    c+=s;
    return c;
  }

  //! scalar - Array
  template< typename T> 
  inline
  multi1d<T> operator-(const T& s, const multi1d<T>& a)
  {
    multi1d<T> c(-a); 
    c+=s;
    return c;
  }
  //! Array - scalar
  template< typename T> 
  inline
  multi1d<T> operator-(const multi1d<T>& a, const T& s)
  {
    multi1d<T> c(a); 
    c-=s;
    return c;
  }

  //! scalar * Array
  template< typename T> 
  inline
  multi1d<T> operator*(const T& s, const multi1d<T>& a)
  {
    multi1d<T> c(a); 
    c*=s;
    return c;
  }

  //! Array * scalar
  template< typename T> 
  inline
  multi1d<T> operator*(const multi1d<T>& a, const T& s)
  {
    multi1d<T> c(a); 
    c*=s;
    return c;
  }

  //! scalar / Array
  template< typename T> 
  inline
  multi1d<T> operator/(const T& s, const multi1d<T>& a)
  {
    multi1d<T> c(a.size());
    c = s;
    c/= a;
    return c;
  }

  //! Array / scalar
  template< typename T> 
  inline
  multi1d<T> operator/(const multi1d<T>& a, const T& s)
  {
    multi1d<T> c(a); 
    c/=s;
    return c;
  }

  //! sqrt
  template< typename T> 
  inline
  multi1d<T> sqrt(const multi1d<T>& a)
  {
    multi1d<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = sqrt(a[i]);
    }
    return c;
  }

  //! log
  template< typename T> 
  inline
  multi1d<T> log(const multi1d<T>& a)
  {
    multi1d<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = log(a[i]);
    }
    return c;
  }

  //! sin
  template< typename T> 
  inline
  multi1d<T> sin(const multi1d<T>& a)
  {
    multi1d<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = sin(a[i]);
    }
    return c;
  }


  //! cos
  template< typename T> 
  inline
  multi1d<T> cos(const multi1d<T>& a)
  {
    multi1d<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = cos(a[i]);
    }
    return c;
  }

  //! tan
  template< typename T> 
  inline
  multi1d<T> tan(const multi1d<T>& a)
  {
    multi1d<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = tan(a[i]);
    }
    return c;
  }

  //! asin
  template< typename T> 
  inline
  multi1d<T> asin(const multi1d<T>& a)
  {
    multi1d<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = asin(a[i]);
    }
    return c;
  }


  //! acos
  template< typename T> 
  inline
  multi1d<T> acos(const multi1d<T>& a)
  {
    multi1d<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = acos(a[i]);
    }
    return c;
  }

  //! atan
  template< typename T> 
  inline
  multi1d<T> atan(const multi1d<T>& a)
  {
    multi1d<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = atan(a[i]);
    }
    return c;
  }



  //! norm2 of an array
  template< typename T> 
  inline
  T norm2(const multi1d<T>& a)
  {
    T nn = a[0]*a[0];  // assumes at least 1 element
    for(int i=1; i < a.size(); ++i)
      nn += a[i]*a[i];

    return nn;
  }


  template<class T> class multi2d
  {
  public:
    multi2d() {F=0;n1=n2=sz=0;}
    multi2d(T *f, int ns2, int ns1) {F=f; n1=ns1; n2=ns2; sz=n1*n2;}
    explicit multi2d(int ns2, int ns1) {F=0;resize(ns2,ns1);}
    ~multi2d() { if (F) delete[] F;}

    //! Copy constructor
    multi2d(const multi2d& s): n1(s.n1), n2(s.n2), sz(s.sz), F(0)
      {
        resize(n2,n1);
        for(int i=0; i < sz; ++i)
  	F[i] = s.F[i];
      }

    //! Allocate mem for the array
    void resize(int ns2, int ns1) {
      if (F)
        delete[] F; 
      n1=ns1; 
      n2=ns2;  
      sz=n1*n2; 
      F = new T[sz];
    }

    //! Size of array
    int size1() const {return n1;}
    int size2() const {return n2;}

    //! Equal operator uses underlying = of T
    multi2d<T>& operator=(const multi2d<T>& s1)
      {
        resize(s1.size2(), s1.size1());   // always resize

        for(int i=0; i < sz; ++i)
  	F[i] = s1.F[i];
        return *this;
      }


    template<class T1>
    multi2d<T>& operator=(const T1& s1)
      {
        if (F == 0)
        {
  	cerr << "multi2d: left hand side not initialized in =" << endl;
  	exit(1);
        }

        for(int i=0; i < sz; ++i)
  	F[i] = s1;
        return *this;
      }


    const T* slice(int j) const {return F+n1*j;}
    T& operator()(int j, int i) {return F[i+n1*j];}
    const T& operator()(int j, int i) const {return F[i+n1*j];}
    multi1d<T> operator[](int j) 
    {
        return multi1d<T>(F+j*n1,n1);
    }
    const multi1d<T> operator[](int j) const 
    {
        return multi1d<T>(F+j*n1,n1);
    }

  private:
    int n1;
    int n2;
    int sz;
    T *F;
  };

}

#endif
