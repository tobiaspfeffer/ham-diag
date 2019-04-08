#ifndef MATRIX_H
#define MATRIX_H


#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <cstdlib>

template<class T>
class matrix 
{
  int n; // size of matrix  (size of first index)
  int m; // size of matrix (size of second index)
  

public:
  std::vector<std::vector<T>> a;
  explicit matrix(const int size, const int size2): n(size),m(size2), a(n,std::vector<T>(m,0)) {}
  
  T& operator()(const int x, const int y) {return a[x][y];}
  const T& operator()(const int x, const int y) const {return a[x][y];}
  
  std::vector<T>& operator()(const int x) {return a[x];}
  const std::vector<T>& operator()(const int x) const {return a[x];}
  
  matrix& operator=(const matrix& rhs)
  {
    if(this != &rhs)
    {
      a = rhs.a;
      n = rhs.n;
      m = rhs.m;
    }
    return *this;
  }
};

#endif
