#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <vector>
#include <stdexcept>
#include <map>
#include <set>
//#include <unordered_map> // need compiler option -std=c++0x
#include <tr1/unordered_map>
#include <queue>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cassert>
#include <ctime>
#include <iomanip>
#include <limits>

using namespace std;

///Implements a type converter via streams
template <class OutType, class InType>
OutType stream_cast(const InType & t)
{
 stringstream ss;
 ss << t; // first insert value to stream
 OutType result; // value will be converted to OutType
 ss >> result; // write value to result
 return result;
}

///Returns a random number between 0 and 1
inline double random01(){
  return (double)rand()/(double)RAND_MAX;
}

///Returns a random integer between 0 and the argument aMax
inline unsigned randomUnsigned(unsigned aMax){
  return (unsigned)(rand() % aMax);
}

///Return an integer hash value for a given input integer in a given domain range
inline int IntHashSimple(int key, int aModulo)
{
  key = ~key + (key << 15); // key = (key << 15) - key - 1;
  key = key ^ (key >> 12);
  key = key + (key << 2);
  key = key ^ (key >> 4);
  key = key * 2057; // key = (key + (key << 3)) + (key << 11);
  key = key ^ (key >> 16);
  return key%aModulo;
}

///Return an integer hash value for a given input integer in a given domain range given an additional seed to select the random hash function
inline int IntHash(int key, int aModulo, unsigned aSeed){
  const double A=sqrt(2)-1;
  return IntHashSimple(key*(aSeed+1)*A,aModulo);
}

unsigned RSHash(const string& str);
unsigned RSHash(const vector<unsigned>& aV);
unsigned APHash(const string& str);
unsigned APHash(const vector<unsigned>& aV);
unsigned HashFunc(const string& str, unsigned aBitMask=2147483647);
unsigned HashFunc(const vector<unsigned>& aList, unsigned aBitMask=2147483647);



//------------------------------------------------------------------------------------------------------------------------
///Returns the time between the creation and destruction of an object of TimerClass
class TimerClass {
public:
  TimerClass();
  ~TimerClass();
private:
  std::time_t start_sec;
  std::clock_t start;
};

//------------------------------------------------------------------------------------------------------------------------
///Plots the increase of an internal counter (which has to be
///explicitely increased by calling the Count member). Returns also
///the time elapsed between the creation and the destruction of an
///object of ProgressBar.
class ProgressBar{
public:
  ProgressBar(unsigned aStep=100);
  ~ProgressBar();
  void Begin();
  void Count();
  unsigned End();
private:
  unsigned mStep;
  unsigned mCounter;
  TimerClass mTimer;
};

//------------------------------------------------------------------------------------------------------------------------
///Implements safe access policies to a vector container and offers
///members to compute various statistical estimators.
class VectorClass{
  friend ostream& operator<<(ostream& out,const VectorClass& aV);
public:
  VectorClass();
  VectorClass(unsigned aSize);
  void operator=(const VectorClass& aVector);
  VectorClass(const VectorClass& aVector);
  VectorClass(const vector<double>& aVector);
  void Init(unsigned aSize);
  void Import(const string& aFileName);
  void Clear();
  unsigned Size()const;
  ostream& Output(ostream& out)const;
  void PushBack(double aValue);
  double& operator[](unsigned i);
  double operator[](unsigned i)const;
  double Sum()const;
  double Mean()const;
  double StandardDeviation()const;
  double Order(double aOrder)const;
  double Median()const;
  double MedianAbsoluteDifference()const;
  double Min()const;
  double Max()const;
  VectorClass RemoveNulls();
  ostream& OutputStatistics(ostream& out);
protected:
  vector<double> mV;
};


#endif
