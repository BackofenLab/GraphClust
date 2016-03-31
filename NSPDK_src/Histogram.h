/* -*- mode:c++ -*- */
#ifndef HISTOGRAM_H
#define HISTOGRAM_H

using namespace std;
///Implements a histogram data structure, i.e. a container that allows either
///1) the insertion of integers (and counts the number of their occurencies)
///or 2) the insertion of pairs of integer (a bin identifier) with an associated real value
///(in this case we sum all the values).
///The other important member function is Add that allows to combine two objects into one
///(inserting all elements of the second object into the first object).
class HistogramClass{
  friend ostream& operator<<(ostream& out, const HistogramClass& aHC);
public:
  HistogramClass();
  HistogramClass(const HistogramClass& aH);
  void Insert(unsigned aBin);
  void Insert(unsigned aBin, double aValue);
  void Add(const HistogramClass& aH);
  ostream& Output(ostream& out)const;
  unsigned Size()const;
public:
  map<unsigned,double> mHistogram;
};

//---------------------------------------------------------------------------------
///Implements a histogram data structure that contains other histograms.
///The idea is to keep counts that are indexed with two levels.
class SecondOrderHistogramClass{
  friend ostream& operator<<(ostream& out, const SecondOrderHistogramClass& aSOH);
public:
  SecondOrderHistogramClass();
  SecondOrderHistogramClass(const SecondOrderHistogramClass& aSOH);
  void Insert(unsigned aSecondOrderBin, unsigned aBin);
  void Insert(unsigned aSecondOrderBin, unsigned aBin, double aValue);
  void Add(const SecondOrderHistogramClass& aSOH);
  void Add(unsigned aSecondOrderBin,const HistogramClass& aH);
  ostream& Output(ostream& out)const;
  unsigned Size()const;
public:
  map<unsigned,HistogramClass> mSecondOrderHistogram;
};

//---------------------------------------------------------------------------------
///Implements a histogram data structure that contains other histograms of histograms.
///The idea is to keep counts that are indexed with three levels.
class ThirdOrderHistogramClass{
  friend ostream& operator<<(ostream& out, const ThirdOrderHistogramClass& aTOH);
public:
  ThirdOrderHistogramClass();
  ThirdOrderHistogramClass(const ThirdOrderHistogramClass& aSOH);
  void Insert(unsigned aThirdOrderBin, unsigned aSecondOrderBin, unsigned aBin);
  void Insert(unsigned aThirdOrderBin, unsigned aSecondOrderBin, unsigned aBin, double aValue);
  void Add(const ThirdOrderHistogramClass& aTOH);
  void Add(unsigned aThirdOrderBin,const SecondOrderHistogramClass& aSOH);
  ostream& Output(ostream& out)const;
  unsigned Size()const;
public:
  map<unsigned,SecondOrderHistogramClass> mThirdOrderHistogram;
};

#endif
