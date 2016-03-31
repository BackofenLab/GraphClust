/* -*- mode:c++ -*- */
#ifndef BASEGRAPHCLASS_H
#define BASEGRAPHCLASS_H

#include <set>
#include <string>
using namespace std;

///Implements a list of vertices + list of edges + adjacency list data structure to represent graphs.
///Properties of vertices and edges are stored in 3 vectors of vectors: 1) symbolic (i.e. strings) 2) numeric (i.e. doubles) 3) status (i.e. boolean).
class BaseGraphClass{
  friend ostream& operator<<(ostream& out, const BaseGraphClass& aSG);

  ///Implements the edge as an ID and a destination vertex ID
  class EdgeClass{
  public:
    EdgeClass();
    EdgeClass(unsigned aDestVertexID,unsigned aEdgeID);
  public:
    unsigned mDestVertexID;
    unsigned mEdgeID;
  };
public:
  BaseGraphClass();
  void ResizeMemory();
  unsigned GetVertexInducedRootedSubGraph(const set<unsigned>& aVertexSet, unsigned aNominalRootIndex, BaseGraphClass& oG)const;
  unsigned InsertVertex();
  unsigned InsertEdge(unsigned aSrcVertexID, unsigned aDestVertexID);
  void SetVertexNumericAttributeList(unsigned aID,const vector<double>& aAttributeList);
  void SetVertexNumericAttributeList(unsigned aID,unsigned aAttributeID, double aValue);
  void SetVertexSymbolicAttributeList(unsigned aID,const vector<string>& aAttributeList);
  void SetVertexSymbolicAttributeList(unsigned aID,unsigned aAttributeID, const string& aValue);
  void SetVertexStatusAttributeList(unsigned aID,const vector<bool>& aAttributeList);
  void SetVertexStatusAttributeList(unsigned aID,unsigned aAttributeID, bool aValue);
  void SetEdgeNumericAttributeList(unsigned aID,const vector<double>& aAttributeList);
  void SetEdgeNumericAttributeList(unsigned aID,unsigned aAttributeID, double aValue);
  void SetEdgeNumericAttributeList(unsigned aSrcID, unsigned aDestID,unsigned aAttributeID, double aValue);
  void SetEdgeSymbolicAttributeList(unsigned aID,const vector<string>& aAttributeList);
  void SetEdgeSymbolicAttributeList(unsigned aID,unsigned aAttributeID, const string& aValue);
  void SetEdgeSymbolicAttributeList(unsigned aSrcID, unsigned aDestID,unsigned aAttributeID, const string& aValue);
  void SetEdgeStatusAttributeList(unsigned aID,const vector<bool>& aAttributeList);
  void SetEdgeStatusAttributeList(unsigned aID,unsigned aAttributeID, bool aValue);
  void SetEdgeStatusAttributeList(unsigned aSrcID, unsigned aDestID,unsigned aAttributeID, bool aValue);
  void SetVertexSymbolicID(unsigned aID,string aSID);
  string GetVertexSymbolicID(unsigned aID) const;
  vector<string> GetVertexSymbolicAttributeList(unsigned aID)const;
  vector<double> GetVertexNumericAttributeList(unsigned aID)const;
  vector<bool> GetVertexStatusAttributeList(unsigned aID)const;
  string GetVertexSymbolicAttributeList(unsigned aID, unsigned aAttributeID)const;
  double GetVertexNumericAttributeList(unsigned aID, unsigned aAttributeID)const;
  bool GetVertexStatusAttributeList(unsigned aID, unsigned aAttributeID)const;
  vector<string> GetEdgeSymbolicAttributeList(unsigned aSrcID, unsigned aDestID)const;
  string GetEdgeSymbolicAttributeList(unsigned aSrcID, unsigned aDestID, unsigned aAttributeID)const;
  vector<double> GetEdgeNumericAttributeList(unsigned aSrcID, unsigned aDestID)const;
  double GetEdgeNumericAttributeList(unsigned aSrcID, unsigned aDestID, unsigned aAttributeID)const;
  vector<bool> GetEdgeStatusAttributeList(unsigned aSrcID, unsigned aDestID)const;
  bool GetEdgeStatusAttributeList(unsigned aSrcID, unsigned aDestID, unsigned aAttributeID)const;
  vector<unsigned> GetVertexAdjacentList(unsigned aID)const;
  unsigned GetEdgeID(unsigned aSrcID, unsigned aDestID)const;
  unsigned GetEdgeSource(unsigned aEdgeID)const;
  unsigned GetEdgeDestination(unsigned aEdgeID)const;
  bool IsEdge(unsigned aSrcID, unsigned aDestID)const;
  ostream& Output(ostream& out)const;
  string Serialize()const;
  unsigned VertexSize()const;
  unsigned EdgeSize()const;
  bool IsEmpty()const;
protected:
  mutable bool mTopologicalChangeOccurrence;

  unsigned mVertexSize;
  unsigned mEdgeSize;

  vector<string> mVertexSymbolicIDList;
  vector< vector<EdgeClass> > mAdjacencyList;

  vector<vector<double> > mVertexNumericAttributeList; 
  vector<vector<string> > mVertexSymbolicAttributeList; 
  vector<vector<bool> > mVertexStatusAttributeList; 

  vector<vector<double> > mEdgeNumericAttributeList; 
  vector<vector<string> > mEdgeSymbolicAttributeList; 
  vector<vector<bool> > mEdgeStatusAttributeList; 
};








#endif
