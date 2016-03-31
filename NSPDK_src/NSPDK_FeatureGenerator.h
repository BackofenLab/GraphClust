/* -*- mode:c++ -*- */
#ifndef NSPDK_FEATUREGENERATOR_H
#define NSPDK_FEATUREGENERATOR_H

#include "Utility.h"
#include "FlagsService.h"
#include "GraphClass.h"
#include "Histogram.h"
#include "FeatureGenerator.h"
#include <limits>

using namespace std;

struct DebugClass {
  void Clear();
  void Output(ostream& out)const;
  void OutputGraphEncoding(ostream& out)const;
  void OutputFeatureEncoding(ostream& out)const;
  void OutputFunctorEncoding(ostream& out)const;
  void OutputBinEncoding(ostream& out)const;
  void OutputPlainEncoding(ostream& out)const;

  map<unsigned, string> mHashToPlainFeatureMap;
  map<unsigned, string> mHashToPlainSubgraphMap;
  map<unsigned, pair<unsigned,unsigned> > mHashToPlainSubgraphPairMap;
  map<unsigned, string> mHashToRadiusDistanceMap;
  map<unsigned, string> mHashToPredicateMap;
  map<unsigned, string> mHashToBinMap;

  multimap<string,string> mEncodingToPlainNeighborhoodMap;
};

class NSPDK_FeatureGenerator : public FeatureGenerator, public FlagsServiceClient {
public:
  NSPDK_FeatureGenerator(const std::string& id="");
  virtual void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList=vector<unsigned>());
  virtual double compute_kernel(const GraphClass& aG, const GraphClass& aH, const vector<unsigned>& aGFirstEndpointList=vector<unsigned>(), const vector<unsigned>& aHFirstEndpointList=vector<unsigned>());
  void OutputParameters(ostream& out)const;
  void Output(ostream& out)const;
  void OutputFeatureMap(ostream& out)const;
  virtual void Clear();
  void ClearCache(string aGraphID){}
protected:
  void ConvertToSparseVector(SVector& x) const;
  void GetFirstEndpoints(const GraphClass& aG, vector<unsigned>& oFirstEndpointRootId)const;
  virtual void GenerateFeatures(const GraphClass& aG, const vector<unsigned>& aFirstEndpointList=vector<unsigned>());
  virtual void GenerateFeatures(const GraphClass& aG, int aDistance, int aRadius, const vector<unsigned>& aFirstEndpointList);
  pair<unsigned,SecondOrderHistogramClass > BoundedRadiusRootedGraphCanonicalForm(int aRootVertexIndex, const GraphClass& aG, int aRadius);
  unsigned HardEncoding(int aRootVertexIndex, const GraphClass& aG, int aRadius);
  SecondOrderHistogramClass SoftEncoding(int aRootVertexIndex, const GraphClass& aG, int aRadius);
  string GetVertexLabel(const GraphClass& aG, unsigned aVertexID)const;
  string GetEdgeLabel(const GraphClass& aG, unsigned aSrcVertexID, unsigned aDestVertexID)const;
  void InsertFeature(const string& aLabel, int aRootVertexIndex, const GraphClass& aG, SecondOrderHistogramClass& oSoftAttributeList);
  unsigned Radius0RootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG);
  SecondOrderHistogramClass Radius0RootedGraphCanonicalFormAttributeList(int aRootVertexIndex, const GraphClass& aG);
  unsigned Radius1RootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG);
  SecondOrderHistogramClass  Radius1RootedGraphCanonicalFormAttributeList(int aRootVertexIndex, const GraphClass& aG);
  virtual unsigned RadiusKRootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG, int aRadius);
  virtual SecondOrderHistogramClass RadiusKRootedGraphCanonicalFormAttributeList(int aRootVertexIndex, const GraphClass& aG, int aRadius);
protected:
  unsigned mRadius;
  unsigned mDistance;
  string mMatchType;
  unsigned mHashBitSize;
  unsigned mBitMask;
  bool mMinKernel;
  bool mNormalization;
  unsigned mDebugVerbosity;

  map<unsigned,ThirdOrderHistogramClass> mFeatureList;
  map<pair<unsigned,unsigned>,pair<unsigned,SecondOrderHistogramClass > > mSubgraphEncodingCache;//semantic: <radius,vertex_id> \mapsto <hash_code,soft_attribute_histogram>
  
  mutable map<unsigned,string> mVertexToLabelMap;
  mutable map<pair<unsigned,unsigned>,string> mEdgeToLabelMap;

  mutable DebugClass mDebugInfo;
};

//----------------------------------------------------------------------------------------------------------------------------------
// Approximate NSPDK
class ANSPDK_FeatureGenerator : public NSPDK_FeatureGenerator{
public:
  ANSPDK_FeatureGenerator(const std::string& id);
protected:
  virtual unsigned RadiusKRootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG, int aRadius);
  virtual SecondOrderHistogramClass RadiusKRootedGraphCanonicalFormAttributeList(int aRootVertexIndex, const GraphClass& aG, int aRadius);
};

//----------------------------------------------------------------------------------------------------------------------------------
//3D NSPDK (suitable for molecules)
class NSPDK3D_FeatureGenerator : public ANSPDK_FeatureGenerator{
public:
  NSPDK3D_FeatureGenerator (const std::string& id);
  virtual void Clear();
  virtual void GenerateFeatures(const GraphClass& aG, int aDistance, int aRadius, const vector<unsigned>& aFirstEndpointList);
  vector<double> GetRootDirection(int aRootVertexIndex, const GraphClass& aG, int aRadius)const;
  vector<double> GetDirection(int aSrcVertexIndex, int aDestVertexIndex, const GraphClass& aG)const;
  unsigned ComputeDirectionAgreement(const vector<double>& aSrcDir,const vector<double>& aDestDir)const;
  vector<double> ComputeAverageDirection(const vector<vector<double> >& aDirectionList)const;
  double ComputeVectorNorm(const vector<double>& aVec)const;
public:
  map<pair<unsigned,unsigned>,vector<double> > mRootDirectionEncodingCache;//semantic: <radius,vertex_id> \mapsto 3D coordinates
};

//----------------------------------------------------------------------------------------------------------------------------------
// Memoized NSPDK
// feature information is stored for each graph-vertex pair; this saves computation time when only viewpoints are changed across feature generation invocations 
class MNSPDK_FeatureGenerator : public NSPDK_FeatureGenerator{
public:
  MNSPDK_FeatureGenerator(const std::string& id="");
  virtual void Clear();
  virtual void ClearCache(string aGraphID);
  virtual void GenerateFeatures(const GraphClass& aG, const vector<unsigned>& aFirstEndpointList);
protected:
  virtual void GenerateFeatures(const GraphClass& aG, int aDistance, int aRadius, const vector<unsigned>& aFirstEndpointList);
protected:
  map<string, map<unsigned , map < unsigned , pair<unsigned,SecondOrderHistogramClass > > > > mMemoizedSubgraphEncodingCache;//semantic: <graph_id,<radius_distance_hash,<vertices_id_hash>>> \mapsto <hash_code,soft_attribute_histogram>, where radius_distance_hash distinguishes between the case of neighbourhood features and path features and vertices_id_hash is only one vertex_id for neighbourhood and is a pair src_dest for path  
};

//----------------------------------------------------------------------------------------------------------------------------------
//Relational NSPDK
//vertex label of viewpoints is marked; this ensures that viewpoint identity influences the feature representation 
class RNSPDK_FeatureGenerator : public NSPDK_FeatureGenerator{
public:
  RNSPDK_FeatureGenerator(const std::string& id);
  virtual void GenerateFeatures(const GraphClass& aG, const vector<unsigned>& aFirstEndpointList=vector<unsigned>());
};

//----------------------------------------------------------------------------------------------------------------------------------
//ALias NSPDK
//the neighborhood of the first endpoint is discarded: this allows different types of vertices to be aliased as they become distinguishable only via their neighborhood
class ALNSPDK_FeatureGenerator : public ANSPDK_FeatureGenerator{
public:
  ALNSPDK_FeatureGenerator(const std::string& id);
  virtual void GenerateFeatures(const GraphClass& aG, int aDistance, int aRadius, const vector<unsigned>& aFirstEndpointList);
};

//----------------------------------------------------------------------------------------------------------------------------------
//Gapped NSPDK
//the distance information is discarded: pairs of neighbourhood subgraphs at different distances are encoded identically
class GNSPDK_FeatureGenerator : public NSPDK_FeatureGenerator{
public:
  GNSPDK_FeatureGenerator(const std::string& id);
  virtual void GenerateFeatures(const GraphClass& aG, int aDistance, int aRadius, const vector<unsigned>& aFirstEndpointList);
};

//----------------------------------------------------------------------------------------------------------------------------------
//Union of Shortest Thick Paths NSPDK
//for each pairs ov vertices within distance D the graph corresponding to the union of the shortest paths with thickness T is created
//the sparse vector representation for each such graph is computed
//the average vector is returned as final repersentation
class USTPNSPDK_FeatureGenerator : public NSPDK_FeatureGenerator{
public:
  USTPNSPDK_FeatureGenerator(const std::string& id);
  virtual void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList=vector<unsigned>());
  virtual void Clear();
protected:
  unsigned mThickness;
  unsigned mThicknessDistance;
  map<pair<unsigned,unsigned>,SVector> mFeatureVectorMap;
};

//----------------------------------------------------------------------------------------------------------------------------------
//Abstract NSPDK  
//out of all first endpoints select only those that are abstraction relations
//for each of them find the is_abstraction vertices A and is_part vertices P
//use A (P) as first endpoints in normal NSPDK and induce feature representation
//reduce feature vector dimension using MinHash technique
//make cartesian product vector V between reduced vector representation Va and reduced vector representation Vp
//return sum of all vectors V as new vector representation
class ABNSPDK_FeatureGenerator : public MNSPDK_FeatureGenerator{
public:
  ABNSPDK_FeatureGenerator(const std::string& id);
  virtual void generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList=vector<unsigned>());
protected:
  SVector MinHash(SVector& aX);
  void GetFirstEndpoints(const GraphClass& aG, vector<unsigned>& oFirstEndpointList)const;
  SVector ComputeCartesianProduct(SVector& x, SVector& z);
protected:
  unsigned mNumMinHashFunctions;
  unsigned mRelationRadius;
  unsigned mRelationDistance;
  double mRelationWeight;
  double mAbstractWeight;
};

#endif /* NSPDK_FEATUREGENERATOR_H */

/*! \mainpage Neighbourhood Subgraph Pairwise Distance Kernel
 *
 * \section intro_sec Introduction
 *
 * The NSPDK library offers a series of procedures to extract an explicit feature representation for generic graphs.
 * The code comprises 3 main parts: 1) Graph managment 2) Feature construction 3) Driver programs.
 *
 * \section part1 Part 1: Graph Managment
 * Classes: BaseGraphClass GraphClass
 *
 * BaseGraphClass provides methods to build a graph and query vertex and edge attributes.
 *
 * GraphClass provides an additional layer of semantic on top of BaseGraphClass. In addition it provides output formatting and pairwise distance query functions. 
 *  \section part2 Part 2: Feature Construction
 * Classes: FeatureGenerator
 *
 * FeatureGenerator provides the method generate_feature_vector that given a graph in input returns a sparse vector of features.
 *  \section part3 Part 3: Driver Programs
 *  \subsection discriminative Discriminative Linear Model: svmsgdnspdk
 *  \subsection clustering Clustering: NSPDK
 *  \subsection hclus Hyerarchical Clustering: KQuickShift
 */
