#include "Utility.h"
#include "BaseGraphClass.h"

using namespace std;
//---------------------------------------------------------------------------------
BaseGraphClass::EdgeClass::EdgeClass():mDestVertexID(numeric_limits<unsigned>::max()),mEdgeID(numeric_limits<unsigned>::max()){}
BaseGraphClass::EdgeClass::EdgeClass(unsigned aDestVertexID,unsigned aEdgeID):mDestVertexID(aDestVertexID),mEdgeID(aEdgeID){}

//---------------------------------------------------------------------------------
ostream& operator<<(ostream& out, const BaseGraphClass& aSG){aSG.Output(out);return out;}
BaseGraphClass::BaseGraphClass():mTopologicalChangeOccurrence(true),mVertexSize(0),mEdgeSize(0){}

void BaseGraphClass::ResizeMemory(){
  //NOTE: using the shrink_to_fit() technique as in http://www.gotw.ca/gotw/054.htm
  //for all vertices
  vector<vector<double> >(mVertexNumericAttributeList).swap(mVertexNumericAttributeList);
  vector<vector<string> >(mVertexSymbolicAttributeList).swap(mVertexSymbolicAttributeList);
  vector<vector<bool> >(mVertexStatusAttributeList).swap(mVertexStatusAttributeList);

  //for all edges
  vector<vector<double> >(mEdgeNumericAttributeList).swap(mEdgeNumericAttributeList);
  vector<vector<string> >(mEdgeSymbolicAttributeList).swap(mEdgeSymbolicAttributeList);
  vector<vector<bool> >(mEdgeStatusAttributeList).swap(mEdgeStatusAttributeList);

  //for adjacency structure
  vector< vector<EdgeClass> > (mAdjacencyList).swap(mAdjacencyList);
}

unsigned BaseGraphClass::GetEdgeSource(unsigned aEdgeID)const{
  unsigned src=0;
  for (unsigned i=0;i<VertexSize();i++)
    for (unsigned j=0; j<mAdjacencyList[i].size();j++)
      if (mAdjacencyList[i][j].mEdgeID==aEdgeID) {
	src=i;
	return src;
      }
  throw range_error("Edge with id: "+stream_cast<string>(aEdgeID)+" does not exist");  
  return src;
}

unsigned BaseGraphClass::GetEdgeDestination(unsigned aEdgeID)const{
  unsigned dest=0;
  for (unsigned i=0;i<VertexSize();i++)
    for (unsigned j=0; j<mAdjacencyList[i].size();j++)
      if (mAdjacencyList[i][j].mEdgeID==aEdgeID) {
	dest=mAdjacencyList[i][j].mDestVertexID;;
	return dest;
      }
  throw range_error("Edge with id: "+stream_cast<string>(aEdgeID)+" does not exist");  
  return dest;  
}

unsigned BaseGraphClass::GetVertexInducedRootedSubGraph(const set<unsigned>& aVertexSet, unsigned aNominalRootIndex, BaseGraphClass& oG)const{
    map<unsigned,unsigned> index_map_nominal_to_real;
    for (set<unsigned>::const_iterator it=aVertexSet.begin();it!=aVertexSet.end();++it){
      unsigned nominal_index=*it;
      unsigned real_index=oG.InsertVertex();
      index_map_nominal_to_real[nominal_index]=real_index;
      oG.SetVertexNumericAttributeList(real_index, mVertexNumericAttributeList[nominal_index]);
      oG.SetVertexSymbolicAttributeList(real_index, mVertexSymbolicAttributeList[nominal_index]);
      oG.SetVertexStatusAttributeList(real_index, mVertexStatusAttributeList[nominal_index]);
    }
    unsigned real_root_index=index_map_nominal_to_real[aNominalRootIndex];
    for (set<unsigned>::const_iterator it=aVertexSet.begin();it!=aVertexSet.end();++it){
      unsigned u=*it;
      unsigned nominal_src_index=u;
      for (unsigned v=0;v<mAdjacencyList[u].size();++v){	
	unsigned nominal_dest_index=mAdjacencyList[u][v].mDestVertexID;
	//if dest vertex is among the input set then add corresponding edge
	if (aVertexSet.count(nominal_dest_index)>0) {
	  unsigned real_src_index=index_map_nominal_to_real[nominal_src_index];
	  unsigned real_dest_index=index_map_nominal_to_real[nominal_dest_index];
	  unsigned nominal_edge_index=mAdjacencyList[u][v].mEdgeID;
	  unsigned real_edge_index=oG.InsertEdge(real_src_index,real_dest_index);
	  oG.SetEdgeNumericAttributeList(real_edge_index,mEdgeNumericAttributeList[nominal_edge_index]);
	  oG.SetEdgeSymbolicAttributeList(real_edge_index,mEdgeSymbolicAttributeList[nominal_edge_index]);
	  oG.SetEdgeStatusAttributeList(real_edge_index,mEdgeStatusAttributeList[nominal_edge_index]);
	}
      }
    }
    return real_root_index;
  }

unsigned BaseGraphClass::InsertVertex(){
    unsigned vertex_size=mVertexSize;
    mVertexSymbolicIDList.push_back("");
    mAdjacencyList.push_back(vector<EdgeClass>());
    mVertexNumericAttributeList.push_back(vector<double>());
    mVertexSymbolicAttributeList.push_back(vector<string>());
    mVertexStatusAttributeList.push_back(vector<bool>());
    mVertexSize++;
    mTopologicalChangeOccurrence=true;
    return vertex_size;
  }

unsigned BaseGraphClass::InsertEdge(unsigned aSrcVertexID, unsigned aDestVertexID){
    if (aSrcVertexID>=mVertexSize || aDestVertexID>=mVertexSize) 
      throw range_error("Edge between non existing vertices: "+stream_cast<string>(aSrcVertexID)+" "+stream_cast<string>(aDestVertexID));
    unsigned edge_size=mEdgeSize;
    mAdjacencyList[aSrcVertexID].push_back(EdgeClass(aDestVertexID,edge_size));
    mEdgeNumericAttributeList.push_back(vector<double>());
    mEdgeSymbolicAttributeList.push_back(vector<string>());
    mEdgeStatusAttributeList.push_back(vector<bool>());
    mEdgeSize++;
    mTopologicalChangeOccurrence=true;
    return edge_size;
  }

void BaseGraphClass::SetVertexSymbolicID(unsigned aID,string aSID) {
  mVertexSymbolicIDList[aID] = aSID;
}
string BaseGraphClass::GetVertexSymbolicID(unsigned aID) const {
  return mVertexSymbolicIDList[aID];
}  

void BaseGraphClass::SetVertexNumericAttributeList(unsigned aID,const vector<double>& aAttributeList) {
    if (aID>=mVertexSize) throw range_error("Setting vertex attribute for non existing vertex: "+stream_cast<string>(aID));
    mVertexNumericAttributeList[aID]=aAttributeList;
  }
void BaseGraphClass::SetVertexNumericAttributeList(unsigned aID,unsigned aAttributeID, double aValue){
  if (aID>=mVertexSize) throw range_error("Setting vertex attribute for non existing vertex: "+stream_cast<string>(aID));
  if (aAttributeID<mVertexNumericAttributeList[aID].size())
    mVertexNumericAttributeList[aID][aAttributeID]=aValue;
  else {
    unsigned size=aAttributeID+1;
    vector<double> attribute_list(size,0);
    for (unsigned i=0;i<mVertexNumericAttributeList[aID].size();++i)
      attribute_list[i]=mVertexNumericAttributeList[aID][i];
    attribute_list[aAttributeID]=aValue;
    mVertexNumericAttributeList[aID]=attribute_list;
  }
}
void BaseGraphClass::SetVertexSymbolicAttributeList(unsigned aID,const vector<string>& aAttributeList){
    if (aID>=mVertexSize) throw range_error("Setting vertex attribute for non existing vertex: "+stream_cast<string>(aID));
    mVertexSymbolicAttributeList[aID]=aAttributeList;
  }
void BaseGraphClass::SetVertexSymbolicAttributeList(unsigned aID,unsigned aAttributeID, const string& aValue){
  if (aID>=mVertexSize) throw range_error("Setting vertex attribute for non existing vertex: "+stream_cast<string>(aID));
  if (aAttributeID<mVertexSymbolicAttributeList[aID].size())
    mVertexSymbolicAttributeList[aID][aAttributeID]=aValue;
  else {
    unsigned size=aAttributeID+1;
    vector<string> attribute_list(size,"");
    for (unsigned i=0;i<mVertexSymbolicAttributeList[aID].size();++i)
      attribute_list[i]=mVertexSymbolicAttributeList[aID][i];
    attribute_list[aAttributeID]=aValue;
    mVertexSymbolicAttributeList[aID]=attribute_list;
  }
}
void BaseGraphClass::SetVertexStatusAttributeList(unsigned aID,const vector<bool>& aAttributeList){
  mTopologicalChangeOccurrence=true; 
  if (aID>=mVertexSize) throw range_error("Setting vertex attribute for non existing vertex: "+stream_cast<string>(aID));
  mVertexStatusAttributeList[aID]=aAttributeList;
  }
void BaseGraphClass::SetVertexStatusAttributeList(unsigned aID,unsigned aAttributeID, bool aValue){
  mTopologicalChangeOccurrence=true;
  if (aID>=mVertexSize) throw range_error("Setting vertex attribute for non existing vertex: "+stream_cast<string>(aID));
  if (aAttributeID<mVertexStatusAttributeList[aID].size())
    mVertexStatusAttributeList[aID][aAttributeID]=aValue;
  else {
    unsigned size=aAttributeID+1;
    vector<bool> attribute_list(size,false);
    for (unsigned i=0;i<mVertexStatusAttributeList[aID].size();++i)
      attribute_list[i]=mVertexStatusAttributeList[aID][i];
    attribute_list[aAttributeID]=aValue;
    mVertexStatusAttributeList[aID]=attribute_list;
  }
}
vector<string> BaseGraphClass::GetVertexSymbolicAttributeList(unsigned aID)const{
  if (aID>=mVertexSize) throw range_error("Getting vertex attributes for non existing vertex: "+stream_cast<string>(aID));
  return mVertexSymbolicAttributeList[aID];
}
vector<double> BaseGraphClass::GetVertexNumericAttributeList(unsigned aID)const{
  if (aID>=mVertexSize) throw range_error("Getting vertex attributes for non existing vertex: "+stream_cast<string>(aID));
    return mVertexNumericAttributeList[aID];
  }
vector<bool> BaseGraphClass::GetVertexStatusAttributeList(unsigned aID)const{
  if (aID>=mVertexSize) throw range_error("Getting vertex attributes for non existing vertex: "+stream_cast<string>(aID));
    return mVertexStatusAttributeList[aID];
  }

string BaseGraphClass::GetVertexSymbolicAttributeList(unsigned aID, unsigned aAttributeID)const{
  if (aID>=mVertexSize) throw range_error("Getting vertex attributes for non existing vertex: "+stream_cast<string>(aID));
  if (aAttributeID>mVertexSymbolicAttributeList[aID].size()) throw range_error("Getting vertex attributes for non existing attribute id: "+stream_cast<string>(aAttributeID));
    return mVertexSymbolicAttributeList[aID][aAttributeID];
  }
double BaseGraphClass::GetVertexNumericAttributeList(unsigned aID, unsigned aAttributeID)const{
  if (aID>=mVertexSize) throw range_error("Getting vertex attributes for non existing vertex: "+stream_cast<string>(aID));
  if (aAttributeID>mVertexNumericAttributeList[aID].size()) throw range_error("Getting vertex attributes for non existing attribute id: "+stream_cast<string>(aAttributeID));
    return mVertexNumericAttributeList[aID][aAttributeID];
  }
bool BaseGraphClass::GetVertexStatusAttributeList(unsigned aID, unsigned aAttributeID)const{
  if (aID>=mVertexSize) throw range_error("Getting vertex attributes for non existing vertex: "+stream_cast<string>(aID));
  if (aAttributeID>mVertexStatusAttributeList[aID].size()) throw range_error("Getting vertex attributes for non existing attribute id: "+stream_cast<string>(aAttributeID));
    return mVertexStatusAttributeList[aID][aAttributeID];
  }


void BaseGraphClass::SetEdgeNumericAttributeList(unsigned aID,const vector<double>& aAttributeList) {
    if (aID>=mEdgeSize) throw range_error("Setting edge attribute for non existing edge: "+stream_cast<string>(aID));
    mEdgeNumericAttributeList[aID]=aAttributeList;
  }

void BaseGraphClass::SetEdgeNumericAttributeList(unsigned aID,unsigned aAttributeID, double aValue){
  if (aID>=mEdgeSize) throw range_error("Setting edge attribute for non existing edge: "+stream_cast<string>(aID));
  if (aAttributeID<mEdgeNumericAttributeList[aID].size())
    mEdgeNumericAttributeList[aID][aAttributeID]=aValue;
  else {
    unsigned size=aAttributeID+1;
    vector<double> attribute_list(size,0);
    for (unsigned i=0;i<mEdgeNumericAttributeList[aID].size();++i)
      attribute_list[i]=mEdgeNumericAttributeList[aID][i];
    attribute_list[aAttributeID]=aValue;
    mEdgeNumericAttributeList[aID]=attribute_list;
  }
}

void BaseGraphClass::SetEdgeNumericAttributeList(unsigned aSrcID, unsigned aDestID,unsigned aAttributeID, double aValue){
    unsigned edge_id=GetEdgeID(aSrcID,aDestID);
    SetEdgeNumericAttributeList(edge_id,aAttributeID, aValue);
  }
void BaseGraphClass::SetEdgeSymbolicAttributeList(unsigned aID,const vector<string>& aAttributeList){
    if (aID>=mEdgeSize) throw range_error("Setting edge attribute for non existing edge: "+stream_cast<string>(aID));
    mEdgeSymbolicAttributeList[aID]=aAttributeList;
  }
void BaseGraphClass::SetEdgeSymbolicAttributeList(unsigned aID,unsigned aAttributeID, const string& aValue){
  if (aID>=mEdgeSize) throw range_error("Setting edge attribute for non existing edge: "+stream_cast<string>(aID));
  if (aAttributeID<mEdgeSymbolicAttributeList[aID].size())
    mEdgeSymbolicAttributeList[aID][aAttributeID]=aValue;
  else {
    unsigned size=aAttributeID+1;
    vector<string> attribute_list(size,"");
    for (unsigned i=0;i<mEdgeSymbolicAttributeList[aID].size();++i)
      attribute_list[i]=mEdgeSymbolicAttributeList[aID][i];
    attribute_list[aAttributeID]=aValue;
    mEdgeSymbolicAttributeList[aID]=attribute_list;
  }
}
void BaseGraphClass::SetEdgeSymbolicAttributeList(unsigned aSrcID, unsigned aDestID, unsigned aAttributeID, const string& aValue){
    unsigned edge_id=GetEdgeID(aSrcID,aDestID);
    SetEdgeSymbolicAttributeList(edge_id,aAttributeID, aValue);
}
void BaseGraphClass::SetEdgeStatusAttributeList(unsigned aID,const vector<bool>& aAttributeList){
  mTopologicalChangeOccurrence=true;
  if (aID>=mEdgeSize) throw range_error("Setting edge attribute for non existing edge: "+stream_cast<string>(aID));
  mEdgeStatusAttributeList[aID]=aAttributeList;
}
void BaseGraphClass::SetEdgeStatusAttributeList(unsigned aID,unsigned aAttributeID, bool aValue){
  mTopologicalChangeOccurrence=true;
  if (aID>=mEdgeSize) throw range_error("Setting edge attribute for non existing edge: "+stream_cast<string>(aID));
  if (aAttributeID<mEdgeStatusAttributeList[aID].size())
    mEdgeStatusAttributeList[aID][aAttributeID]=aValue;
  else {
    unsigned size=aAttributeID+1;
    vector<bool> attribute_list(size,false);
    for (unsigned i=0;i<mEdgeStatusAttributeList[aID].size();++i)
      attribute_list[i]=mEdgeStatusAttributeList[aID][i];
    attribute_list[aAttributeID]=aValue;
    mEdgeStatusAttributeList[aID]=attribute_list;
  }
}
void BaseGraphClass::SetEdgeStatusAttributeList(unsigned aSrcID, unsigned aDestID, unsigned aAttributeID, bool aValue){
  unsigned edge_id=GetEdgeID(aSrcID,aDestID);
  SetEdgeStatusAttributeList(edge_id,aAttributeID,aValue);
}
vector<string> BaseGraphClass::GetEdgeSymbolicAttributeList(unsigned aSrcID, unsigned aDestID)const{
  unsigned edge_id=GetEdgeID(aSrcID,aDestID);
  return mEdgeSymbolicAttributeList[edge_id];
}
string BaseGraphClass::GetEdgeSymbolicAttributeList(unsigned aSrcID, unsigned aDestID, unsigned aAttributeID)const{
  unsigned edge_id=GetEdgeID(aSrcID,aDestID);
  if (aAttributeID>=mEdgeSymbolicAttributeList[edge_id].size())
    throw range_error("Getting edge attribute for non existing attribute id: "+stream_cast<string>(aAttributeID));
  return mEdgeSymbolicAttributeList[edge_id][aAttributeID];
}
vector<double> BaseGraphClass::GetEdgeNumericAttributeList(unsigned aSrcID, unsigned aDestID)const{
  unsigned edge_id=GetEdgeID(aSrcID,aDestID);
  return mEdgeNumericAttributeList[edge_id];
  }
double BaseGraphClass::GetEdgeNumericAttributeList(unsigned aSrcID, unsigned aDestID, unsigned aAttributeID)const{
  unsigned edge_id=GetEdgeID(aSrcID,aDestID);
  if (aAttributeID>=mEdgeNumericAttributeList[edge_id].size())
    throw range_error("Getting edge attribute for non existing attribute id: "+stream_cast<string>(aAttributeID));
  return mEdgeNumericAttributeList[edge_id][aAttributeID];
}
vector<bool> BaseGraphClass::GetEdgeStatusAttributeList(unsigned aSrcID, unsigned aDestID)const{
  unsigned edge_id=GetEdgeID(aSrcID,aDestID);
  return mEdgeStatusAttributeList[edge_id];
}
bool BaseGraphClass::GetEdgeStatusAttributeList(unsigned aSrcID, unsigned aDestID, unsigned aAttributeID)const{
  unsigned edge_id=GetEdgeID(aSrcID,aDestID);
  if (aAttributeID>=mEdgeStatusAttributeList[edge_id].size())
    throw range_error("Getting edge attribute for non existing attribute id: "+stream_cast<string>(aAttributeID));
  return mEdgeStatusAttributeList[edge_id][aAttributeID];
}
bool BaseGraphClass::IsEdge(unsigned aSrcID, unsigned aDestID)const{
  for (unsigned j=0; j<mAdjacencyList[aSrcID].size();j++)
    if (mAdjacencyList[aSrcID][j].mDestVertexID==aDestID) return true;
  return false;
}
unsigned BaseGraphClass::GetEdgeID(unsigned aSrcID, unsigned aDestID)const{
    for (unsigned j=0; j<mAdjacencyList[aSrcID].size();j++)
      if (mAdjacencyList[aSrcID][j].mDestVertexID==aDestID) {
	unsigned edge_id=mAdjacencyList[aSrcID][j].mEdgeID;
	return edge_id;
      }
    throw range_error("Edge between "+stream_cast<string>(aSrcID)+" and "+stream_cast<string>(aDestID)+" does not exist");
    return 0;
  }
vector<unsigned> BaseGraphClass::GetVertexAdjacentList(unsigned aID)const{
    vector<unsigned> adjacent_list;
    for (unsigned i=0;i<mAdjacencyList[aID].size();++i)
      adjacent_list.push_back(mAdjacencyList[aID][i].mDestVertexID);
    return adjacent_list;
  }
ostream& BaseGraphClass::Output(ostream& out)const{
    out<<"Graph adjacency list ("<<mAdjacencyList.size()<<")"<<endl;
    for (unsigned i=0;i<mAdjacencyList.size();++i){
      out<<i<<" ";
      for (unsigned j=0;j<mAdjacencyList[i].size();++j)
	out<<mAdjacencyList[i][j].mDestVertexID<<" ";
      out<<endl;      
    }

    out<<"Vertex Numeric Attribute List ("<<mVertexNumericAttributeList.size()<<")"<<endl;
    for (unsigned i=0;i<mVertexNumericAttributeList.size();++i){
      out<<i<<" ";
      for (unsigned j=0;j<mVertexNumericAttributeList[i].size();++j)
	out<<j<<":"<<mVertexNumericAttributeList[i][j]<<" ";
      out<<endl;
    }
      
    out<<"Vertex Symbolic Attribute List ("<<mVertexSymbolicAttributeList.size()<<")"<<endl;
    for (unsigned i=0;i<mVertexSymbolicAttributeList.size();++i){
      out<<i<<" ";
      for (unsigned j=0;j<mVertexSymbolicAttributeList[i].size();++j)
	out<<j<<":"<<mVertexSymbolicAttributeList[i][j]<<" ";
      out<<endl;
    }
      
    out<<"Vertex Status Attribute List ("<<mVertexStatusAttributeList.size()<<")"<<endl;
    for (unsigned i=0;i<mVertexStatusAttributeList.size();++i){
      out<<i<<" ";
      for (unsigned j=0;j<mVertexStatusAttributeList[i].size();++j)
	out<<j<<":"<<mVertexStatusAttributeList[i][j]<<" ";
      out<<endl;
    }
      
    out<<"Edge Numeric Attribute List ("<<mEdgeNumericAttributeList.size()<<")"<<endl;
    for (unsigned i=0;i<mEdgeNumericAttributeList.size();++i){
      out<<i<<" ";
      for (unsigned j=0;j<mEdgeNumericAttributeList[i].size();++j)
	out<<j<<":"<<mEdgeNumericAttributeList[i][j]<<" ";
      out<<endl;
    }
      
    out<<"Edge Symbolic Attribute List ("<<mEdgeSymbolicAttributeList.size()<<")"<<endl;
    for (unsigned i=0;i<mEdgeSymbolicAttributeList.size();++i){
      out<<i<<" ";
      for (unsigned j=0;j<mEdgeSymbolicAttributeList[i].size();++j)
	out<<j<<":"<<mEdgeSymbolicAttributeList[i][j]<<" ";
      out<<endl;
    }
      
    out<<"Edge Status Attribute List ("<<mEdgeStatusAttributeList.size()<<")"<<endl;
    for (unsigned i=0;i<mEdgeStatusAttributeList.size();++i){
      out<<i<<" ";
      for (unsigned j=0;j<mEdgeStatusAttributeList[i].size();++j)
	out<<j<<":"<<mEdgeStatusAttributeList[i][j]<<" ";
      out<<endl;
    }
      
    return out;
  }

string BaseGraphClass::Serialize()const{
    string encoding;
    encoding+=stream_cast<string>(mVertexSize)+" ";
    for (unsigned i=0;i<mVertexSize;++i){
      string vlabel=mVertexSymbolicAttributeList[i][0];//NOTE: output only first symbolic attribute as vertex label 
      encoding+=vlabel+" ";
    }
    encoding+=stream_cast<string>(mEdgeSize)+" ";
    for (unsigned u=0;u<mVertexSize;++u){
      for (unsigned j=0; j<mAdjacencyList[u].size();j++){
	unsigned v=mAdjacencyList[u][j].mDestVertexID;
	unsigned eid=mAdjacencyList[u][j].mEdgeID;
	string elabel=mEdgeSymbolicAttributeList[eid][0];//NOTE: output only first symbolic attribute as edge label 
	encoding+=stream_cast<string>(u)+" "+stream_cast<string>(v)+" "+elabel+" ";
      }
    }
    return encoding;
  }
unsigned BaseGraphClass::VertexSize()const{return mVertexSize;}
unsigned BaseGraphClass::EdgeSize()const{return mEdgeSize;}
bool BaseGraphClass::IsEmpty()const{if (VertexSize()==0 && EdgeSize()==0) return true; else return false;}

//---------------------------------------------------------------------------------









