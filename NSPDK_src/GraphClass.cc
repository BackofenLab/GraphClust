//* -*- mode:c++ -*- */
#include "GraphClass.h"

std::string& getcolor(const string& s) {
  static std::map<std::string,int> colorize;
  static std::string colors[] = {
    "#ffaa66", "#aaff66", "#ff66aa", "#66ffaa", "#aa66ff", "#66aaff",
    "#ff6666", "#66ff66", "#6666ff", "#ffff66", "#ff66ff", "#66ffff",
    "#ffaaaa", "#aaffaa", "#aaaaff", "#ffffaa", "#ffaaff", "#aaffff",
    "#ffdddd", "#ddffdd", "#ddddff", "#ffffdd", "#ffddff", "#ddffff",
    "#aaaa66", "#ff9900"
  };
  static int global_idx = 0;
  int idx;
  std::map<std::string,int>::const_iterator found = colorize.find(s);
  if (found != colorize.end()) {
    idx = found->second;
  } else {
    idx = global_idx = ((global_idx >= 25) ? 0 : (global_idx +1));
    colorize[s] = idx;
  }
  return colors[idx];
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------

ostream& operator<<(ostream& out, const GraphClass& aG){return aG.Output(out);}

string GraphClass::GetGraphID()const{return mGraphID;}
void GraphClass::SetGraphID(string aGraphID){mGraphID=aGraphID;}

ostream& GraphClass::Output(ostream& out)const{
  for (unsigned i=0;i<mAdjacencyList.size();++i){
    out<<"id: "<<i<<" ";
    out << "SymID: " << GetVertexSymbolicID(i) << std::endl;
    out<<"  symb[ ";
    for (unsigned k=0;k<mVertexSymbolicAttributeList[i].size();++k)
      out<<mVertexSymbolicAttributeList[i][k]<<" ";
    out<<"]" << std::endl;
    out<<"  num[";
    for (unsigned k=0;k<mVertexNumericAttributeList[i].size();++k)
      out<<mVertexNumericAttributeList[i][k]<<" ";
    out<<"]" << std::endl;
    out<<"  status[ ";
    out<<"kernel point:"<<(GetVertexKernelPoint(i)?"yes":"no")<<" ";
    out<<"vertex kind:"<<(GetVertexKind(i)?"entity":"relation")<<" ";
    out<<"first endpoint:"<<(GetVertexViewPoint(i)?"yes":"no")<<" ";
    out<<"dead:"<<(GetVertexDead(i)?"yes":"no")<<" ";
    out<<"]" << std::endl;
    out<<"  Adjacency list: ";
    for (unsigned k=0;k<mAdjacencyList[i].size();++k)
      out<<mAdjacencyList[i][k].mDestVertexID<<" ";
    out<<endl;
  }
  
  out<<"Map (src,dest) -> distance ("<<mSrcDestMaptoDistance.size()<<"):"<<endl;
  for(map<pair<unsigned,unsigned>,int>::const_iterator it=mSrcDestMaptoDistance.begin();it!=mSrcDestMaptoDistance.end();++it)
    out<<"(src id:"<<it->first.first<<" , dest id:"<<it->first.second<<") -> distance:"<<it->second<<endl;
  
  out<<"Map (src,distance) -> dest list ("<<mSrcDistanceMaptoDestList.size()<<"):"<<endl;
  for(  map<pair<unsigned,int>,vector<unsigned> >::const_iterator it=mSrcDistanceMaptoDestList.begin();it!=mSrcDistanceMaptoDestList.end();++it){
    out<<"(src id:"<<it->first.first<<" , distance:"<<it->first.second<<") -> dest id list: ";
    const vector<unsigned>& dest_list=it->second;
    for (unsigned i=0;i<dest_list.size();++i)	out<<dest_list[i]<<" ";
    out<<endl;
  }
  return out;
}

// void GraphClass::SaveAsUnionDotFile(const string& aGid, const string& aFilename)const{
//   std::ofstream dot_stream;
//   dot_stream.open(aFilename.c_str(),ios::out);
//   dot_stream << "graph \"" << aGid << "\"{" << std::endl;
//   for (unsigned v=0; v<mVertexSize; ++v) {
//     dot_stream << v << " [label=\"";
//     for (unsigned j=0;j<mVertexSymbolicAttributeList[v].size();++j) {
//       if (j==0)
// 	dot_stream << mVertexSymbolicAttributeList[v][j] << (mVertexSymbolicAttributeList[v].size()>1 ? "(" : "");
//       else if (j<mVertexSymbolicAttributeList[v].size()-1)
// 	dot_stream 
// 	  //<< j << ":" 
// 	  << mVertexSymbolicAttributeList[v][j] << ",";
//       else
// 	dot_stream 
// 	  //<< j << ":" 
// 	  << mVertexSymbolicAttributeList[v][j] << ")";
//     }
//     for (unsigned j=0;j<mVertexNumericAttributeList[v].size();++j) {
//       if (j==0) dot_stream << "(";
//       dot_stream 
// 	//<< j << ":" 
// 	<< mVertexNumericAttributeList[v][j]
// 	<< (j==mVertexNumericAttributeList[v].size()-1 ? ")" :",");
//     }
//     dot_stream << "\", shape=\"circle\""
// 	       << ",width=\""<< stream_cast<string>(mVertexNumericAttributeList[v][0])<<"\""
// 	       << ",color=\"" << getcolor(mVertexSymbolicAttributeList[v][0])<<"\""
// 	       << ",style=\"filled\""
// 	       <<" , penwidth=\""<<stream_cast<string>(mVertexNumericAttributeList[v][0]*3)<<"\""
//                << "]" << std::endl;
//   }
//   for (unsigned u=0; u<mVertexSize; ++u) {
//     for (unsigned vpos=0; vpos < mAdjacencyList[u].size(); ++vpos) {
//       unsigned int v = mAdjacencyList[u][vpos].mDestVertexID;
//       if (u<v){ // internally there are both u-v and v-u but we want to draw an undirected graph!
// 	unsigned edge_id = mAdjacencyList[u][vpos].mEdgeID;
// 	dot_stream << u << " -- " << v << " [label=\"";
// 	for (unsigned j=0;j<mEdgeSymbolicAttributeList[edge_id].size();++j){
// 	  dot_stream
// 	    //<< j << ":"
// 	    << mEdgeSymbolicAttributeList[edge_id][j]<<" ";
// 	}
// 	for (unsigned j=0;j<mEdgeNumericAttributeList[edge_id].size();++j){
// 	  if (j==0) dot_stream << "(";
// 	  dot_stream
// 	    //<< j << ":" 
// 	    << mEdgeNumericAttributeList[edge_id][j]
// 	    << (j==mEdgeNumericAttributeList[edge_id].size()-1 ? ")" :",");
// 	}
// 	dot_stream << "\"";
// 	dot_stream<<" , weight=\""<<mEdgeNumericAttributeList[edge_id][0]<<"\"";
// 	dot_stream<<" , penwidth=\""<<mEdgeNumericAttributeList[edge_id][0]<<"\"";
// 	dot_stream << "]" << std::endl;
//       } else {}
//     }
//   }
//   dot_stream << "}" << std::endl;
//   dot_stream.close();
// }






std::string GraphClass::vertex_label_serialize(unsigned v, std::string separator) const {
  stringstream os;
  os << (GetVertexKernelPoint(v) == true ? "*" : "");
  for (unsigned j=0;j<mVertexSymbolicAttributeList[v].size();++j) {
    if (j==0)
      os << mVertexSymbolicAttributeList[v][j] << (mVertexSymbolicAttributeList[v].size()>1 ? "(" : "");
    else if (j<mVertexSymbolicAttributeList[v].size()-1)
      os 
        //<< j << ":" 
        << mVertexSymbolicAttributeList[v][j] << ",";
    else
      os 
        //<< j << ":" 
        << mVertexSymbolicAttributeList[v][j] << ")";
  }
  for (unsigned j=0;j<mVertexNumericAttributeList[v].size();++j) {
    if (j==0) os << "(";
    os 
      //<< j << ":" 
      << mVertexNumericAttributeList[v][j]
      << (j==mVertexNumericAttributeList[v].size()-1 ? ")" :",");
  }
  if (GetVertexKind(v) == true) os << separator << GetVertexSymbolicID(v);
  if (IsSliced()) os << separator << "{" << GetSliceID(v) << "}";
  return os.str();
}

std::string GraphClass::edge_label_serialize(unsigned edge_id) const {
  stringstream os;
  for (unsigned j=0;j<mEdgeSymbolicAttributeList[edge_id].size();++j)
    os
      //<< j << ":"
      << mEdgeSymbolicAttributeList[edge_id][j] << " ";
  for (unsigned j=0;j<mEdgeNumericAttributeList[edge_id].size();++j)
    os
      //<< j << ":" 
      << mEdgeNumericAttributeList[edge_id][j] << " ";
  return os.str();
}



 
/**
   Create a dot file of a given interpretation for
   debugging/visualization. E-vertices are represented as boxes,
   R-vertices as diamonds, dead vertices are dashed. Symbolic IDs and
   slice identifiers if meaningful are printed for debugging purposes.
 */
void GraphClass::SaveAsDotFile(const string& aFilename)const{
  std::ofstream dot_stream;
  dot_stream.open(aFilename.c_str(),ios::out);
  dot_stream << "graph \"" << mGraphID << "\"{" << std::endl;
  for (unsigned v=0; v<mVertexSize; ++v) {
    dot_stream << v << " [label=\"" 
               << vertex_label_serialize(v) << "\", shape="
               << (GetVertexKind(v) == true ? "box" : "diamond")
      	       << ",style=" << (GetVertexDead(v) == false ? "filled" : "dashed")
               << ",fillcolor=\"" << getcolor(IsSliced()? GetSliceID(v) : mVertexSymbolicAttributeList[v][0])
               << "\"]" << std::endl;
  }
  for (unsigned u=0; u<mVertexSize; ++u) {
    for (unsigned vpos=0; vpos < mAdjacencyList[u].size(); ++vpos) {
      unsigned int v = mAdjacencyList[u][vpos].mDestVertexID;
      if (v<u){ // internally there are both u-v and v-u but we want to draw an undirected graph!
	unsigned edge_id = mAdjacencyList[u][vpos].mEdgeID;
	dot_stream << u << " -- " 
                   << v << " [label=\"" 
                   << edge_label_serialize(edge_id) << "\"]" << std::endl;
      } else {}
    }
  }
  dot_stream << "}" << std::endl;
  dot_stream.close();
}


/**
   Create a GML file of a given interpretation for
   debugging/visualization. E-vertices are represented as boxes,
   R-vertices as diamonds, dead vertices are dashed. Symbolic IDs and
   slice identifiers if meaningful are printed for debugging purposes.
 */
void GraphClass::SaveAsGMLFile(const string& aFilename)const{
  std::ofstream gml_stream;
  gml_stream.open(aFilename.c_str(),ios::out);
  gml_stream << "graph [" << std::endl
             << "\tcomment \"" << mGraphID << "\"" << std::endl
             << "\tdirected 0"  << std::endl
             << "\tid 0"  << std::endl
             << "\tlabel \"" << mGraphID << "\"" << std::endl;
  for (unsigned v=0; v<mVertexSize; ++v) {
    gml_stream << "\tnode ["  << std::endl
               << "\t\tid " << v << std::endl
               << "\t\tLabelGraphics [ text \""  << vertex_label_serialize(v," - ") << "\" ]" << std::endl;
    gml_stream << "\t\tgraphics [" << std::endl
               << "\t\t\ttype \""<<(GetVertexKind(v) == true ? "rectangle" : "diamond") << "\"" << std::endl
      	       << "\t\t\tfill \""<<(getcolor(IsSliced()?GetSliceID(v):mVertexSymbolicAttributeList[v][0]))<<"\""<<std::endl
               << "\t\t]" << std::endl;
    gml_stream << "\t]"  << std::endl; // end node
  }
  for (unsigned u=0; u<mVertexSize; ++u) {
    for (unsigned vpos=0; vpos < mAdjacencyList[u].size(); ++vpos) {
      unsigned int v = mAdjacencyList[u][vpos].mDestVertexID;
      if (v<u){ // internally there are both u-v and v-u but we want to draw an undirected graph!
	unsigned edge_id = mAdjacencyList[u][vpos].mEdgeID;
	gml_stream << "\tedge [" << std::endl
                   << "\t\tlabel \"" << edge_label_serialize(edge_id) << "\"" << std::endl;
        gml_stream << "\t\tsource " << u << std::endl
                   << "\t\ttarget " << v << std::endl;
        gml_stream << "\t]"  << std::endl; // end edge
      } else {}
    }
  }
  gml_stream << "]" << std::endl;
  gml_stream.close();
}

/**
   Create a GDL file of a given interpretation for
   debugging/visualization (using aiSee). E-vertices are represented
   as boxes, R-vertices as diamonds, dead vertices are
   dashed. Symbolic IDs and slice identifiers if meaningful are
   printed for debugging purposes.
 */
void GraphClass::SaveAsGDLFile(const string& aFilename)const{
  std::ofstream gdl_stream;
  gdl_stream.open(aFilename.c_str(),ios::out);
  gdl_stream << "graph: { title: " << "\""<< mGraphID << "\"" << std::endl
             << "\tlayoutalgorithm: forcedir" << std::endl
             << std::endl;
  for (unsigned v=0; v<mVertexSize; ++v) {
    gdl_stream << "node: {"  
      //               << " color: \"" << (getcolor(IsSliced()?GetSliceID(v):mVertexSymbolicAttributeList[v][0]))<<"\""
               << " title: \"" << v << "\"" 
               << " shape: " <<(GetVertexKind(v) == true ? "box" : "rhomb") 
               << " label: \""  << vertex_label_serialize(v,"\n") << "\" }" << std::endl;
  }
  for (unsigned u=0; u<mVertexSize; ++u) {
    for (unsigned vpos=0; vpos < mAdjacencyList[u].size(); ++vpos) {
      unsigned int v = mAdjacencyList[u][vpos].mDestVertexID;
      if (v<u){ // internally there are both u-v and v-u but we want to draw an undirected graph!
	//unsigned edge_id = mAdjacencyList[u][vpos].mEdgeID;
	gdl_stream << "edge: { " 
                   << "sourcename: \"" << u << "\""
                   << " targetname: \"" << v << "\" }" << std::endl;
      } else {}
    }
  }
  gdl_stream << "}" << std::endl;
  gdl_stream.close();
}



/**
   Create a CSV file of a given interpretation for
   debugging/visualization. Only interactions (between entities and
   relationships) are saved.
 */
void GraphClass::SaveAsCSVFile(const string& aFilename)const{
  std::ofstream csv_stream;
  csv_stream.open(aFilename.c_str(),ios::out);
  for (unsigned u=0; u<mVertexSize; ++u) {
    for (unsigned vpos=0; vpos < mAdjacencyList[u].size(); ++vpos) {
      unsigned int v = mAdjacencyList[u][vpos].mDestVertexID;
      if (v<u){ // internally there are both u-v and v-u but we want to draw an undirected graph!
        unsigned edge_id = mAdjacencyList[u][vpos].mEdgeID;
        csv_stream << u << "," << vertex_label_serialize(u, " - ") << ","
                   << v << "," << vertex_label_serialize(v, " - ") << ","
                   << edge_id << "," << edge_label_serialize(edge_id) << std::endl;
      } else {}
    }
  }
  csv_stream.close();
}


void GraphClass::SaveAsGspanFile(const string& aFilename)const{
  std::ofstream gspan_stream;
  gspan_stream.open(aFilename.c_str(),ios::out);
  gspan_stream << "t # " << mGraphID << std::endl;
  for (unsigned v=0; v<mVertexSize; ++v) {
    gspan_stream << "v " << v << " ";
    for (unsigned i=0;i<mVertexSymbolicAttributeList[v].size();i++)
      gspan_stream << mVertexSymbolicAttributeList[v][i]<<" ";
    for (unsigned i=0;i<mVertexNumericAttributeList[v].size();i++)
      gspan_stream << mVertexNumericAttributeList[v][i]<<" ";
    gspan_stream << endl;
  }
  for (unsigned u=0; u<mVertexSize; ++u) {
    for (unsigned vpos=0; vpos < mAdjacencyList[u].size(); ++vpos) {
      unsigned int v = mAdjacencyList[u][vpos].mDestVertexID;
      unsigned e = mAdjacencyList[u][vpos].mEdgeID;
      if (u<v){ // internally there are both u-v and v-u but we want to draw an undirected graph!
	gspan_stream << "e " << u << " " << v << " ";
	for (unsigned i=0;i<mEdgeSymbolicAttributeList[e].size();i++)
	  gspan_stream << mEdgeSymbolicAttributeList[e][i] << " ";
	for (unsigned i=0;i<mEdgeNumericAttributeList[e].size();i++)
	  gspan_stream << mEdgeNumericAttributeList[e][i] << " ";
	gspan_stream << endl;
      }
      else {}
    }
  }
  gspan_stream.close();
}

void GraphClass::ExportGraph(const string& aFilename, const string& aFormat) const {
  if (aFormat == "dot")
    SaveAsDotFile(aFilename + "." + aFormat);
  else if (aFormat == "gml")
    SaveAsGMLFile(aFilename + "." + aFormat);
  else if (aFormat == "gdl")
    SaveAsGDLFile(aFilename + "." + aFormat);
  else if (aFormat == "csv")
    SaveAsCSVFile(aFilename + "." + aFormat);
  else if (aFormat == "gspan")
    SaveAsGspanFile(aFilename + "." + aFormat);
  else
    throw std::out_of_range("Invalid export format in GraphClass::ExportGraph()");
}


void GraphClass::SetSliceID(unsigned aID, const string& aSliceID) {
  if (aID>=mSliceIdList.size())
    mSliceIdList.resize(2*mSliceIdList.size());
  mSliceIdList[aID] = aSliceID;
}

string GraphClass::GetSliceID(unsigned aID) const {
  return mSliceIdList[aID];
}

void GraphClass::SetVertexKernelPoint(unsigned aID, bool aState){
  mTopologicalChangeOccurrence=true;
  SetVertexStatusAttributeList(aID,KERNEL_POINT_ID,aState);
}
bool GraphClass::GetVertexKernelPoint(unsigned aID) const {
  return GetVertexStatusAttributeList(aID,KERNEL_POINT_ID);
  }
void GraphClass::SetVertexKind(unsigned aID, bool aKind){
  SetVertexStatusAttributeList(aID,KIND_ID,aKind);
  }
bool GraphClass::GetVertexKind(unsigned aID) const {
  return GetVertexStatusAttributeList(aID,KIND_ID);
  }
void GraphClass::SetVertexViewPoint(unsigned aID, bool aState){
  SetVertexStatusAttributeList(aID,VIEWPOINT_ID,aState);
  }
bool GraphClass::GetVertexViewPoint(unsigned aID) const {
  return GetVertexStatusAttributeList(aID,VIEWPOINT_ID);
  }
void GraphClass::SetVertexAlive(unsigned aID, bool aState){
  SetVertexDead(aID,! aState);
}
bool GraphClass::GetVertexAlive(unsigned aID) const {
  return ! GetVertexDead(aID);
  }
void GraphClass::SetVertexDead(unsigned aID, bool aState){
  mTopologicalChangeOccurrence=true;
  SetVertexStatusAttributeList(aID,DEAD_ID,aState);
}
void GraphClass::KillVertices(std::string aLabel){
  mTopologicalChangeOccurrence=true;
  for (unsigned v = 0; v < mVertexSize; ++v) {
    if (GetVertexLabel(v)==aLabel) {
      SetVertexDead(v,true);
      cout << "Graph_id:" << mGraphID << " killed_vertex_id:" << v <<" label:"<< vertex_label_serialize(v," - ") << std::endl;////FIXME: only for debugging purposes
    }
  }
}

bool GraphClass::GetVertexDead(unsigned aID) const {
  return GetVertexStatusAttributeList(aID,DEAD_ID);
  }

bool GraphClass::GetVertexAbstraction(unsigned aID) const {
  return GetVertexStatusAttributeList(aID,ABSTRACTION_ID);
  }
void GraphClass::SetVertexAbstraction(unsigned aID, bool aStatus){
  SetVertexStatusAttributeList(aID,ABSTRACTION_ID,aStatus);
  }

void GraphClass::SetVertexLabel(unsigned aID, const string& aLabel){
  SetVertexSymbolicAttributeList(aID,LABEL_VERTEX_ATTRIBUTE_ID,aLabel);
  }
string GraphClass::GetVertexLabel(unsigned aID)const{
  return GetVertexSymbolicAttributeList(aID,LABEL_VERTEX_ATTRIBUTE_ID);
}
string GraphClass::GetEdgeLabel(unsigned aSrcID, unsigned aDestID)const{
  return GetEdgeSymbolicAttributeList(aSrcID,aDestID,EDGE_ATTRIBUTE_ID);
  }
void GraphClass::SetEdgeLabel(unsigned aSrcID, unsigned aDestID, const string& aLabel){
  SetEdgeSymbolicAttributeList(aSrcID,aDestID,EDGE_ATTRIBUTE_ID,aLabel);
}
bool GraphClass::GetEdgeAbstractionOf(unsigned aSrcID, unsigned aDestID)const{
  return GetEdgeStatusAttributeList(aSrcID,aDestID,EDGE_ABSTRACTIONOF_ID);
  }
void GraphClass::SetEdgeAbstractionOf(unsigned aSrcID, unsigned aDestID, bool aStatus){
  SetEdgeStatusAttributeList(aSrcID,aDestID,EDGE_ABSTRACTIONOF_ID,aStatus);
}
bool GraphClass::GetEdgePartOf(unsigned aSrcID, unsigned aDestID)const{
  return GetEdgeStatusAttributeList(aSrcID,aDestID,EDGE_PARTOF_ID);
  }
void GraphClass::SetEdgePartOf(unsigned aSrcID, unsigned aDestID, bool aStatus){
  SetEdgeStatusAttributeList(aSrcID,aDestID,EDGE_PARTOF_ID,aStatus);
}
bool GraphClass::Check()const{
  if (VertexSize()==0) throw logic_error("Empty data structure: vertex set");
  //if (EdgeSize()==0) throw logic_error("Empty data structure: edge set");
  if (mSrcDestMaptoDistance.size()==0) throw logic_error("Empty data structure: (src,dest) -> distance map");
  if (mSrcDistanceMaptoDestList.size()==0) throw logic_error("Empty data structure: (src,distance) -> dest list map");
  return true;
}

/**
   Returns a vector of ids of vertices adjacent to aID that are
   ALIVE. Overriden function from
   BaseGraphClass::GetVertexAdjacentList(unsigned aID)const.
 */
vector<unsigned> GraphClass::GetVertexAdjacentList(unsigned aID)const{
  vector<unsigned> adjacent_list;
  for (unsigned i=0;i<mAdjacencyList[aID].size();++i){
    unsigned vid=mAdjacencyList[aID][i].mDestVertexID;
    if (GetVertexAlive(vid)) 
      adjacent_list.push_back(vid);
  }
  return adjacent_list;
}

void GraphClass::ComputePairwiseDistanceInformation(int aMaxDistance,int aMaxRadius, vector<unsigned> aViewPointList)const{
  if (aMaxDistance!=-1) mMaxDistance=aMaxDistance;
  if (aMaxRadius!=-1) mMaxRadius=aMaxRadius;
  if (mTopologicalChangeOccurrence==true) {
    unsigned distance_bound=max(mMaxRadius,mMaxDistance);
    mSrcDestMaptoDistance.clear();
    if (aViewPointList.size()==0){
      for (unsigned i=0;i<VertexSize();i++){
	if (GetVertexDead(i)==false) // Paolo
	  if (GetVertexAbstraction(i)==false)
	    SingleVertexBoundedBreadthFirstVisit(i,distance_bound,mSrcDestMaptoDistance);
      } 
    } else {
      for (unsigned i=0;i<aViewPointList.size();i++){
	unsigned id=aViewPointList[i];
	if (GetVertexDead(id)==false) // Paolo
	  if (GetVertexAbstraction(i)==false)
	    SingleVertexBoundedBreadthFirstVisit(id,distance_bound,mSrcDestMaptoDistance);
      }
    }
    mSrcDistanceMaptoDestList.clear();
    for (map<pair<unsigned,unsigned>,int>::iterator it=mSrcDestMaptoDistance.begin();it!=mSrcDestMaptoDistance.end();++it){
      unsigned src_id=it->first.first;
      unsigned dest_id=it->first.second;
      int distance=it->second;
      mSrcDistanceMaptoDestList[make_pair(src_id,distance)].push_back(dest_id);
    }

    mTopologicalChangeOccurrence=false;
  } else {}
}

void GraphClass::SingleVertexBoundedBreadthFirstVisit(unsigned aRootVertexIndex, int aRadius, map<pair<unsigned,unsigned>,int>& oSrcDestMaptoDistance)const{
  map<int,int> dest_mapto_distance;
  dest_mapto_distance[aRootVertexIndex]=0;
  map<int,bool> already_explored;
  already_explored[aRootVertexIndex]=true;
  queue<int> q;
  q.push(aRootVertexIndex);//initialize queue with the root vertex
  while (q.empty()==false){
    int u=q.front();
    for (unsigned j=0; j<mAdjacencyList[u].size();j++){
      int v=mAdjacencyList[u][j].mDestVertexID;
      if (already_explored[v]==true){}
      else if (GetVertexDead(v)==true){}
      else if (GetVertexAbstraction(v)==true){}
      else {
	if (dest_mapto_distance[u]+1<=aRadius){
	  dest_mapto_distance[v]=dest_mapto_distance[u]+1;
	  already_explored[v]=true;
	  q.push(v);
	}
      }
    }
    q.pop();
  }
  //compute (src,dest) \mapto distance
  for (map<int,int>::const_iterator it=dest_mapto_distance.begin();it!=dest_mapto_distance.end();++it){
    unsigned dest_vertex_id=(unsigned)(it->first);
    int distance=it->second;
    oSrcDestMaptoDistance.insert(make_pair(make_pair(aRootVertexIndex,dest_vertex_id),distance));
  }
}

vector<unsigned> GraphClass::GetFixedDistanceVertexIDList(unsigned aSrcID, int aDistance)const{
  if (mTopologicalChangeOccurrence==true) ComputePairwiseDistanceInformation(); 
  if (mSrcDistanceMaptoDestList.count(make_pair(aSrcID,aDistance))==0) return vector<unsigned>(0);
  else return mSrcDistanceMaptoDestList.find(make_pair(aSrcID,aDistance))->second;
}
int GraphClass::PairwiseDistance(unsigned aSrcID, unsigned aDestID)const{
  if (mTopologicalChangeOccurrence==true) ComputePairwiseDistanceInformation();
  if (mSrcDestMaptoDistance.count(make_pair(aSrcID,aDestID))==0) return -1;
  else return (mSrcDestMaptoDistance.find(make_pair(aSrcID,aDestID)))->second;
}

vector<unsigned> GraphClass::GetUnionShortestPathsVertexIDList(unsigned aSrcID, unsigned aDestID, unsigned aDistance)const{
  set<unsigned> union_shortest_paths_set;
  for (unsigned d=0;d<=aDistance;d++){
    //get all vertices at distance d from src 
    vector<unsigned> src_neighbors=GetFixedDistanceVertexIDList(aSrcID, d);
    set<unsigned> src_neighbors_set;
    for (unsigned i=0;i<src_neighbors.size();i++) src_neighbors_set.insert(src_neighbors[i]);
    //get all vertices at distance aDistance-d from dest
    vector<unsigned> dest_neighbors=GetFixedDistanceVertexIDList(aDestID, aDistance-d);
    set<unsigned> dest_neighbors_set;
    for (unsigned i=0;i<dest_neighbors.size();i++) dest_neighbors_set.insert(dest_neighbors[i]);
    //store intersection set
    set_intersection( src_neighbors_set.begin(), src_neighbors_set.end(), dest_neighbors_set.begin(), dest_neighbors_set.end(), inserter(union_shortest_paths_set,union_shortest_paths_set.begin()) );	
  } 
  vector<unsigned> union_shortest_paths_list;
  for (set<unsigned>::iterator it=union_shortest_paths_set.begin();it!=union_shortest_paths_set.end();++it)
    union_shortest_paths_list.push_back(*it);
  return union_shortest_paths_list;
}

set<unsigned> GraphClass::GetUnionThickShortestPathsVertexIDSet(unsigned aSrcID, unsigned aDestID, unsigned aDistance, unsigned aThickness)const{
  vector<unsigned> path=GetUnionShortestPathsVertexIDList(aSrcID, aDestID, aDistance);
  set<unsigned> thick_path;
  for (unsigned i=0;i<path.size();i++){
    //for each vertex extract neighborhood subgraph
    vector<unsigned> neighborhood=GetNeighborhoodVertexIDList(path[i],aThickness);
    for (unsigned i=0;i<neighborhood.size();++i)
      thick_path.insert(neighborhood[i]);
  }
  return thick_path;
}

vector<unsigned> GraphClass::GetNeighborhoodVertexIDList(unsigned aSrcID, unsigned aRadius)const{
  vector<unsigned> neighborhood;
  for (unsigned d=0;d<=aRadius;d++){
    vector<unsigned> circle=GetFixedDistanceVertexIDList(aSrcID,d);
    for (unsigned i=0;i<circle.size();++i)
      neighborhood.push_back(circle[i]);
  }
  return neighborhood;
}
