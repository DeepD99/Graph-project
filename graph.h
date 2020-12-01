/*graph.h*/

//
// <<Deep Diganvker>>
//
// Basic graph class using adjacency matrix representation.  Currently
// limited to a graph with at most 100 vertices.
//
// original author: Prof. Joe Hummel
// U. of Illinois, Chicago
// CS 251: Spring 2020
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <vector>
#include <set>
#include <map>

using namespace std;


template<typename VertexT, typename WeightT>
class graph
{
private:
  struct adjList
  {
    bool     EdgeExists;
	VertexT  From;
	VertexT  To;
    WeightT  Weight;
	
    adjList()
    {
      EdgeExists = false;
    }
  };
  
  int numEdge;
  map <VertexT, vector <adjList>> Data;
  vector<VertexT>  Vertices;
  vector<adjList> edges;
  adjList adjInst;

  
  
  bool _LookupVertex(VertexT v) const
  {
    if(Data.find(v) != Data.end()){
		return true;
	}
    return false;
  }


public:
  
  
  //default constructor
  graph(){
	  numEdge = 0;
  }
  
  
  graph(const graph& other){
	  this->numEdge = other.numEdge;
	  this->Vertices = other.Vertices;
	  this->Data = other.Data;
  }
  
  
  graph& operator=(const graph& other){
	  this->numEdge = other.numEdge;
	  this->Vertices = other.Vertices;
	  this->Data = other.Data;
	  return *this;
  }
	

  // NumVertices
  int NumVertices() const
  {
	return static_cast<int>(this->Vertices.size());
  }

  
  //NumEdges
  int NumEdges() const
  {
    return numEdge;
  }

  
  //addVertex
  bool addVertex(VertexT v)
  {
  
	if(_LookupVertex(v)){
		return false;
	}
	else{
		this->Vertices.push_back(v);
		Data.insert(make_pair(v, edges));
	}
    return true;
  }


  //addEdge
  bool addEdge(VertexT from, VertexT to, WeightT weight)
  {
    if(_LookupVertex(from) == false){
		return false;
	}
	
	if(_LookupVertex(to) == false){
		return false;
	}
	
	else{
		for(unsigned int i = 0; i < Data.at(from).size(); i++){
			if(Data.at(from).at(i).To == to){
				Data.at(from).at(i).Weight = weight;
				return true;
			}
		}
		
		adjInst.From = from;
		adjInst.To = to;
		adjInst.Weight = weight;
		Data.at(from).push_back(adjInst);
		numEdge++;
	}
    return true;
  }


  //getWeight
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const
  {
    if(_LookupVertex(from) == false){
		return false;
	}
	
	if(_LookupVertex(to) == false){
		return false;
	}
	
	else{
		for(unsigned int i = 0; i < Data.at(from).size(); i++){
			if(to == Data.at(from).at(i).To){
				weight = Data.at(from).at(i).Weight;
				return true;
			}
		}
	}
    return false;
  }


  //neighbors
  set<VertexT> neighbors(VertexT v) const
  {
    set<VertexT>  S;
	
	if(!_LookupVertex(v)){
		return S;
	}
	
	else{
		for(unsigned int i = 0; i < Data.at(v).size(); i++){
			VertexT dest = Data.at(v).at(i).To;
			S.insert(dest);
		}
	}
    return S;
  }

  
  //getVertices
  vector<VertexT> getVertices() const
  {
    return this->Vertices;  // returns a copy:
  }
  

  //
  // dump
  // 
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const
  {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:" << endl;
	
	vector <adjList> vec;
	
	for(auto graph: Data){
		output << graph.first << ": ";
		vec = this->Data.find(graph.first)->second;
		for(adjList adj: vec){
			output << "(" << graph.first << ", " << adj.To << ", " << adj.Weight << ")";
		}
		output << endl;
	}
    output << "**************************************************" << endl;
  }

};






