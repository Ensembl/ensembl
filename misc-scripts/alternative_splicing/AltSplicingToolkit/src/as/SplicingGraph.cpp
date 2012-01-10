/*
 *    AltSplicingToolkit
 *    Author: Gautier Koscielny <koscieln@ebi.ac.uk>
 *
 *    Copyright (c) 1999-2012 The European Bioinformatics Institute and 
 *    Genome Research Limited, and others.  All rights reserved. 
 *
 *    Redistribution and use in source and binary forms, with or without 
 *    modification, are permitted provided that the following conditions 
 *    are met: 
 *
 *    1. Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer. 
 *
 *    2. Redistributions in binary form must reproduce the above copyright 
 *       notice, this list of conditions and the following disclaimer 
 *       in the documentation and/or other materials provided with the 
 *       distribution. 
 *
 *    3. The name "Ensembl" must not be used to endorse or promote products 
 *       derived from this software without prior written permission.  For 
 *       written permission, please contact helpdesk@ensembl.org 
 *
 *    4. Products derived from this software may not be called "Ensembl" 
 *       nor may "Ensembl" appear in their names without prior written 
 *       permission of the Ensembl developers. 
 *
 *    5. Redistributions in any form whatsoever must retain the following 
 *       acknowledgement: 
 *
 *         "This product includes software developed by Ensembl 
 *         (http://www.ensembl.org/)" 
 *
 *    THIS SOFTWARE IS PROVIDED BY THE ENSEMBL GROUP ``AS IS'' AND ANY 
 *    EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
 *    PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE ENSEMBL GROUP OR ITS 
 *    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
 *    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
 *    OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 *    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 *    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 *    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 */

#include "SplicingGraph.h"
#include "SplicingVertex.h"
#include "Exon.h"
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <cassert>

namespace as {

  int SplicingGraph::objectCount = 0;

  SplicingGraph::SplicingGraph()
  {
    objectCount++;
    size = 0;
    printMsg("SplicingGraph()");
  }

  SplicingGraph::SplicingGraph(std::vector<int> genomicLocations,
			       SplicingVertexMap vMap,
			       Edges e) : genomicLocations(genomicLocations), vertexMap(vMap), edges(e)
  {
    objectCount++;

    this->size = (unsigned int) this->vertexMap.size();

    std::cout << " Number of vertices " << vertexMap.size() << std::endl;
    std::cout << " Number of edges " << edges.size() << std::endl;
    std::cout << " Create the adjacency matrix " << std::endl;

    // resize the adjacency matrix
    create_matrix(&adjacencyMatrix, vertexMap.size());

    std::cout << " adjacency matrix done " << std::endl;

    //
    // loop over the edges and fill the adjacency matrix
    //

    std::vector<Edge>::const_iterator edgeIt;

    for(edgeIt=edges.begin(); edgeIt!=edges.end(); edgeIt++)
      {
	//std::cout << edgeIt->first << "\t" << edgeIt->second << std::endl;
	// get the index for the first
	int firstIndex = vertexMap[edgeIt->first].getIndex();
	int secondIndex = vertexMap[edgeIt->second].getIndex();

	adjacencyMatrix[firstIndex][secondIndex] = true;

      }

    //display_matrix(&adjacencyMatrix, vertexMap.size());

    // compute the input and outputs
    computeInputsOutputs();

    std::cout << "I/O\t" << inputs.size() << "/" << outputs.size() << std::endl;

    //adjacencyMatrix = matrix;

    std::cout << " Graph representation: " << std::endl;

    //
    // loop over the vertexMap
    //

    for( SplicingVertexMap::iterator ii=vertexMap.begin(); ii!=vertexMap.end(); ++ii)
      {
	int index = (*ii).second.getIndex();
	std::cout << (*ii).first << ": " <<  index;

	// for each of them, look at the in and out_degree
	std::cout << "\tin: " << this->getInDegree(index);
	std::cout << "\tout: " << this->getOutDegree(index) << std::endl;

      }



  }

  SplicingGraph::~SplicingGraph()
  {

    objectCount--;
    delete_matrix(&adjacencyMatrix, size);

    printMsg("~SplicingGraph()");
  }

  void SplicingGraph::create_matrix(bool ***m, int size)
  {

    (*m) = new bool*[size];

    for (int i = 0 ; i < size; i++)
    {
      (*m)[i] = new bool[size];

      for (int j = 0; j < size; j++)
	(*m)[i][j] = false;

    }
  }

  void SplicingGraph::delete_matrix(bool ***m, int size)
  {
    for (int i = 0 ; i < size; i++)
      {
	delete [] (*m)[i];
      }

    delete [] (*m);
  }

  void SplicingGraph::display_matrix(bool ***m, int size)
  {
    for (int i = 0 ; i < size; i++)
      {
	for (int j = 0; j < size; j++)
	  std::cout << (*m)[i][j] << " ";
	std::cout << std::endl;
      }
  }

  void SplicingGraph::printMsg(const std::string& msg)
  {
    if (msg.size() != 0) std::cout << msg << ": ";
    std::cout << "objectCount = "
	      << objectCount << std::endl;
  }

  void SplicingGraph::computeInputsOutputs()
  {
    //
    // loop over the vertexMap
    //

    for( SplicingVertexMap::iterator ii=vertexMap.begin(); ii!=vertexMap.end(); ++ii)
      {
	int index = (*ii).second.getIndex();

	if (this->getInDegree(index) == 0)
	  {
	    inputs.push_back((*ii).second);
	  }

	if (this->getOutDegree(index) == 0)
	  {
	    outputs.push_back((*ii).second);
	  }
      }
  }

  const SplicingVertexVector SplicingGraph::intersect(SplicingVertexVector setA, SplicingVertexVector setB)
  {

    SplicingVertexVector inter;

    SplicingVertexVector::const_iterator cItA;
    SplicingVertexVector::const_iterator cItB;

    for(cItA=setA.begin(); cItA!=setA.end(); cItA++)
      {

	for(cItB=setB.begin(); cItB!=setB.end(); cItB++)
	  {

	    if (cItA->getIndex() == cItB->getIndex())
	      inter.push_back(*cItA);

	  }
      }

    return inter;
  }

  SplicingVertexVector SplicingGraph::getInputs()
  {
    // get the inputs of the graph
    return inputs;
  }

  SplicingVertexVector SplicingGraph::getOutputs()
  {
    // get the outputs of the graph
    return outputs;
  }

  int SplicingGraph::getInDegree(unsigned int vertex)
  {
    int inDegree = 0;
    for (BooleanIndex i=0; i<this->size; i++)
      {
	if (adjacencyMatrix[i][vertex])
	  inDegree++;
      }
    return inDegree;
  }

  int SplicingGraph::getOutDegree(unsigned int vertex)
  {
    int outDegree = 0;
    for (BooleanIndex i=0; i<this->size; i++)
      {
	if (adjacencyMatrix[vertex][i])
	  outDegree++;
      }
    return outDegree;
  }

  SplicingVertexVector SplicingGraph::getInVertices(const SplicingVertex& vertex)
  {
    SplicingVertexVector v;
    unsigned int vertexIndex = vertex.getIndex();

    for (BooleanIndex i=0; i<this->size; i++)
      {
	if (adjacencyMatrix[i][vertexIndex])
	  {
	    int gLoc = genomicLocations[i];
	    v.push_back(vertexMap[gLoc]);
	  }
      }
    return v;

  }

  SplicingVertexVector SplicingGraph::getOutVertices(const SplicingVertex& vertex)
  {
    SplicingVertexVector v;
    unsigned int vertexIndex = vertex.getIndex();

    for (BooleanIndex i=0; i<this->size; i++)
      {
	if (adjacencyMatrix[vertexIndex][i])
	  {
	    int gLoc = genomicLocations[i];
	    v.push_back(vertexMap[gLoc]);
	  }
      }
    return v;
  }

  SplicingGraph SplicingGraph::buildGraph(const Transcript& t1, const Transcript& t2)
  {

    int vertexCount = 0;
    int lastDonorSite = 0;

    std::vector<int> genomicLocations;
    SplicingVertexMap vertexMap;
    Edges edges;

    //
    // list all exons from transcript t1
    //

    ExonVector exonsT1 = t1.getExons();

    //
    // put all t1 exons in the SplicingVertex Map
    // there is always an acceptor/donor site pairs
    // hence 2 vertices
    // add them to the vertexMap
    //

    // no donor site for the moment
    lastDonorSite = 0;

    std::vector<Exon>::const_iterator cii;
    for(cii=exonsT1.begin(); cii!=exonsT1.end(); cii++)
      {

	// acceptor site
	SplicingVertex acceptor(cii->getStart(), ACCEPTOR_SITE, vertexCount++);
	genomicLocations.push_back(cii->getStart());
	acceptor.addExon(*cii);

	// donor site
	SplicingVertex donor(cii->getEnd(), DONOR_SITE, vertexCount++);
	genomicLocations.push_back(cii->getEnd());
	donor.addExon(*cii);

	// add the donor/acceptor site to the vertices map
	vertexMap[cii->getStart()] = acceptor;
	vertexMap[cii->getEnd()] = donor;

	// link to the previous donor site
	// we create an edge between the previous exon
	// and the current exon

	if (lastDonorSite > 0) {

	  Edge e(lastDonorSite, cii->getStart());
	  edges.push_back(e);
	}

	// link acceptor site and donor site

	Edge e(cii->getStart(), cii->getEnd());
	edges.push_back(e);

	// finally the last donor site is the end of the current exon

	lastDonorSite = cii->getEnd();

      }

    //
    // we should have twice the number of exons in the vertexMap
    //

    assert(vertexMap.size() == 2*exonsT1.size());

    //
    // we should always have the number of vertices minus 1
    // edges in the graph
    //

    assert(vertexMap.size() == edges.size()+1);

    //
    // list all exons from transcript t2
    //

    std::vector<Exon> exonsT2 = t2.getExons();


    //
    // add all the exons of t2 to the vertexMap taking care
    // of not overwriting the exons of t1
    //


    // precondition
    // reset lastDonorSite

    lastDonorSite = 0;
    bool bLastDonorIsNovel = false;

    for(cii=exonsT2.begin(); cii!=exonsT2.end(); cii++)
      {

	bool bNewAcceptor = vertexMap.find(cii->getStart()) == vertexMap.end();
	bool bNewDonor = vertexMap.find(cii->getEnd()) == vertexMap.end();


	// look if the acceptor site exists:

	if(bNewAcceptor)
	  {

	    // new acceptor site
	    SplicingVertex acceptor(cii->getStart(), ACCEPTOR_SITE, vertexCount++);
	    genomicLocations.push_back(cii->getStart());
	    acceptor.addExon(*cii);
	    vertexMap[cii->getStart()] = acceptor;
	    //std::cout<< cii->getStart() << std::endl;

	  } else {

	    vertexMap[cii->getStart()].addExon(*cii);
	    assert(vertexMap[cii->getStart()].getExons().size() > 1);

	  }

	// look if the donor site exists:

	if(bNewDonor)
	  {

	    // new donor site
	    SplicingVertex donor(cii->getEnd(), DONOR_SITE, vertexCount++);
	    genomicLocations.push_back(cii->getEnd());
	    donor.addExon(*cii);
	    vertexMap[cii->getEnd()] = donor;
	    //std::cout << cii->getEnd() << std::endl;

	  } else {

	    vertexMap[cii->getEnd()].addExon(*cii);
	    assert(vertexMap[cii->getEnd()].getExons().size() > 1);

	  }

	//
	// if the lastDonor site was a new donor
	// or if the current acceptor site is new
	// link the last donor to this acceptor site
	//

	if (bLastDonorIsNovel || (bNewAcceptor && lastDonorSite > 0)) {

	  Edge e(lastDonorSite, cii->getStart());
	  edges.push_back(e);

	}

	if (bNewAcceptor || bNewDonor) {

	  // link the 2 nodes because one of them is new
	  Edge e(cii->getStart(), cii->getEnd());
	  edges.push_back(e);

	}

	// keep the last donor position
	lastDonorSite = cii->getEnd();
	bLastDonorIsNovel = bNewDonor;

      }

    //
    // now we have all the vertices and all the edges at disposal
    // we create a splicing graph
    // based on the number of vertex needed.

    std::cout << "Create graph of size " << (vertexMap.size()) << std::endl;

    return SplicingGraph(genomicLocations,vertexMap,edges);

  }

  SplicingEventVector SplicingGraph::computeSplicingEvents()
  {

    SplicingEventVector v;

    //
    // compute alternative initiation
    //

    if (inputs.size() == 2)
      {
	// create an alternative initiation event
	// if and only if the next donor site is
	// the same
	// if they are different, this is an alternative
	// first exon.

	//
	// check whether the donor site is the same
	//

	SplicingVertex v1 = inputs[0];
	SplicingVertexVector v1Out = getOutVertices(v1);

	SplicingVertex v2 = inputs[1];
	SplicingVertexVector v2Out = getOutVertices(v2);


	AlternativeInitiation ai();
	//ai.setStart
	//v.push_back(ai);

      }


    return v;
  }


}
