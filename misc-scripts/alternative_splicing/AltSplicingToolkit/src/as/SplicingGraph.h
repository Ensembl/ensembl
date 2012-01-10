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

#ifndef _SPLICING_GRAPH_H
#define _SPLICING_GRAPH_H

#include <utility>                   // for std::pair
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <map>
#include <utility>                   // for std::pair
#include <algorithm>

#include "Exon.h"
#include "Transcript.h"
#include "SplicingVertex.h"
#include "SplicingEvent.h"

#include <boost/multi_array.hpp>

using namespace boost;
using namespace as;

namespace as {

  //class SplicingVertex;

  // graph definition
  //typedef adjacency_matrix<directedS> Graph;
	/*
	 * An edge is simple a pair of integer position in the
	 * adjacency matrix.
	 */
  typedef std::pair<int,int> Edge;

  /*
   * A set of edges to return to the caller
   * is a vector of Edge instances.
   */
  typedef std::vector<Edge> Edges;
  typedef std::map<int, SplicingVertex> SplicingVertexMap;
  typedef std::vector<SplicingVertex> SplicingVertexVector;
  typedef boost::multi_array<bool, 2> BooleanMatrix;
  typedef BooleanMatrix::index BooleanIndex;
  //BooleanMatrix::extent_gen booleanExtents;

  class SplicingGraph {

    static int objectCount;

  public:
    SplicingGraph();
    ~SplicingGraph();

  protected:

    SplicingGraph(std::vector<int> gLocations, SplicingVertexMap vMap, Edges edges);
    void computeInputsOutputs();
    const SplicingVertexVector intersect(SplicingVertexVector setA, SplicingVertexVector setB);

  public:

    static SplicingGraph buildGraph(const Transcript& t1, const Transcript& t2);
    static void printMsg(const std::string& msg = "");

    SplicingVertexVector getInputs();
    SplicingVertexVector getOutputs();

    int getInDegree(unsigned int vertex);
    int getOutDegree(unsigned int vertex);

    SplicingVertexVector getInVertices(const SplicingVertex& vertex);
    SplicingVertexVector getOutVertices(const SplicingVertex& vertex);


    static void create_matrix(bool ***m, int size);
    static void delete_matrix(bool ***m, int size);
    static void display_matrix(bool ***m, int size);

    SplicingEventVector computeSplicingEvents();

  protected:

    std::vector<int> genomicLocations; // genomic locations of the vertices
    SplicingVertexMap vertexMap;       // map of all the vertices
    Edges edges;

    bool **adjacencyMatrix;            // 2-dimension array (matrix)
    unsigned int size;                 // number of vertices

    SplicingVertexVector inputs;       // inputs of the graph
    SplicingVertexVector outputs;      // outputs of the graph

  };

}

#endif /* not defined _SPLICING_GRAPH_H */
