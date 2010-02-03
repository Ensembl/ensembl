/*
 *    AltSplicingToolkit
 *    Author: Gautier Koscielny <koscieln@ebi.ac.uk>
 *
 *    Copyright (c) 1999-2010 The European Bioinformatics Institute and 
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

#include "Feature.h"
#include "Transcript.h"
//#include "Exon.h"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

namespace as {

  //: Feature::Feature()

  Transcript::Transcript()
  {
    type = TRANSCRIPT_TYPE;
  }

  Transcript::Transcript(std::string featureIdentifier) : GeneFeature(featureIdentifier)
  {
    type = TRANSCRIPT_TYPE;
  }

  void Transcript::addExon(const shared_ptr<Exon> &exon)
  {
    //if (strand == 1) {
      exons.push_back(exon);
    //} else {
   //   exons.push_front(exon);
   // }

    if (exon->getStart() < start) {
      start = exon->getStart();
    }
    if (exon->getEnd() > end) {
      end = exon->getEnd();
    }

  }

  vector< shared_ptr< Exon > >& Transcript::getExons()
  {
    vExons.clear();
    vector< shared_ptr< Exon > >::iterator it;
    vExons.insert(it, exons.begin(), exons.end());
    return vExons;
  }

  /**
   * Create a full structure including the introns.
   */
  vector< shared_ptr< TranscriptFeature> >& Transcript::getTranscriptFeatures() {

    int lastExonStart = 0;
    int lastExonEnd = 0;
    //int exonIndex = (strand == 1) ? 0 : exons.size()+1;
    int exonIndex = 0;

    // create a vector of exon - intron on the heap
    features.clear();

    list< shared_ptr< Exon > >::const_iterator exonIterator;

    for(exonIterator=exons.begin(); exonIterator!=exons.end(); exonIterator++)
      {

        // current exon


        int exonStart = (**exonIterator).getStart();
        int exonEnd = (**exonIterator).getEnd();

        // if there is an exon before the current exon
        if (lastExonStart > 0) {

          ostringstream osstream;
          osstream << "intron" << exonIndex << "-" << (exonIndex+1);
          std::string my_id = osstream.str();

          // create the intron on the heap
          Intron *intron =  new Intron(my_id);
          // get the donor site of the previous exon
          intron->setStart((strand == 1) ? lastExonEnd + 1 : exonEnd + 1 );
          // get the acceptor site of th current exon
          intron->setEnd((strand == 1) ? exonStart - 1 : lastExonStart - 1 );
          // get the strand and chromosome of the previous feature
          intron->setStrand(getStrand());
          intron->setChromosome(getChromosome());

          shared_ptr< TranscriptFeature > pIntron(intron);
          features.push_back(pIntron);
        }

        // push the current exon shared pointer
        features.push_back(*exonIterator);

        lastExonStart = exonStart;
        lastExonEnd = exonEnd;

        exonIndex++; // = exonIndex + ((strand == 1) ? 1 : -1);

      }
    return features;
  }

  Transcript::~Transcript()
  {
    //cout << "Calling transcript destructor: " << identifier << endl;
  }

}
