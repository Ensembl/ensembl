/*
 *    AltSplicingToolkit
 *    Author: Gautier Koscielny <koscieln@ebi.ac.uk>
 *
 *    Copyright (c) 1999-2011 The European Bioinformatics Institute and 
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

#ifndef REGIONCHUNK_H_
#define REGIONCHUNK_H_

#include <vector>
#include <list>
#include <string>

#include <boost/shared_ptr.hpp>

#include "Coordinates.h"
#include "Transcript.h"
#include "ExonChunk.h"
#include "SplicingEvent.h"
#include "Gene.h"

namespace as
{

  class RegionChunk : public as::Coordinates
  {
  public:
    RegionChunk();
    virtual
    ~RegionChunk();

    //void addExonChunk(const shared_ptr<ExonChunk> &exonChunk);
    void mergeTranscript(const shared_ptr<Transcript> &transcript);
    vector< shared_ptr< SplicingEvent > >&  getConstitutiveExonEvents();
    void checkConstitutiveExon(unsigned int upperBound);
    void getEventsGffOutput(ostream &oStream, string datasource) const;
    void getSummaryOutput(ostream &oStream) const;

    int computeFeatureStart(int start, int end) const;
    int computeFeatureEnd(int start, int end) const;

    void setGene(const shared_ptr<Gene> &gene);
    const shared_ptr<Gene> &getGene() const;

  protected:
    shared_ptr<Gene> gene;
    vector< shared_ptr<ExonChunk> > exonChunks;
    int strand;
    int referentialPosition;

    vector< shared_ptr< SplicingEvent > > constitutiveExonEvents;

  };

}

#endif /* REGIONCHUNK_H_ */
