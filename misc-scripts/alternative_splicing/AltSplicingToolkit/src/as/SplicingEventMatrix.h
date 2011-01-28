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

#ifndef SPLICINGEVENTMATRIX_H_
#define SPLICINGEVENTMATRIX_H_

#include <iostream>
#include <algorithm>

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_list.hpp>

#include "TranscriptFeature.h"
#include "Transcript.h"
#include "SplicingEvent.h"
#include "SplicingEventContainer.h"

using namespace std;
using namespace boost;

namespace as {

  // code for overlaps type
  const short int NO_OVERLAP = 0x0; // 0        0000
  const short int ID5P_ID3P = 0xf; //15         1111
  const short int DIFF5P_ID3P = 0xe; //14       1110
  const short int ID5P_DIFF3P = 0xd; //13       1101
  const short int PART_OF = 0xc; // 12          1100
  const short int OVERLAP = 0x4; // 4           0100

  // bitwise constants
  const short int ID5P = 0x1; //1               0001
  const short int ID3P = 0x2; //2               0010

  class SplicingEventMatrix: public SplicingEventContainer
  {
  public:
    SplicingEventMatrix(const shared_ptr<Transcript> & t1, const shared_ptr<Transcript> & t2);
    virtual
    ~SplicingEventMatrix();

    static short int overlapCode(int start1, int end1, int start2, int end2);

    void computeSplicingEvents(bool bRELAX);

  protected:
    void buildOverlapMatrix();


  private:
    static void create_matrix(short int ***m, int xSize, int ySize);
    static void delete_matrix(short int ***m, int xSize, int ySize);
    static void display_matrix(short int ***m, int xSize, int ySize);

    void checkAlternativeFirstLastExon(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index);
    void checkIntronRetention(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index);
    void checkIntronIsoform(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index);
    void checkExonIsoform(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index, bool bRELAX);
    void checkCassetteExon(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index, bool bRELAX);
    void checkMutualExclusion(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index, bool bRELAX);

    inline int computeFeatureStart(int start, int end) const;
    inline int computeFeatureEnd(int start, int end) const;

  protected:
    short int **overlapMatrix; // is allocated here and released here
    shared_ptr<Transcript> t1; // pointer to a transcript defined elsewhere, no worry about the allocation
    shared_ptr<Transcript> t2; // pointer to a transcript defined elsewhere, no worry about the allocation
    vector< shared_ptr<TranscriptFeature> > &t1Features;
    vector< shared_ptr<TranscriptFeature> > &t2Features;
    int xSize;
    int ySize;
    int strand;
    int referentialPosition;

  };

}

#endif /* SPLICINGEVENTMATRIX_H_ */
