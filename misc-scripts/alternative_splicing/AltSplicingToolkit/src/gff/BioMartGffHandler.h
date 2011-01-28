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

#ifndef BIOMARTGFFHANDLER_H_
#define BIOMARTGFFHANDLER_H_

#include <string>
#include <fstream>
#include <iostream>
#include <map>

#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/shared_ptr.hpp>

#include <log4cpp/Category.hh>

#include "GffParser.h"
#include "GffEventModel.h"
#include "as/Feature.h"
#include "as/TranscriptFeature.h"
#include "as/Transcript.h"
#include "as/Gene.h"
#include "util/StringUtil.h"

using namespace std;
using namespace boost;
using namespace as;

namespace gff {

  const boost::regex eTab("\\t");
  const boost::regex eComma(";\\s*");
  const boost::regex eFeatureId("[a-z_]+\\s+\"(.+)\"");

  class BioMartGffHandler: public GffHandler, public GffNewGeneEventHandler {

  public:
    BioMartGffHandler();
    BioMartGffHandler(int limit);
    virtual ~BioMartGffHandler();

  public:
    void start();
    bool newline(string & str);
    void end();

    const map<string, shared_ptr<Exon> >&getExons();
    const map<string, shared_ptr<Transcript> >&getTranscripts();

  private:
    void fireNewgeneEvent();

  protected:
    int countExons;
    int countTranscripts;
    int column;

    int limit;
    int countGenes;

    string chr;
    unsigned int exonStart;
    unsigned int exonEnd;
    short int strand;

    string geneIdentifier;
    Gene *gene;
    bool newGene;

    string transcriptIdentifier;
    string exonIdentifier;

    string currentTranscript;
    string previousExon;

    map<string, shared_ptr<Exon> > exons;
    map<string, shared_ptr<Transcript> > transcripts;

    boost::sregex_token_iterator noMoreTokens;


  };

}

#endif /* BIOMARTGFFHANDLER_H_ */
