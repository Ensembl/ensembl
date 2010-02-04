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

#ifndef GFFLISTENER_H_
#define GFFLISTENER_H_

#include <list>
#include <map>
#include <iostream>

#include <boost/shared_ptr.hpp>

#include "as/TranscriptFeature.h"
#include "as/Transcript.h"


using namespace std;
using namespace boost;
using namespace as;

namespace gff
{

  class GffEventListener
  {
    public: virtual ~GffEventListener();
  };


  class GffEventHandler
  {
    public: virtual ~GffEventHandler();
  };

  class GffNewGeneEvent
  {
  public:
    GffNewGeneEvent(map<string, shared_ptr<Transcript> > &transcripts, map<string, shared_ptr<Exon> > &exons);
    virtual ~GffNewGeneEvent();

    map<string, shared_ptr<Transcript> > &getTranscripts();
    map<string, shared_ptr<Exon> > &getExons();

  protected:

    map<string, shared_ptr<Transcript> > &transcripts;
    map<string, shared_ptr<Exon> > &exons;

  };

  /**
   * ABC (Abstract Base Class)
   */
  class GffNewGeneEventListener
  {
  public:
    GffNewGeneEventListener();
    virtual ~GffNewGeneEventListener();
    virtual void command(const shared_ptr<GffNewGeneEvent> &event ) = 0; // {cout << "GffNewGeneEventListener::command" << endl; } // Pure virtual function, no body.
  };

  class GffNewGeneEventHandler
  {
  public:
    GffNewGeneEventHandler();
    virtual
    ~GffNewGeneEventHandler();

    //void registerNewGeneEventListener(const shared_ptr<GffNewGeneEventListener> & listener);
    void registerNewGeneEventListener(GffNewGeneEventListener* listener);
    void triggerEvent(const shared_ptr<GffNewGeneEvent> &event);

  private:
    //list< shared_ptr<GffNewGeneEventListener> > listeners;
    list< GffNewGeneEventListener* > listeners;

  };


}

#endif /* GFFLISTENER_H_ */
