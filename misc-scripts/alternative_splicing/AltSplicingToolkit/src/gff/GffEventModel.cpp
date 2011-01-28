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

#include "GffEventModel.h"

namespace gff
{

  /**
   * GffEventListener methods
   */

  GffEventListener::~GffEventListener() {}


  /**
   * GffEventHandler methods
   */

  GffEventHandler::~GffEventHandler() {}


  GffNewGeneEvent::GffNewGeneEvent(map<string, shared_ptr<Transcript> > &transcripts, map<string, shared_ptr<Exon> > &exons) : transcripts(transcripts), exons(exons)
  {

  }

  GffNewGeneEvent::~GffNewGeneEvent()
  {
    //cout << "destroy GffNewGeneEvent" << endl;
  }

  map<string, shared_ptr<Transcript> > &GffNewGeneEvent::getTranscripts()
  {
    return transcripts;
  }

  map<string, shared_ptr<Exon> > &GffNewGeneEvent::getExons()
  {
    return exons;
  }

  /**
   * GffNewGeneEventListener methods
   */

  GffNewGeneEventListener::GffNewGeneEventListener() {}
  GffNewGeneEventListener::~GffNewGeneEventListener()
  {
    //cout << "destroy GffNewGeneEventListener" << endl;
 }

  /**
   * GffNewGeneEventHandler methods
   */
  GffNewGeneEventHandler::GffNewGeneEventHandler()
  {
  }

  GffNewGeneEventHandler::~GffNewGeneEventHandler()
  {
    // clear the listener list
    //cout << "destroy GffNewGeneEventHandler start: " << listeners.size() << endl;
    listeners.clear();
    //cout << "destroy GffNewGeneEventHandler end: " << listeners.size() << endl;
  }

  //void GffNewGeneEventHandler::registerNewGeneEventListener(const shared_ptr<GffNewGeneEventListener> & listener)
  void GffNewGeneEventHandler::registerNewGeneEventListener(GffNewGeneEventListener* listener)
  {
    listeners.push_back(listener);
  }

  void GffNewGeneEventHandler::triggerEvent(const shared_ptr<GffNewGeneEvent> &event)
  {
    //list< shared_ptr < GffNewGeneEventListener> >::iterator ii;
    list<  GffNewGeneEventListener* >::iterator ii;

    for(ii=listeners.begin(); ii!=listeners.end(); ii++) {

      (*ii)->command(event);
    }

  }
}
