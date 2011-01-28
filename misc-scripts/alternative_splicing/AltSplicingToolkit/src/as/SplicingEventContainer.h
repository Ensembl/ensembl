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

#ifndef SPLICINGEVENTCONTAINER_H_
#define SPLICINGEVENTCONTAINER_H_

#include <iostream>
#include <string>
#include <algorithm>

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "SplicingEvent.h"

using namespace boost;

namespace as
{

  class SplicingEventContainer
  {
  public:
    SplicingEventContainer();
    virtual
    ~SplicingEventContainer();

    vector< shared_ptr< SplicingEvent > >&  getIntronRetentionEvents();
    vector< shared_ptr< SplicingEvent > >&  getIntronIsoformEvents();
    vector< shared_ptr< SplicingEvent > >&  getExonIsoformEvents();
    vector< shared_ptr< SplicingEvent > >&  getMutuallyExclusiveEvents();
    vector< shared_ptr< SplicingEvent > >&  getCassetteExonEvents();

    vector< shared_ptr< SplicingEvent > >&  getAlternativeInitiationEvents();
    vector< shared_ptr< SplicingEvent > >&  getAlternativeTerminationEvents();

    vector< shared_ptr< SplicingEvent > >&  getAlternativeFirstExonEvents();
    vector< shared_ptr< SplicingEvent > >&  getAlternativeLastExonEvents();

    void addNewEvent(const shared_ptr<SplicingEvent> &event);
    void mergeSplicingEvents(SplicingEventContainer &eventContainer);

    void getSummaryOutput(ostream &oStream) const;
    void getGffOutput(ostream &oStream, string datasource) const;

    int getEventCount() const;

  protected:

    inline void getEventsGffOutput(ostream &oStream, string datasource, const vector< shared_ptr< SplicingEvent > >& events) const;
    void mergeSplicingEventVectors(vector< shared_ptr<SplicingEvent> > &eventSetA, vector< shared_ptr<SplicingEvent> > &eventSetB);

    vector< shared_ptr< SplicingEvent > > intronRetentionEvents;
    vector< shared_ptr< SplicingEvent > > intronIsoformEvents;
    vector< shared_ptr< SplicingEvent > > exonIsoformEvents;
    vector< shared_ptr< SplicingEvent > > mutuallyExclusiveEvents;
    vector< shared_ptr< SplicingEvent > > cassetteExonEvents;

    vector< shared_ptr< SplicingEvent > > alternativeInitiationEvents;
    vector< shared_ptr< SplicingEvent > > alternativeTerminationEvents;

    vector< shared_ptr< SplicingEvent > > alternativeFirstExonEvents;
    vector< shared_ptr< SplicingEvent > > alternativeLastExonEvents;

  };

}

#endif /* SPLICINGEVENTCONTAINER_H_ */
