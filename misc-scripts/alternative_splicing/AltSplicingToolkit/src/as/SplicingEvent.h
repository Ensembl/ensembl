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

#ifndef _SPLICING_EVENT_H
#define _SPLICING_EVENT_H

#include <utility> // for pair
#include <algorithm> // to merge 2 vectors
#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/utility.hpp>

#include "GeneFeature.h"
#include "Transcript.h"
#include "TranscriptFeature.h"

using namespace std;
using namespace boost;

namespace as {

  const int NO_EVENT = 2;
  const int MXE_EVENT = 4;   // mutually exclusive
  const int CE_EVENT = 5;    // cassette exon / skipped exon
  const int IR_EVENT = 6;    // intron retention
  const int II_EVENT = 7;    // intron isoform
  const int EI_EVENT = 8;    // exon isoform
  const int A3SS_EVENT = 9;  // Alternative 3' splice site
  const int A5SS_EVENT = 10;  // Alternative 5' splice site

  const int ALT_FIRST_LAST_EXON_EVENT = 16;  // first last exon events
  const int AI_EVENT = 17;  // alternative initiation
  const int AT_EVENT = 18;  // alternative termination
  const int AFE_EVENT = 19; // alternative first exon
  const int ALE_EVENT = 20; // alternative last exon

  const int CONSTITUTIVE_EVENT = 32; // constitutive region/exon event
  const int CNE_EVENT = 33; // constitutive exon
  const int CNR_EVENT = 34; // constitutive region

  class SplicingEvent : public GeneFeature
 {

  public:

    SplicingEvent();
    SplicingEvent(int type, unsigned int start, unsigned int end, string chr, short int strand);
    virtual ~SplicingEvent();

    void setSetA(const list< shared_ptr< TranscriptFeature > > &a);
    void setSetB(const list< shared_ptr< TranscriptFeature > > &b);
    void setConstitutiveExons(const list< shared_ptr< TranscriptFeature > > &constitutives);

    void setSitesA(const list< int > &a);
    void setSitesB(const list< int > &b);
    void setConstitutiveSites(const list< int > &c);

    bool equals(const SplicingEvent &event) const;
    bool compareFeatures(const list< shared_ptr< TranscriptFeature > > &set1, const list< shared_ptr< TranscriptFeature > > &set2 ) const;
    bool compareFeatureSites(const list< int > &sites1, const list< int > &sites2 ) const;

    const list< shared_ptr< TranscriptFeature > >& getSetA() const;
    const list< shared_ptr< TranscriptFeature > >& getSetB() const;
    const list< shared_ptr< TranscriptFeature > >& getConstitutiveExons() const;

    const list< int >& getSitesA() const;
    const list< int >& getSitesB() const;
    const list< int >& getConstitutiveSites() const;

    void outputAsGFF(ostream &oStream, string datasource, int number) const;
    string getTypeAsString() const;

    bool contains(vector< shared_ptr< SplicingEvent > > &events) const;
    vector< shared_ptr< SplicingEvent > >::iterator find(vector< shared_ptr< SplicingEvent > > &events) const;
    void addTranscriptPair(const shared_ptr< Transcript > & t1, const shared_ptr< Transcript > & t2);
    void mergeTranscriptPairs(const vector< pair< shared_ptr< Transcript >, shared_ptr< Transcript > > > &pairs);
    const vector< pair< shared_ptr< Transcript >, shared_ptr< Transcript > > > &getTranscriptPairs() const;
    void mergeCoordinates(unsigned int start, unsigned int end);
    void mergeSetA(const list< shared_ptr< TranscriptFeature > > &a);
    void mergeSetB(const list< shared_ptr< TranscriptFeature > > &a);
    void mergeConstitutiveExons(const list< shared_ptr< TranscriptFeature > > &constitutives);

  private:

    void mergeSets(list< shared_ptr< TranscriptFeature > > &x, const list< shared_ptr< TranscriptFeature > > &y);

    string findInvolvedTranscriptIdentifier(ostream &oStream, list<string> &transcriptList, TranscriptFeature &feature) const;

  protected:

    //TranscriptSets sets;
    list< shared_ptr< TranscriptFeature > > setA;
    list< shared_ptr< TranscriptFeature > > setB;
    list< shared_ptr< TranscriptFeature > > constitutiveExons;

    list< int > sitesA;
    list< int > sitesB;
    list< int > constitutiveSites;

    vector< pair< shared_ptr< Transcript >, shared_ptr< Transcript > > > transcriptPairs;

  };


  class SplicingEventNotFoundException : public std::exception
  {
    public:
      SplicingEventNotFoundException(string m="Splicing event not found") : msg(m) {}
      ~SplicingEventNotFoundException() throw() {}
      const char* what() const throw() { return msg.c_str(); }

    private:
      string msg;

  };

  class AlternativeInitiation : public SplicingEvent
    {

    public:
      AlternativeInitiation();
      ~AlternativeInitiation();

    };

  class AlternativeTermination : public SplicingEvent
    {

    public:
      AlternativeTermination();
      ~AlternativeTermination();

    };

  typedef vector<SplicingEvent> SplicingEventVector;

}



#endif /* not defined _SPLICING_EVENT_H */
