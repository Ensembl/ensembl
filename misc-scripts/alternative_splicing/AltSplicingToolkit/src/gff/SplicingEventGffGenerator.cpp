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

#include <log4cpp/Category.hh>

#include "SplicingEventGffGenerator.h"
#include "as/RegionChunk.h"


namespace gff
{

  SplicingEventGffGenerator::SplicingEventGffGenerator(ostream &oStream, string datasource, bool bCNE, bool bRELAX) : oStream(oStream), datasource(datasource)
  {
    countAI = 0;
    countAT = 0;
    countAFE = 0;
    countALE = 0;
    countIR = 0;
    countII = 0;
    countEI = 0;
    countCE = 0;
    countMXE = 0;

    genes = 0;
    genesWithEvents = 0;
    genesWithSeveralTranscripts = 0;

    this->bCNE = bCNE;
		this->bRELAX = bRELAX;

  }

  SplicingEventGffGenerator::~SplicingEventGffGenerator()
  {
    //cout << "destroy SplicingEventGffGenerator" << endl;
  }

  void SplicingEventGffGenerator::command(const shared_ptr<GffNewGeneEvent> &triggeredEvent)
  {

    log4cpp::Category& root = log4cpp::Category::getRoot();

    map<string, shared_ptr<Transcript> > transcripts = triggeredEvent->getTranscripts();
    map<string, shared_ptr<Exon> > exons = triggeredEvent->getExons();

    /***********************************************************************************/
    /* FIND CONSTITUTIVE EXONS                                                         */
    /***********************************************************************************/

    /**
     * The following data structure will contain a compact representation of
     * the region convenient to determine the list of constitutive exon
     * This could also work with a splicing graph structure
     */
    RegionChunk chunk;

    /**
     * Loop over all the transcripts and compare them to find the constitutive pieces of the gene
     * Find also the constitutive exons / region
     */
    for( map<string, shared_ptr<Transcript> >::const_iterator ii=transcripts.begin(); ii!=transcripts.end(); ++ii) {

      shared_ptr<Transcript> transcript = ii->second;
      chunk.mergeTranscript(transcript);

    }

    chunk.checkConstitutiveExon(transcripts.size());
    chunk.getEventsGffOutput(oStream, datasource);
    countCNE += chunk.getConstitutiveExonEvents().size();

    //return;



    /***********************************************************************************/
    /* FIND ALTERNATIVE SPLICING EVENTS                                                */
    /***********************************************************************************/

    /**
     * The following data structure will contain all the
     * splicing events classified by type (Exon Isoform, Cassette Exon, etc.)
     * allocated on the stack
     */
    SplicingEventContainer container;

    if (!bCNE) {

    /**
     * Loop over all the transcripts and compare them to find alternative events

     */
    int max2 = 0;

    for( map<string, shared_ptr<Transcript> >::const_iterator ii=transcripts.begin(); ii!=transcripts.end(); ++ii)
      {
        // reinit index2
        int index2 = 0;


        for( map<string, shared_ptr<Transcript> >::const_iterator ij=transcripts.begin(); ij!=transcripts.end(); ++ij)
          {

            // skip the comparisons already done (half the space)
            while (index2 <= max2) {
              index2++;
              ++ij;
            }

            if (ij == transcripts.end())
              break;

            shared_ptr<Transcript> t1 = ii->second;
            shared_ptr<Transcript> t2 = ij->second;

            if (t1->getIdentifier().compare(t2->getIdentifier()) != 0)
               // && t1->getIdentifier().compare("ENST00000353540") == 0
               // && t2->getIdentifier().compare("ENST00000357654") == 0)
              {
                root.infoStream() << t1->getIdentifier() + " <=> " + t2->getIdentifier() << log4cpp::eol;

                // build a splicing matrix
                //cout << "build the splicing matrix" << endl;
                SplicingEventMatrix matrix(t1, t2);

                // compute splicing events
                //cout << "compute splicing events" << endl;
                matrix.computeSplicingEvents(bRELAX);

                // store the new events in the container.
                //cout << "store the new events in the container" << endl;
                container.mergeSplicingEvents(matrix);

		//container.getSummaryOutput(cerr);

              }

          }

          max2++;

      }

    container.getGffOutput(oStream, datasource);

    }

    // cumulative stats

    genes++; // number of genes parsed from the input stream

    // number of genes with several transcripts

    if (transcripts.size() > 1) {

      genesWithSeveralTranscripts++;

    }

    // genes with events

    if (container.getEventCount() > 0) {

      genesWithEvents++;

      countAI += container.getAlternativeInitiationEvents().size();
      countAT += container.getAlternativeTerminationEvents().size();
      countAFE += container.getAlternativeFirstExonEvents().size();
      countALE += container.getAlternativeLastExonEvents().size();
      countIR += container.getIntronRetentionEvents().size();
      countII += container.getIntronIsoformEvents().size();
      countEI += container.getExonIsoformEvents().size();
      countCE += container.getCassetteExonEvents().size();
      countMXE += container.getMutuallyExclusiveEvents().size();

    }
  }

  int SplicingEventGffGenerator::getCountAI() const { return countAI; }
  int SplicingEventGffGenerator::getCountAT() const { return countAT; }
  int SplicingEventGffGenerator::getCountAFE() const { return countAFE; }
  int SplicingEventGffGenerator::getCountALE() const { return countALE; }
  int SplicingEventGffGenerator::getCountIR() const { return countIR; }
  int SplicingEventGffGenerator::getCountII() const { return countII; }
  int SplicingEventGffGenerator::getCountEI() const { return countEI; }
  int SplicingEventGffGenerator::getCountCE() const { return countCE; }
  int SplicingEventGffGenerator::getCountMXE() const { return countMXE; }
  int SplicingEventGffGenerator::getCountCNE() const { return countCNE; }
  int SplicingEventGffGenerator::getGeneCount() const { return genes; }
  int SplicingEventGffGenerator::getGenesWithEventsCount() const { return genesWithEvents; }
  int SplicingEventGffGenerator::getGenesWithSeveralTranscriptsCount() const { return genesWithSeveralTranscripts; }
  int SplicingEventGffGenerator::getEventCount() const
  {
    // don't sum the constitutive exon events
    return (countAI+countAT+countAFE+countALE+countIR+countII+countEI+countCE+countMXE);
  }

}
