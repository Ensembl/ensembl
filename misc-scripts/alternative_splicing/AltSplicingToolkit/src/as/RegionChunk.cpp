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

#include <sstream>

#include <log4cpp/Category.hh>

#include "RegionChunk.h"

namespace as
{

  RegionChunk::RegionChunk() : as::Coordinates(0,0)
  {
  }

  RegionChunk::~RegionChunk()
  {
    // TODO Auto-generated destructor stub
  }

  void RegionChunk::setGene(const shared_ptr<Gene> &gene)
  {
    this->gene = gene;
  }

  const shared_ptr<Gene> &RegionChunk::getGene() const
  {
    return gene;
  }

  int RegionChunk::computeFeatureStart(int start, int end) const
  {
    return (strand == 1) ? start - referentialPosition + 1 : referentialPosition - end + 1;
  }

  int RegionChunk::computeFeatureEnd(int start, int end) const
  {
    return (strand == 1) ? end - referentialPosition + 1 : referentialPosition - start + 1;
  }

  void RegionChunk::mergeTranscript(const shared_ptr<Transcript> &transcript)
  {

    /*
     * Of course, strand and gene information ought to be the same each time
     * we merge a transcript.
     */

    this->strand = transcript->getStrand();
    this->gene = transcript->getGene();

    /*
     * Computes the reference position if this is the first time we insert a transcript
     * or the new reference position if we insert a latter transcript.
     */

    if (this->start == 0) {
      referentialPosition = (strand == 1) ? transcript->getStart() : transcript->getEnd();
    } else {
      referentialPosition = (strand == 1) ? min(transcript->getStart(), this->start) : max(transcript->getEnd(), this->end);
    }

    log4cpp::Category& root = log4cpp::Category::getRoot();
    root.infoStream() << "Merging transcript '" << transcript->getIdentifier() << "' to the region chunk." << log4cpp::eol;

    //
    // trivial case:
    // RegionChunk is empty, so we create an ExonChunk instance for each Exon
    //

    if (exonChunks.size() == 0) {

      /**
       * just insert each exon in turn without looking at the coordinates.
       * This is ordered by exon start positions.
       */

      vector< shared_ptr< TranscriptFeature> >&tFeatures = (*transcript).getTranscriptFeatures();
      vector< shared_ptr<TranscriptFeature> >::const_iterator tFeaturesIt;

      for( tFeaturesIt=tFeatures.begin(); tFeaturesIt!=tFeatures.end(); tFeaturesIt++ ) {

        if ((**tFeaturesIt).getType() == EXON_TYPE) {

          int eChunkStart = computeFeatureStart((**tFeaturesIt).getStart(), (**tFeaturesIt).getEnd()) ;
          int eChunkEnd = computeFeatureEnd((**tFeaturesIt).getStart(), (**tFeaturesIt).getEnd()) ;

          ExonChunk *eChunk = new ExonChunk( (**tFeaturesIt).getStart(), (**tFeaturesIt).getEnd() );
          eChunk->addExon(*tFeaturesIt);
          shared_ptr< ExonChunk > peChunk(eChunk);
          exonChunks.push_back(peChunk);

        }

      }

    } else {

      /**
       * We have to merge this new transcript with the existing region.
       * Conditions:
       * - no overlap => insert a new Exon Chunk at the right position (ordering by exon chunk start)
       * - overlap => merge the current exon with one or several Exon chunks and update the exon chunk coordinates
       */

      vector< shared_ptr< TranscriptFeature> >&tFeatures = (*transcript).getTranscriptFeatures();
      vector< shared_ptr<TranscriptFeature> >::const_iterator tFeaturesIt;

      for( tFeaturesIt=tFeatures.begin(); tFeaturesIt!=tFeatures.end(); tFeaturesIt++ ) {



        // look at the exon

        if ((**tFeaturesIt).getType() == EXON_TYPE) {

          // creates and prepares the next chunk

          unsigned int chunkStart =(**tFeaturesIt).getStart();
          unsigned int chunkEnd =(**tFeaturesIt).getEnd();

          ExonChunk *eChunk = new ExonChunk( chunkStart, chunkEnd ) ;
          eChunk->addExon(*tFeaturesIt);
          shared_ptr< ExonChunk > peChunk(eChunk);

          // transform the coordinates
          int featureStart = computeFeatureStart((**tFeaturesIt).getStart(), (**tFeaturesIt).getEnd()) ;
          int featureEnd = computeFeatureEnd((**tFeaturesIt).getStart(), (**tFeaturesIt).getEnd()) ;


          // initialize interval positions

          int startPos = -1;
          int endPos = -1;
          int index = 0;

          //cout << "Will merge the current exon" << endl;

	  /*
	   * Now iterate on the existing exonic chunks
	   */

          vector< shared_ptr<ExonChunk> >::const_iterator exonChunkIt;

          for( exonChunkIt=exonChunks.begin(); exonChunkIt!=exonChunks.end(); exonChunkIt++ ) {

            // transform the coordinates:
            int exonChunkStart = computeFeatureStart((**exonChunkIt).getStart(), (**exonChunkIt).getEnd()) ;
            int exonChunkEnd = computeFeatureEnd((**exonChunkIt).getStart(), (**exonChunkIt).getEnd()) ;

            // check whether we went upstream of the current featureEnd

            if (exonChunkStart > featureEnd) {

              // if we had already a startPos, we don't need to modify it
              // the condition is enough to break the iteration.
              if (startPos < 0) {
                startPos = index;
              }

              break;
            }

            // if not, look for overlaps

            if (!(exonChunkStart > featureEnd || exonChunkEnd < featureStart) ) {

	      root.infoStream() << "overlap between chunk " << exonChunkStart << "-" <<  exonChunkEnd << " and current feature " << featureStart << "-" << featureEnd ;

              if (startPos < 0) {
                startPos = index;
              }

              endPos = index;

              // because there is an overlap we can precompute the chunk extension in advance
              // but of course it depends on the strand

              if (this->strand == 1) {

                chunkStart = min(chunkStart, (**exonChunkIt).getStart());
                chunkEnd = max(chunkEnd, (**exonChunkIt).getEnd());

              } else {

                chunkStart = min(chunkStart, (**exonChunkIt).getStart());
                chunkEnd = max(chunkEnd, (**exonChunkIt).getEnd());

              }

              // Merge the existing chunk with the new one
	      // We will reinitialize the start and end later
              eChunk->mergeExons((**exonChunkIt).getExons());

            } // overlap

            index++;

           } // iterate on the exon chunks

          // ok erase first

          if (endPos >= 0) {

	    root.infoStream() << "erase from " << startPos << "to " << endPos << "(" << exonChunks.size() << ")";
            // erase and replace
	    if (startPos == endPos) {
	      exonChunks.erase (exonChunks.begin()+startPos);
	    } else {
	      exonChunks.erase (exonChunks.begin()+startPos,exonChunks.begin()+endPos+1);
	    }
	    root.infoStream() << "new vector size " << exonChunks.size();
            // reinit start end
            eChunk->setStart(chunkStart);
            eChunk->setEnd(chunkEnd);

          }

          // then insert a new exon chunk

          if (startPos >= 0) {

            // insert the exonChunk at the position startPos

            exonChunks.insert ( exonChunks.begin() + startPos , peChunk );
	    root.infoStream() << "insert at " << startPos << "(" << exonChunks.size() << ")";

          } else {

            // add a new exonChunk at the end of the structure because we couldn't find a position
	    root.infoStream() << "push back " << startPos << "(" << exonChunks.size() << ")";
            exonChunks.push_back ( peChunk );

          }

        } // EXON_TYPE

    } // for

  } // merge

    // now display all

    vector< shared_ptr<ExonChunk> >::const_iterator exonChunkIt;

    for( exonChunkIt=exonChunks.begin(); exonChunkIt!=exonChunks.end(); exonChunkIt++ ) {

      // transform the coordinates:
      int exonChunkStart = computeFeatureStart((**exonChunkIt).getStart(), (**exonChunkIt).getEnd()) ;
      int exonChunkEnd = computeFeatureEnd((**exonChunkIt).getStart(), (**exonChunkIt).getEnd()) ;

      root.infoStream() << "[" << exonChunkStart << "-" <<  exonChunkEnd << "]";

    }
    root.infoStream() << log4cpp::eol;

  }

  /**
   * Returns the computed alternative initiation events.
   */
  vector< shared_ptr< SplicingEvent > >& RegionChunk::getConstitutiveExonEvents()
  {
    return constitutiveExonEvents;
  }

  void  RegionChunk::checkConstitutiveExon(unsigned int upperBound) {

    log4cpp::Category& root = log4cpp::Category::getRoot();

    if (!constitutiveExonEvents.empty()) {
      constitutiveExonEvents.clear();
        }

      if (exonChunks.empty()) {
        root.error("Constitutive exon check: Empty data structure.");
      } else if (this->gene == 0) {
        root.error("Constitutive exon check: No gene information available.");
        return;
      }

      root.infoStream() << "Looking for constitutive exons on " << upperBound << " transcripts of gene '" <<  this->gene->getIdentifier() << "'." << log4cpp::eol;

      vector< shared_ptr<ExonChunk> >::const_iterator exonChunkIt;

      /**
       * Loop on all exon chunks from the region
       */

      for( exonChunkIt=exonChunks.begin(); exonChunkIt!=exonChunks.end(); exonChunkIt++ ) {

        root.infoStream() << "Checking exon chunk\t" <<  (**exonChunkIt).getStart() << "-" <<  (**exonChunkIt).getEnd() << log4cpp::eol;

        // rank each pair of exon coordinates from the chunk in a map
        map<string,int> rankMap;
        // push all exons in the indexMap
        map<string, list<int> > indexMap;

        vector< shared_ptr<TranscriptFeature> >::const_iterator exonIt;
        int index = 0;

        /**
         * Loop on all exons
         */
        for( exonIt=(**exonChunkIt).getExons().begin(); exonIt!=(**exonChunkIt).getExons().end(); exonIt++ ) {

          // add the current exon coordinates to the map if not existing and set the rank to 0

          std::ostringstream oss;
          oss << (**exonIt).getStart() << "-" << (**exonIt).getEnd();
          string mapKey = oss.str();

          if (rankMap.find(mapKey) == rankMap.end()) {
            rankMap[mapKey]=0;
          }
          rankMap[mapKey]++;
          indexMap[mapKey].push_back(index);

          /**
           * This is the condition for a constitutive event to occur:
           * There must be as many exons covering a position from the exon chunk 
	   * as many known transcripts in the genomic region.
           */
          if (rankMap[mapKey] == upperBound) {

            root.infoStream() << "Found " << upperBound << " times a constitutive exon location " << mapKey << log4cpp::eol;

            list<shared_ptr < TranscriptFeature > > constitutives;

            list< int >::const_iterator exonIndexIt = indexMap[mapKey].begin();
            int pos = (*exonIndexIt);
            shared_ptr< TranscriptFeature > exon = ((**exonChunkIt).getExons().at(pos));

            SplicingEvent *event =
              new SplicingEvent(
                  CNE_EVENT,
                  exon->getStart(),
                  exon->getEnd(),
                  exon->getChromosome(),
                  exon->getStrand()
              );

            event->setGene(this->gene);

            /* set constitutive exon sites */
            list<int> sites;
            sites.push_back(1);
            sites.push_back(exon->getStart());
            sites.push_back(exon->getEnd());
            event->setConstitutiveSites(sites);

            /**
             * Loop on the index of the exon list
             */


            for ( exonIndexIt=indexMap[mapKey].begin(); exonIndexIt!=indexMap[mapKey].end(); exonIndexIt++ ) {

              pos = (*exonIndexIt);
              exon= ((**exonChunkIt).getExons().at(pos));

              // look if it is not already in the constitutive list.
              list< shared_ptr< TranscriptFeature > >::const_iterator constitutiveIndexIt;
              bool alreadyPushed = false;
              for (constitutiveIndexIt=constitutives.begin(); constitutiveIndexIt!=constitutives.end(); constitutiveIndexIt++ ) {
                if ((**constitutiveIndexIt).getIdentifier().compare(exon->getIdentifier()) == 0) {
                  alreadyPushed = true;
                  break;
                }
              }

              if (!alreadyPushed) {
                // add the constitutive exons (check for identifier in Ensembl)
                constitutives.push_back(exon);
              }

            }

            event->setConstitutiveExons(constitutives);

            constitutiveExonEvents.push_back(shared_ptr<SplicingEvent>(event));

          } // condition to create a constitutive event

          index++;
        } // for each exon

      } // for each exon chunk

    }

  void RegionChunk::getEventsGffOutput(ostream &oStream, string datasource) const
    {
      int number = 1;
      vector< shared_ptr< SplicingEvent > >::const_iterator mii;
      for(mii=constitutiveExonEvents.begin(); mii!=constitutiveExonEvents.end(); mii++) {
          (*mii)->outputAsGFF(oStream, datasource, number);
          number++;
      }

    }

  void RegionChunk::getSummaryOutput(ostream &oStream) const
    {
      oStream << "constitutiveExonEvents:\t\t" << constitutiveExonEvents.size() << endl;
    }
}
