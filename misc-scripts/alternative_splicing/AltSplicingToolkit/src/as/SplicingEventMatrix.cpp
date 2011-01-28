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

#include "SplicingEventMatrix.h"
#include <log4cpp/Category.hh>
#include <sstream>

namespace as {

  SplicingEventMatrix::SplicingEventMatrix(const shared_ptr<Transcript> & t1, const shared_ptr<Transcript> & t2) : t1(t1), t2(t2), t1Features((*t1).getTranscriptFeatures()), t2Features((*t2).getTranscriptFeatures())
  {
    xSize = 0;
    ySize = 0;
    strand = t1->getStrand();
    referentialPosition = (strand == 1) ? min(t1->getStart(), t2->getStart()) : max(t1->getEnd(), t2->getEnd());

    buildOverlapMatrix();

  }

  SplicingEventMatrix::~SplicingEventMatrix()
  {
    delete_matrix(&overlapMatrix, xSize, ySize);
  }

  int SplicingEventMatrix::computeFeatureStart(int start, int end) const
  {
    return (strand == 1) ? start - referentialPosition + 1 : referentialPosition - end + 1;
  }

  int SplicingEventMatrix::computeFeatureEnd(int start, int end) const
  {
    return (strand == 1) ? end - referentialPosition + 1 : referentialPosition - start + 1;
  }

  void SplicingEventMatrix::buildOverlapMatrix() {

    xSize = t1Features.size();
    ySize = t2Features.size();

    //cout << "xSize=" << xSize << ";ySize=" << ySize << endl;

    create_matrix(&overlapMatrix, xSize, ySize);

    //display_matrix(&overlapMatrix, xSize, ySize);

    // now compute overlaps O(n2), could be reduce
    // loop over the transcript features


    int t1Index = 0;
    vector< shared_ptr<TranscriptFeature> >::const_iterator t1FeaturesIt;

    for( t1FeaturesIt=t1Features.begin(); t1FeaturesIt!=t1Features.end(); t1FeaturesIt++ )
      {
        int t1FeatureStart = computeFeatureStart((**t1FeaturesIt).getStart(), (**t1FeaturesIt).getEnd()) ;
        int t1FeatureEnd = computeFeatureEnd((**t1FeaturesIt).getStart(), (**t1FeaturesIt).getEnd()) ;
        string t1Identifier = (**t1FeaturesIt).getIdentifier();

        //cout << t1Index << " " << t1->getIdentifier() << "::" << t1Identifier  << " " << t1FeatureStart << "-" << t1FeatureEnd << endl;

        int t2Index = 0;
        vector< shared_ptr<TranscriptFeature> >::const_iterator t2FeaturesIt;

        for(t2FeaturesIt=t2Features.begin(); t2FeaturesIt!=t2Features.end(); t2FeaturesIt++)
          {
            // compute the overlap code with the current exon
            int t2FeatureStart = computeFeatureStart((**t2FeaturesIt).getStart(), (**t2FeaturesIt).getEnd()) ;
            int t2FeatureEnd = computeFeatureEnd((**t2FeaturesIt).getStart(), (**t2FeaturesIt).getEnd()) ;
            string t2Identifier = (**t2FeaturesIt).getIdentifier();

            short int code = overlapCode(t1FeatureStart, t1FeatureEnd, t2FeatureStart, t2FeatureEnd);

            //cout << "\t\t" << t2Index << " " << t2->getIdentifier() << "::" <<  t2Identifier << " " << t2FeatureStart << "-" << t2FeatureEnd << " " << code << endl;
            if (code != NO_OVERLAP) {
              // update the matrix
              overlapMatrix[t1Index][t2Index] = code;
            }
            t2Index++;
          }

        t1Index++;
      }

    // display the matrix again
    display_matrix(&overlapMatrix, xSize, ySize);
  }

  void SplicingEventMatrix::computeSplicingEvents(bool bRELAX) {

    log4cpp::Category& root = log4cpp::Category::getRoot();

    // first get strand and chromosome
    //std::string chr = t1.getChromosome();
    //int strand = t1.getStrand();

    int t1Index = 0;


    vector< shared_ptr<TranscriptFeature> >::iterator t1FeaturesIt;

    for(t1FeaturesIt=t1Features.begin(); t1FeaturesIt!=t1Features.end(); t1FeaturesIt++) {

      int t2Index = 0;
      vector< shared_ptr<TranscriptFeature> >::iterator t2FeaturesIt;

      for(t2FeaturesIt=t2Features.begin(); t2FeaturesIt!=t2Features.end(); t2FeaturesIt++) {

        if (overlapMatrix[t1Index][t2Index] != NO_OVERLAP) {

          root.infoStream() << "* [" << (**t1FeaturesIt).getIdentifier() << "<=>" << (**t2FeaturesIt).getIdentifier() << "]";

          if (overlapMatrix[t1Index][t2Index] != ID5P_ID3P) {

            root.infoStream() << " OVERLAP" << log4cpp::eol;
            /*
             * look at the type of the first and the second feature compared
             */

            if ((**t1FeaturesIt).getType() == EXON_TYPE && (**t2FeaturesIt).getType() == EXON_TYPE) {

              /*
               * exon isoform/alternative initiation/alternative polyadenylation
               * exon isoform: two exons flanked by overlapping introns:
               *
               */
              checkExonIsoform(*t1FeaturesIt, *t2FeaturesIt, t1Index, t2Index, bRELAX);

            } else if ((**t1FeaturesIt).getType() == INTRON_TYPE && (**t2FeaturesIt).getType() == INTRON_TYPE) {

              /*
               * intron isoform
               * two introns overlapping with flanking exons overlapping
               */
              checkIntronIsoform(*t1FeaturesIt, *t2FeaturesIt, t1Index, t2Index);

            } else if (((**t1FeaturesIt).getType() == EXON_TYPE && (**t2FeaturesIt).getType() == INTRON_TYPE) || ((**t1FeaturesIt).getType() == INTRON_TYPE && (**t2FeaturesIt).getType() == EXON_TYPE)) {

              /*
               * one of the feature is part of the other
               * 3 possible outcomes: intron retention / cassette exon / mutually exclusive
               */
              if (overlapMatrix[t1Index][t2Index] == PART_OF) {

                int t1Start = computeFeatureStart((**t1FeaturesIt).getStart(), (**t1FeaturesIt).getEnd()) ;
                int t2Start = computeFeatureStart((**t2FeaturesIt).getStart(), (**t2FeaturesIt).getEnd()) ;

                if (t1Start < t2Start && (**t1FeaturesIt).getType() == EXON_TYPE ||
                    t2Start < t1Start && (**t2FeaturesIt).getType() == EXON_TYPE) {

                      checkIntronRetention(*t1FeaturesIt, *t2FeaturesIt, t1Index, t2Index);

                } else {

                  /*
                   * seek cassette exons or mutually exclusive exons
                   * condition for cassette exon:
                   * there exist a 5p/3p flanking exon overlapping with
                   */

                  checkCassetteExon(*t1FeaturesIt, *t2FeaturesIt, t1Index, t2Index, bRELAX);

                  if ((**t1FeaturesIt).getType() == EXON_TYPE) {

                    checkMutualExclusion(*t1FeaturesIt, *t2FeaturesIt, t1Index, t2Index, bRELAX);

                  }


                }

              } // part of

            } // exon / intron comparison
          }

          else  {

            root.infoStream() << " IDENTICAL" << log4cpp::eol;

            // track AFE/ALE events
            if ((**t1FeaturesIt).getType() == EXON_TYPE && (**t2FeaturesIt).getType() == EXON_TYPE) {

              /*
               * alternative first / alternative last exon
               */
              checkAlternativeFirstLastExon(*t1FeaturesIt, *t2FeaturesIt, t1Index, t2Index);
            }
          } // the features are identical

        } // there is an overlap between 2 features

        // move downstream
        t2Index++;
      }

      // move downstream
      t1Index++;
    }

    //getGffOutput(cout, "Ensembl");

  }

  short int SplicingEventMatrix::overlapCode(int start1, int end1, int start2, int end2) {
    // compute min start1 start2
    short int code = NO_OVERLAP;
    int maxStart = max(start1, start2);
    int minEnd = min(end1, end2);
    if (maxStart <= minEnd) {

      if (start1 == start2) {
        if (end1 == end2) {
          code = ID5P_ID3P;
        } else {
          code = ID5P_DIFF3P;
        }
      } else {
        if (end1 == end2) {
          code = DIFF5P_ID3P;
        } else {
          if ((start1 < start2 && end1 > end2) || (start2 < start1 && end2 > end1)) {
            code = PART_OF;
          } else {
            code = OVERLAP;
          }
        }
      }
    }
    return code;
  }

  void SplicingEventMatrix::create_matrix(short int ***m, int xSize, int ySize)
  {

    (*m) = new short int*[xSize];

    for (int i = 0 ; i < xSize; i++)
      {
        (*m)[i] = new short int[ySize];

        for (int j = 0; j < ySize; j++)
          (*m)[i][j] = NO_OVERLAP;

      }
  }

  void SplicingEventMatrix::delete_matrix(short int ***m, int xSize, int ySize)
  {
    for (int i = 0 ; i < xSize; i++)
      {
        delete [] (*m)[i];
      }

    delete [] (*m);
  }

  void SplicingEventMatrix::display_matrix(short int ***m, int xSize, int ySize)
  {
    log4cpp::Category& root = log4cpp::Category::getRoot();

    std::ostringstream oss;

    oss << endl << " X ";

    for (int j = 0 ; j < ySize; j++)
    {
      oss << ((j < 10) ? " " : "") << j << " ";
    }

    oss << endl;


    for (int i = 0 ; i < xSize; i++)
      {
        oss << ((i < 10) ? " " : "") << i << " ";

        for (int j = 0; j < ySize; j++)
          if ((*m)[i][j] < 10) {
            oss << " " << (*m)[i][j] << " ";
          } else {
            oss << (*m)[i][j] << " ";
          }

	oss << endl;

      }

    root.info(oss.str());

  }

  /**
   * <p>This method will detect any intron isoform event from 2 introns.</p>
   * <p>Intron isoform events rules:<p>
   * <ul>
   * <li>2 introns overlaps</li>
   * <li>There is an upstream exon</li>
   * <li>There is a downstream exon</li>
   * <li>their respective features on 5p and 3p overlap</li>
   * </ul>
   * <p>Schematic representation</p>
   * [====]-------[====]
   * [====]---------[==]
   * there is at least one pair of exons (5p or/and 3p) where the donor site is identical
   */
  void SplicingEventMatrix::checkIntronIsoform(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index)
  {
    log4cpp::Category& root = log4cpp::Category::getRoot();
    root.info("checking intron isoform event...");

    list<shared_ptr < TranscriptFeature > > list1;
    list<shared_ptr < TranscriptFeature > > list2;

    list< int > sites1;
    list< int > sites2;

    /**
     * Pre-conditions:
     *  - there must be one exon before and one exon after
     *  - Exons must overlap
     */
    if ( ( t1Index - 1 ) >= 0 &&
        ( t1Index + 1 ) < xSize &&
        ( t2Index - 1 ) >= 0 &&
        ( t2Index + 1 ) < ySize &&
        overlapMatrix[t1Index - 1][t2Index - 1] != NO_OVERLAP &&
        overlapMatrix[t1Index + 1][t2Index + 1] != NO_OVERLAP) {

          // allocate on the heap
          SplicingEvent *event =
            new SplicingEvent(
                II_EVENT,
                min(f1->getStart(), f2->getStart()),
                max(f1->getEnd(), f2->getEnd()),
                f1->getChromosome(),
                f1->getStrand()
            );

          /*
           * It's important to have a clear definition
           * of an intron isoform.
           * We will store the flanking exons.
           * the longest intron will always be stored in A and
           * the shortest intron in B.
           */

          //list1.push_back(f1);
          list1.push_back(t1Features[t1Index - 1]);
          list1.push_back(t1Features[t1Index + 1]);

          //list2.push_back(f2);
          list2.push_back(t2Features[t2Index - 1]);
          list2.push_back(t2Features[t2Index + 1]);

          sites1.push_back(0);
          sites1.push_back(f1->getStart());
          sites1.push_back(f1->getEnd());

          sites2.push_back(0);
          sites2.push_back(f2->getStart());
          sites2.push_back(f2->getEnd());

          if (f1->getLength() > f2->getLength()) {

            event->setSetA(list1);
            event->setSetB(list2);

            event->setSitesA(sites1);
            event->setSitesB(sites2);

            event->addTranscriptPair(t1, t2);

          } else {

            event->setSetA(list2);
            event->setSetB(list1);

            event->setSitesA(sites2);
            event->setSitesB(sites1);

            event->addTranscriptPair(t2, t1);

          }

          event->setGene(t1->getGene());

          intronIsoformEvents.push_back(shared_ptr<SplicingEvent>(event));

          //event->outputAsGFF(cout, "Ensembl");
    }

  }

  void SplicingEventMatrix::checkAlternativeFirstLastExon(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index)
  {
    log4cpp::Category& root = log4cpp::Category::getRoot();
    root.info("checking AFE/ALE...");

    list<shared_ptr < TranscriptFeature > > list1;
    list<shared_ptr < TranscriptFeature > > list2;
    list<shared_ptr < TranscriptFeature > > constitutives;

    list< int > sites1;
    list< int > sites2;
    list< int > constitutiveSites;

    /*
     * AFE, alternative first exon: constitutive second exon and mutually exclusive first exon.
     * @todo extends to multiple alternative first exons with relaxed constraints.
     */


    if (t1Index == 2 && t2Index == 2 &&
        overlapMatrix[t1Index - 2][t2Index - 2] == NO_OVERLAP &&
        overlapMatrix[t1Index][t2Index] == ID5P_ID3P
    ) {

      root.infoStream() << "checking alternative first exon: " << t1Index << "," << t2Index << " " << overlapMatrix[t1Index][t2Index] << "," << overlapMatrix[t1Index - 2][t2Index - 2] << log4cpp::eol;

      // allocate a new splicing event on the heap
      SplicingEvent *event = new SplicingEvent(AFE_EVENT,
          min(t1Features[t1Index - 2]->getStart(), t2Features[t2Index - 2]->getStart()),
          max(t1Features[t1Index - 2]->getEnd(), t2Features[t2Index - 2]->getEnd()),
          f1->getChromosome(),
          f1->getStrand());

      // add the alternative first exon to the SetA and SetB
      // but check the coordinates before:

      list1.push_back(t1Features[t1Index - 2]);
      list2.push_back(t2Features[t2Index - 2]);

      sites1.push_back(1);
      sites1.push_back(t1Features[t1Index - 2]->getStart());
      sites1.push_back(t1Features[t1Index - 2]->getEnd());

      sites2.push_back(1);
      sites2.push_back(t2Features[t2Index - 2]->getStart());
      sites2.push_back(t2Features[t2Index - 2]->getEnd());

      if (t1Features[t1Index - 2]->getEnd() < t2Features[t2Index - 2]->getStart()) {

        event->setSetA(list1);
        event->setSetB(list2);
        event->setSitesA(sites1);
        event->setSitesB(sites2);
        event->addTranscriptPair(t1, t2);

      } else {

        event->setSetA(list2);
        event->setSetB(list1);
        event->setSitesA(sites2);
        event->setSitesB(sites1);
        event->addTranscriptPair(t2, t1);

      }

      // add the constitutive exons (check for identifier in Ensembl)
      constitutives.push_back(f1);
      constitutiveSites.push_back(1);
      constitutiveSites.push_back(f1->getStart());
      constitutiveSites.push_back(f1->getEnd());

      if (f1->getIdentifier().compare(f2->getIdentifier()) != 0) {

        constitutives.push_back(f2);

      }

      event->setConstitutiveExons(constitutives);
      event->setConstitutiveSites(constitutiveSites);


      event->setGene(t1->getGene());

      alternativeFirstExonEvents.push_back(shared_ptr<SplicingEvent>(event));
    }

    /*
     * ALE, alternative last exon: constitutive penultimate exon and mutually exclusive last exon.
     * @todo extends to multiple alternative first exons with relaxed constraints.
     */

    if (( t1Index + 3) == xSize && (t2Index + 3) == ySize &&
        overlapMatrix[t1Index][t2Index] == ID5P_ID3P &&
        overlapMatrix[t1Index + 2][t2Index + 2] == NO_OVERLAP
    ) {

      root.infoStream() << "checking alternative last exon: " << t1Index << "," << t2Index << " " << overlapMatrix[t1Index][t2Index] << "," << overlapMatrix[t1Index + 2][t2Index + 2] << log4cpp::eol;

      // allocate a new splicing event on the heap
      SplicingEvent *event = new SplicingEvent(ALE_EVENT,
          min(t1Features[t1Index + 2]->getStart(), t2Features[t2Index + 2]->getStart()),
          max(t1Features[t1Index + 2]->getEnd(), t2Features[t2Index + 2]->getEnd()),
          f1->getChromosome(),
          f1->getStrand());

      // add the alternative last exon to the SetA and SetB
      // but check the coordinates before:

      list1.push_back(t1Features[t1Index + 2]);
      list2.push_back(t2Features[t2Index + 2]);

      sites1.push_back(1);
      sites1.push_back(t1Features[t1Index + 2]->getStart());
      sites1.push_back(t1Features[t1Index + 2]->getEnd());

      sites2.push_back(1);
      sites2.push_back(t2Features[t2Index + 2]->getStart());
      sites2.push_back(t2Features[t2Index + 2]->getEnd());

      if (t1Features[t1Index + 2]->getEnd() < t2Features[t2Index + 2]->getStart()) {

        event->setSetA(list1);
        event->setSetB(list2);
        event->setSitesA(sites1);
        event->setSitesB(sites2);
        event->addTranscriptPair(t1, t2);

      } else {

        event->setSetA(list2);
        event->setSetB(list1);
        event->setSitesA(sites2);
        event->setSitesB(sites1);
        event->addTranscriptPair(t2, t1);

      }

      // add the constitutive exons (check for identifier in Ensembl)
      constitutives.push_back(f1);
      constitutiveSites.push_back(1);
      constitutiveSites.push_back(f1->getStart());
      constitutiveSites.push_back(f1->getEnd());

      if (f1->getIdentifier().compare(f2->getIdentifier()) != 0) {

        constitutives.push_back(f2);

      }


      event->setConstitutiveExons(constitutives);
      event->setConstitutiveSites(constitutiveSites);


      event->setGene(t1->getGene());

      alternativeLastExonEvents.push_back(shared_ptr<SplicingEvent>(event));
    }

  }

  /**
   * exon isoform events rules:
   * 2 exons overlap
   * their respective features on 5p and 3p overlap
   * there is at least one pair of exons (5p or/and 3p) where the donor site is identical
   */
  void SplicingEventMatrix::checkExonIsoform(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index, bool bRELAX)
  {
    log4cpp::Category& root = log4cpp::Category::getRoot();
    root.info("checking exon isoform event...");

    list<shared_ptr < TranscriptFeature > > list1;
    list<shared_ptr < TranscriptFeature > > list2;

    list<int> sites1;
    list<int> sites2;

    /**
     * alternative initiation: same exon ends but different exon starts.
     * It could be the same exon but it's good to track for Alternative TSS later
     */

    if (t1Index == 0 && t2Index == 0 && overlapMatrix[t1Index][t2Index] == DIFF5P_ID3P) {

      // allocate a new splicing event on the heap
      SplicingEvent *event = new SplicingEvent(AI_EVENT,
          min(f1->getStart(), f2->getStart()),
          max(f1->getEnd(), f2->getEnd()),
          f1->getChromosome(),
          f1->getStrand());

      list1.push_back(f1);

      sites1.push_back(1);
      sites1.push_back(f1->getStart());
      sites1.push_back(f1->getEnd());

      list1.push_back(f2);

      sites1.push_back(1);
      sites1.push_back(f2->getStart());
      sites1.push_back(f2->getEnd());

      event->setSetA(list1);

      event->setSitesA(sites1);

      event->addTranscriptPair(t1, t2);
      event->setGene(t1->getGene());

      alternativeInitiationEvents.push_back(shared_ptr<SplicingEvent>(event));
    }

    /**
     * alternative termination: same exon starts but different exon ends.
     * It could be the same exon but it's good to track for Alternative PolyA later
     */

    if (t1Index == xSize-1 && t2Index == ySize-1 && overlapMatrix[t1Index][t2Index] == ID5P_DIFF3P) {

      // allocate a new splicing event on the heap
      SplicingEvent *event = new SplicingEvent(AT_EVENT,
          min(f1->getStart(), f2->getStart()),
          max(f1->getEnd(), f2->getEnd()),
          f1->getChromosome(),
          f1->getStrand());

      list1.push_back(f1);

      sites1.push_back(1);
      sites1.push_back(f1->getStart());
      sites1.push_back(f1->getEnd());

      list1.push_back(f2);

      sites1.push_back(1);
      sites1.push_back(f2->getStart());
      sites1.push_back(f2->getEnd());

      event->setSetA(list1);

      event->setSitesA(sites1);

      event->addTranscriptPair(t1, t2);
      event->setGene(t1->getGene());

      alternativeTerminationEvents.push_back(shared_ptr<SplicingEvent>(event));
    }



    /**
     * Exon Isoform event
     * Pre-conditions:
     *  - same end : one intron before
     *  - same start: one intron after
     *  - intron overlaps
     *  For EXON ISOFORM EVENTS
     *  - identical 5p exon donor sites
     *  - identical 3p exon acceptor sites
     */


    bool a3ss = 
			(
			 f1->getStrand() == 1 &&
			 f1->getEnd() == f2->getEnd() && f1->getStart() != f2->getStart() &&
			 ( t1Index - 1 ) >= 0 && ( t2Index - 1 ) >= 0 && 
			 (t1Features[t1Index - 1]->getStart() == t2Features[t2Index - 1]->getStart() ||
				(bRELAX == 1 && (overlapMatrix[t1Index - 1][t2Index - 1] & OVERLAP) == OVERLAP))
			 ) ||
			(
			 f1->getStrand() == -1 &&
			 f1->getEnd() == f2->getEnd() && f1->getStart() != f2->getStart() &&
			 ( t1Index + 1 ) < xSize && ( t2Index + 1 ) < ySize && 
			 (t1Features[t1Index + 1]->getStart() == t2Features[t2Index + 1]->getStart() ||
			 (bRELAX == 1 && (overlapMatrix[t1Index + 1][t2Index + 1] & OVERLAP) == OVERLAP))
			 );

			bool a5ss = 
			(
			 f1->getStrand() == 1 &&
			 f1->getStart() == f2->getStart() && f1->getEnd() != f2->getEnd() &&
			 ( t1Index + 1 ) < xSize && ( t2Index + 1 ) < ySize && 
			 (t1Features[t1Index + 1]->getEnd() == t2Features[t2Index + 1]->getEnd() ||
			 (bRELAX == 1 && (overlapMatrix[t1Index + 1][t2Index + 1] & OVERLAP) == OVERLAP))
			 ) ||
			(f1->getStrand() == -1 &&
			 f1->getStart() == f2->getStart() && f1->getEnd() != f2->getEnd() &&
			 ( t1Index - 1 ) >= 0 && ( t2Index - 1 ) >= 0 && 
			 (t1Features[t1Index - 1]->getEnd() == t2Features[t2Index - 1]->getEnd() ||
				(bRELAX == 1 && (overlapMatrix[t1Index - 1][t2Index - 1] & OVERLAP) == OVERLAP))
			 );

    bool exonIsoform =
      (f1->getStart() != f2->getStart() && f1->getEnd() != f2->getEnd() &&
      ( t1Index - 1 ) >= 0 && ( t2Index - 1 ) >= 0 && overlapMatrix[t1Index - 1][t2Index - 1] != NO_OVERLAP &&
      ( t1Index + 1 ) < xSize && ( t2Index + 1 ) < ySize && overlapMatrix[t1Index + 1][t2Index + 1] != NO_OVERLAP);


    int eventType = NO_EVENT;

    /**
     * Ok, check now that we have an Exon Isoform Event.
     * But we can refine the definition:
     * If the 3' splice site vary but not the 5' then we have a competing 3' splice site
     * If the 5' splice site vary but not the 3' then we have a competing 5' splice site
     * If both vary then just let it as an Exon Isoform Event if and only if the flanking exon
     * have the same splice sites.
     * watch out the strand
     */

    if (a3ss) {

      eventType = A3SS_EVENT;


    } else if (a5ss) {

      eventType = A5SS_EVENT;

    } else if (exonIsoform) {

      // check the intron boundary
       bool checkPreviousIntron =
             ((f1->getStrand() == 1) ? (t1Features[t1Index - 1]->getStart()
                 == t2Features[t2Index - 1]->getStart()) : (t1Features[t1Index
                 - 1]->getEnd() == t2Features[t2Index - 1]->getEnd()));


      bool checkNextIntron =
              ((f1->getStrand() == -1) ?
                  (t1Features[t1Index + 1]->getStart() == t2Features[t2Index + 1]->getStart()) :
                    (t1Features[t1Index + 1]->getEnd() == t2Features[t2Index + 1]->getEnd()));

            if (checkPreviousIntron && checkNextIntron) {
              eventType = EI_EVENT;
            }

    }

    if (eventType != NO_EVENT) {

    // allocate a new splicing event on the heap
    SplicingEvent *event = new SplicingEvent(eventType,
              min(f1->getStart(), f2->getStart()),
              max(f1->getEnd(), f2->getEnd()),
              f1->getChromosome(),
              f1->getStrand());

          /*
           * It's important to have a clear definition
           * of exon isoforms. We consider here deletion of a part
           * of one exon as a specific event. Thus, the longest exon
           * will always be stored in A and the shortest exon in B
           *
           */

          if (f1->getLength() > f2->getLength()) {

            list1.push_back(f1);
            list2.push_back(f2);
            event->setSetA(list1);
            event->setSetB(list2);
            event->addTranscriptPair(t1, t2);

            sites1.push_back(1);
            sites1.push_back(f1->getStart());
            sites1.push_back(f1->getEnd());

            sites2.push_back(1);
            sites2.push_back(f2->getStart());
            sites2.push_back(f2->getEnd());

            event->setSitesA(sites1);
            event->setSitesB(sites2);

          } else {

            list1.push_back(f2);
            list2.push_back(f1);
            event->setSetA(list1);
            event->setSetB(list2);
            event->addTranscriptPair(t2, t1);

            sites1.push_back(1);
            sites1.push_back(f2->getStart());
            sites1.push_back(f2->getEnd());

            sites2.push_back(1);
            sites2.push_back(f1->getStart());
            sites2.push_back(f1->getEnd());

            event->setSitesA(sites1);
            event->setSitesB(sites2);

          }

          event->setGene(t1->getGene());

          exonIsoformEvents.push_back(shared_ptr<SplicingEvent>(event));

          //event->outputAsGFF(cout, "Ensembl");
    }
  }



  /*
   * Compute intron retention events.
   * intron retention event
   * seek overlapping flanking exons on 5p and 3p
   *
   */
  void SplicingEventMatrix::checkIntronRetention(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index)
    {
      log4cpp::Category& root = log4cpp::Category::getRoot();
      root.info("checking intron retention event...");

      list<shared_ptr < TranscriptFeature > > list1;
      list<shared_ptr < TranscriptFeature > > list2;

      list< int > sites1;
      list< int > sites2;

      if (f1->getType() == EXON_TYPE) {

        //cout << " f1 is an exon:" << (( t2Index - 1 ) >= 0) << " " << (( t2Index + 1 ) < ySize) << endl;
        //cout << "t1Index:" << t1Index << " t2Index:" << t2Index << endl;

        /*
         * Look at the corresponding flanking exons on the other patterns.
         * We have to look in each direction
         * that is t2Index -1 and t2Index + 1 exists
         * and possibly more overlaps
         * [===================================]
         *   [===]----[==]--------[====]----[====]
         */

        if (( t2Index - 1 ) >= 0 &&
            ( t2Index + 1 ) < ySize &&
            (overlapMatrix[t1Index][t2Index - 1] & OVERLAP) == OVERLAP  &&
            (overlapMatrix[t1Index][t2Index + 1] & OVERLAP) == OVERLAP) {

              int t2Pos = 1; // points to the 5' exon
              list<shared_ptr < TranscriptFeature > > exonFeatures; // five prime features

              /* look at all the 5' exons and loop while there is an overlap */
              while (( t2Index - t2Pos ) >= 0 && (overlapMatrix[t1Index][t2Index - t2Pos] & OVERLAP) == OVERLAP) {

                // cerr << "[-]t2Index:" << t2Index << ", t2Pos:" << t2Pos << endl;

                exonFeatures.push_front(t2Features[t2Index - t2Pos]);
                t2Pos += 2;

              }

              t2Pos = 1; // points to the 3' exon
              while (( t2Index + t2Pos ) < ySize && (overlapMatrix[t1Index][t2Index + t2Pos] & OVERLAP) == OVERLAP) {

                // cerr << "[+]t2Index:" << t2Index << ", t2Pos:" << t2Pos << ", ySize:" << ySize << endl;
                // cerr << (t2Index + t2Pos) << endl;
                ////cerr << t2Features[t2Index - t2Pos]->getIdentifier() << ":"
                // cerr << t2Features[t2Index + t2Pos]->getStart() << "-" << t2Features[t2Index + t2Pos]->getEnd() << endl;

                exonFeatures.push_back(t2Features[t2Index + t2Pos]);
                t2Pos += 2;
              }

              // cerr << "-=-==-=--=-=-=-=-=-=-==-" << endl;

              SplicingEvent *event =
                new SplicingEvent(
                    IR_EVENT,
                    min(f1->getStart(), exonFeatures.front()->getStart()),
                    max(f1->getEnd(), exonFeatures.back()->getEnd()),
                    f1->getChromosome(),
                    f1->getStrand()
                );

              /*
               * Store the exons
               */

              list< shared_ptr< TranscriptFeature > >::const_iterator featureIterator;

              for(featureIterator=exonFeatures.begin(); featureIterator!=exonFeatures.end(); featureIterator++) {

                sites2.push_back(1);
                sites2.push_back((**featureIterator).getStart());
                sites2.push_back((**featureIterator).getEnd());

              }

              /*
               * Longest intron
               */

              list1.push_back(f1);
              event->setSetA(list1);

              sites1.push_back(1);
              sites1.push_back(f1->getStart());
              sites1.push_back(f1->getEnd());

              event->setSitesA(sites1);

              event->setSetB(exonFeatures);
              event->setSitesB(sites2);

              event->addTranscriptPair(t1, t2);
              event->setGene(t1->getGene());

              addNewEvent(shared_ptr<SplicingEvent>(event));
              //intronRetentionEvents.push_back(shared_ptr<SplicingEvent>(event));

              //event->outputAsGFF(cout, "Ensembl");

        }

      } else {

        /*
         * Look at the corresponding flanking exons on the other patterns.
         * that is t1Index -1 and t1Index + 1 exists
         *
         */

        //cout << ( t1Index - 1 ) << " " << ( t1Index + 1 ) << " ";
        //cout << overlapMatrix[t1Index - 1][t2Index] << " " << (overlapMatrix[t1Index - 1][t2Index] & ID5P) << " - ";
        //cout << overlapMatrix[t1Index + 1][t2Index] << " " << (overlapMatrix[t1Index + 1][t2Index] & ID3P) << endl;

        if (( t1Index - 1 ) >= 0 &&
            ( t1Index + 1 ) < xSize &&
            (overlapMatrix[t1Index - 1][t2Index] & OVERLAP) == OVERLAP  &&
            (overlapMatrix[t1Index + 1][t2Index] & OVERLAP) == OVERLAP) {

              int t1Pos = 1; // points to the 5' exon
              list<shared_ptr < TranscriptFeature > > exonFeatures; // five prime features

              /* look at all the 5' exons and loop while there is an overlap */
              while (( t1Index - t1Pos ) >= 0 && (overlapMatrix[t1Index - t1Pos][t2Index] & OVERLAP) == OVERLAP) {

                //cerr << "[-]t1Index:" << t1Index << ", t1Pos:" << t1Pos << endl;

                exonFeatures.push_front(t1Features[t1Index - t1Pos]);
                t1Pos += 2;

              }

              t1Pos = 1; // points to the 3' exon
              while (( t1Index + t1Pos ) < xSize && (overlapMatrix[t1Index + t1Pos][t2Index] & OVERLAP) == OVERLAP) {

                //cerr << "[+]t1Index:" << t1Index << ", t1Pos:" << t1Pos << ", xSize:" << xSize << endl;
                //cerr << (t1Index + t1Pos) << endl;
                ////cerr << t2Features[t2Index - t2Pos]->getIdentifier() << ":"
                //cerr << t1Features[t1Index + t1Pos]->getStart() << "-" << t1Features[t1Index + t1Pos]->getEnd() << endl;

                exonFeatures.push_back(t1Features[t1Index + t1Pos]);
                t1Pos += 2;
              }

              //cerr << "-----------------------------" << endl;

              SplicingEvent *event =
                new SplicingEvent(
                    IR_EVENT,
                    min(f2->getStart(), exonFeatures.front()->getStart()),
                    max(f2->getEnd(), exonFeatures.back()->getEnd()),
                    f2->getChromosome(),
                    f2->getStrand());

              /*
               * Store the exons
               */

              list< shared_ptr< TranscriptFeature > >::const_iterator featureIterator;

              for(featureIterator=exonFeatures.begin(); featureIterator!=exonFeatures.end(); featureIterator++) {

                sites1.push_back(1);
                sites1.push_back((**featureIterator).getStart());
                sites1.push_back((**featureIterator).getEnd());

              }

              /*
               * Longest intron
               */

              list2.push_back(f2);
              event->setSetA(list2);

              sites2.push_back(1);
              sites2.push_back(f2->getStart());
              sites2.push_back(f2->getEnd());

              event->setSitesA(sites2);

              event->setSetB(exonFeatures);
              event->setSitesB(sites1);

              event->addTranscriptPair(t2, t1);
              event->setGene(t2->getGene());

              addNewEvent(shared_ptr<SplicingEvent>(event));
              //intronRetentionEvents.push_back(shared_ptr<SplicingEvent>(event));

              //event->outputAsGFF(cout, "Ensembl");
        }
      }
    }

  void SplicingEventMatrix::checkCassetteExon(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index, bool bRELAX)
  {
    log4cpp::Category& root = log4cpp::Category::getRoot();
    root.info("checking cassette exon event...");

    list<shared_ptr < TranscriptFeature > > cassetteExons;

    /*
     * Scan 5p region for cassette exons
     */

    bool b5p = false;
    bool b3p = false;

    int index5p = 0;
    bool seek5p = false;

    if (f1->getType() == EXON_TYPE) {

      index5p = t1Index;

      do {

        /*
         * we look upstream to find other exons covered by the same intron.
         * if the current feature is an intron and has the same 5P end as the intron this is ok
         */

        index5p--;
        seek5p = index5p >= 0 && (overlapMatrix[index5p][t2Index] == PART_OF || overlapMatrix[index5p][t2Index] == ID5P_DIFF3P);

        if (seek5p && t1Features[index5p]->getType() == EXON_TYPE) {
          cassetteExons.push_front(t1Features[index5p]);
        }

      } while (seek5p);

      // push the exon in the vector
      cassetteExons.push_back(f1);

      /*
       * condition for 5p: no more features or overlapping exons with identical 3p ends
       * if we have reached the end of the feature and all the exons are inside the intron
       * we can't really decide what is happening. We take a conservative approach, uncomment
       * and replace the following line if you want to widen the type of events.
       */

      //b5p = (index5p < 0 || (overlapMatrix[index5p][t2Index - 1] & ID3P == ID3P));
      b5p = (index5p >= 0 && 
						 (
							(bRELAX == 1 && (overlapMatrix[index5p][t2Index - 1] & OVERLAP) == OVERLAP) ||
							(overlapMatrix[index5p][t2Index - 1] & ID3P) == ID3P)
						 );
			
    } else {

      index5p = t2Index;
      do {
        // look upstream
        index5p--;
        seek5p = index5p >= 0 && (overlapMatrix[t1Index][index5p] == PART_OF || overlapMatrix[t1Index][index5p] == ID5P_DIFF3P);
        if (seek5p && t2Features[index5p]->getType() == EXON_TYPE ) {
          cassetteExons.push_front(t2Features[index5p]);
        }

      } while (seek5p);

      // push the exon in the vector
      cassetteExons.push_back(f2);

      /*
       * condition for 5p: no more features or overlapping exons with identical 3p ends
			 * if we relax the conditions (bRELAX == 1), we accept that flanking exons overlaps
			 * but don't have necessarly the same 3p ends.
			 * Here, we still require a flanking exon.
       */
      b5p = (index5p >= 0 && 
						 ( 
							(bRELAX == 1 && (overlapMatrix[t1Index - 1][index5p] & OVERLAP) == OVERLAP) || 
							(overlapMatrix[t1Index - 1][index5p] & ID3P) == ID3P)
						 );

    }

    if (b5p) {

      /*
       * Scan 3p region
       */
      int index3p = 0;

      bool seek3p = false;

      if (f1->getType() == EXON_TYPE) {

        index3p = t1Index;

        do {
          // look upstream
          index3p++;
          seek3p = index3p < xSize && (overlapMatrix[index3p][t2Index] == PART_OF || overlapMatrix[index3p][t2Index] == DIFF5P_ID3P);
          if (seek3p && t1Features[index3p]->getType() == EXON_TYPE) {
            cassetteExons.push_back(t1Features[index3p]);
          }
        } while (seek3p);


        /*
         * condition for 3p: no more features or overlapping exons with identical 3p ends
         */
        //b3p = (index3p >= xSize || (overlapMatrix[index3p][t2Index + 1] & ID5P == ID5P));
        b3p = (index3p < xSize && 
							 ((bRELAX == 1 && (overlapMatrix[index3p][t2Index + 1] & OVERLAP) == OVERLAP) ||
								(overlapMatrix[index3p][t2Index + 1] & ID5P) == ID5P)
							 );
				
      } else {

        index3p = t2Index;

        do {
          // look upstream
          index3p++;
          seek3p = index3p < ySize && (overlapMatrix[t1Index][index3p] == PART_OF || overlapMatrix[t1Index][index3p] == DIFF5P_ID3P);
          if (seek3p && t2Features[index3p]->getType() == EXON_TYPE) {
            cassetteExons.push_back(t2Features[index3p]);
          }

        } while (seek3p);

        /*
         * condition for 3p: no more features or overlapping exons with identical 3p ends
				 * Again, we can relax the constraints with the bRELAX option to accept overlapping 
				 * exons on the 3p end.
         */
        //b3p = (index3p >= ySize || (overlapMatrix[t1Index + 1][index3p] & ID5P == ID5P));
        b3p = (index3p < ySize && 
							 ((bRELAX == 1 && (overlapMatrix[t1Index + 1][index3p] & OVERLAP) == OVERLAP) ||
								(overlapMatrix[t1Index + 1][index3p] & ID5P) == ID5P)
							 );

      }
    }

    /*
     * Cassette exon conditions fulfilled
     */
    //cout << "----" << endl;
    if (b3p && b5p) {

      shared_ptr<TranscriptFeature> intron = (f1->getType() == INTRON_TYPE) ? f1 : f2;

      list<int> sites1;
      //list<int> sites2;

      SplicingEvent *event =
        new SplicingEvent(
            CE_EVENT,
            intron->getStart(),
            intron->getEnd(),
            t1->getChromosome(),
            t1->getStrand()
        );

      //list< shared_ptr< TranscriptFeature > > intronList;
      //intronList.push_back(intron);
      event->setSetA(cassetteExons);
      //event->setSetB(intronList);

      list< shared_ptr< TranscriptFeature > >::const_iterator featureIterator;
      for(featureIterator=cassetteExons.begin(); featureIterator!=cassetteExons.end(); featureIterator++) {
          sites1.push_back(1);
          sites1.push_back((**featureIterator).getStart());
          sites1.push_back((**featureIterator).getEnd());
      }

      //sites2.push_back(0);
      //sites2.push_back(intron->getStart());
      //sites2.push_back(intron->getEnd());

      event->setSitesA(sites1);
      //event->setSitesB(sites2);

      if ((f1->getType() == INTRON_TYPE)) {

        event->addTranscriptPair(t1, t2);

      } else {

        event->addTranscriptPair(t2, t1);

      }
      event->setGene(t1->getGene());

      addNewEvent(shared_ptr<SplicingEvent>(event));

    }

    //cout << "++++" << endl;
  }

  void SplicingEventMatrix::checkMutualExclusion(const shared_ptr<TranscriptFeature> & f1, const shared_ptr<TranscriptFeature> & f2, int t1Index, int t2Index, bool bRELAX)
  {

    /*
     * Scan 5p region for mutually exclusive exons
     */



      log4cpp::Category& root = log4cpp::Category::getRoot();
      root.info("check mutually exclusive event...");

      list<shared_ptr < TranscriptFeature > > mutualExonT1;
      list<shared_ptr < TranscriptFeature > > mutualExonT2;

      list<int> sites1;
      list<int> sites2;

      int index5pT1;
      int index5pT2;

      bool b5pExclusive = false;
      bool b3pExclusive = false;

      bool b5pSeek = true;
      bool b3pSeek = true;

      int b5pPosition = 0;
      int b3pPosition = 0;


      /*
       * look for the next 5P exon that overlap with another exon
       */


      index5pT1 = t1Index;
      index5pT2 = t2Index;


      while (index5pT1 >= 0 && index5pT2 >= 0 && b5pSeek) {

        /*
         * look at the current exon
         */

        if (t1Features[index5pT1]->getType() == INTRON_TYPE) {

          index5pT1--; // look upstream for the next exon

        } else
          if (t1Features[index5pT1]->getType() == EXON_TYPE) {

            /*
             * is this exon it part of the current intron ?
             */
            root.infoStream() << "[b5pSeek] index5pT1: " << index5pT1 << " "<< t1Features[index5pT1]->getIdentifier() << " <=> " << t2Features[index5pT2]->getIdentifier() << log4cpp::eol;

            int t1Start = computeFeatureStart(t1Features[index5pT1]->getStart(), t1Features[index5pT1]->getEnd()) ;
            int t2Start = computeFeatureStart(t2Features[index5pT2]->getStart(), t2Features[index5pT2]->getEnd()) ;

            if (overlapMatrix[index5pT1][index5pT2] == PART_OF && t2Features[index5pT2]->getType() == INTRON_TYPE && t2Start < t1Start) {

              root.infoStream() << "[b5pSeek] insert " <<  t1Features[index5pT1]->getIdentifier() << " as mutual exon" << log4cpp::eol;
              mutualExonT1.push_front(t1Features[index5pT1]);
              index5pT1--; //then look upstream

            } else {

              /* No, this exon is not part of the current intron
               * so  move t2 pointer upstream for an overlap instead
               * and check that we still point onto the pattern (index5pT2 >= 0)
               */

              index5pT2--; // look upstream on the second transcript

              do {

                //cout << "index5pT2:" << index5pT2 << endl;

                /**
                 * The exon overlap with the previous feature.
                 */
                if ((overlapMatrix[index5pT1][index5pT2] & OVERLAP) == OVERLAP) {

                  //cout << "there is an overlap" << endl;

                  /*
                   * found overlap with another exon (we have to stop)
                   */

                  if (t2Features[index5pT2]->getType() == EXON_TYPE) {

                    //cout << "t2 feature is an exon ";

                    /*
                     * EXON <=> EXON overlap = STOP CONDITION
                     */
                    if ((overlapMatrix[index5pT1][index5pT2] & ID3P) == ID3P || bRELAX == 1) {

                      //cout << "with same 3P" << endl;

                      b5pSeek = false;
                      b5pExclusive = true;
                      b5pPosition = (strand == 1) ? t2Features[index5pT2]->getEnd() : t2Features[index5pT2]->getStart();

                    } else {

                      //cout << "with different 3P" << endl;
                      b5pSeek = false;
                      b5pExclusive = false;

                    }
                  } else {

                    /*
                     * EXON <=> INTRON overlap = GO BACK TO MAIN LOOP TO TREAT THIS CASE
                     * we make sure that it's always an intron in the main loop
                     */
                    break;

                  }
                } else {

                  /*
                   * NO OVERLAP
                   */

                  int t1End = computeFeatureEnd(t1Features[index5pT1]->getStart(), t1Features[index5pT1]->getEnd()) ;
                  int t2Start = computeFeatureStart(t2Features[index5pT2]->getStart(), t2Features[index5pT2]->getEnd()) ;

                  if (t2Features[index5pT2]->getType() == EXON_TYPE && t2Start > t1End) {

                    mutualExonT2.push_front(t2Features[index5pT2]);
                  }
                } // no overlap (so this exon is downstream and is mutual)

                index5pT2--;

              } while (index5pT2 >= 0 && b5pSeek);
            }
          } // t1 feature is an exon

      } // seek 5p



      /*
       * look for the next 5P exon that overlap with another exon
       */

      //cout << "now look for the next 5P exon that overlap with another exon" << endl;

      index5pT1 = t1Index;
      index5pT2 = t2Index;


      while (index5pT1 < xSize && index5pT2 < ySize && b3pSeek) {

        /*
         * look at the current exon
         */

        if (t1Features[index5pT1]->getType() == INTRON_TYPE) {

          index5pT1++;

        } else
          if (t1Features[index5pT1]->getType() == EXON_TYPE) {

            /*
             * is this exon it part of the current intron ?
             */

            root.infoStream() << "[b3pSeek] index5pT1: " << index5pT1 << " " << t1Features[index5pT1]->getIdentifier() << " <=> " << t2Features[index5pT2]->getIdentifier() << log4cpp::eol;


            int t1Start = computeFeatureStart(t1Features[index5pT1]->getStart(), t1Features[index5pT1]->getEnd()) ;
            int t2Start = computeFeatureStart(t2Features[index5pT2]->getStart(), t2Features[index5pT2]->getEnd()) ;

            if ( overlapMatrix[index5pT1][index5pT2] == PART_OF && t2Start < t1Start ) {

              /* we don't want to store the same exon twice */
              if (t1Features[index5pT1]->getIdentifier().compare(f1->getIdentifier()) != 0) {
                mutualExonT1.push_back(t1Features[index5pT1]);
              }
              index5pT1++;

            } else {

              /* No, the t2 intron feature does not cover t1 exon feature
               * so instead move t2 pointer downstream for an overlap
               * we can move downstream because t2 feature is an intron
               * the next feature will be an exon.
               */

              index5pT2++; // move downstream to the next exon.

              do {

                root.infoStream() << "scan forward index5pT2:" << index5pT2 << " " << t2Features[index5pT2]->getIdentifier() << log4cpp::eol;
                root.infoStream() << overlapMatrix[index5pT1][index5pT2] << log4cpp::eol;
                int o = (overlapMatrix[index5pT1][index5pT2] & OVERLAP);
		root.infoStream() << o << log4cpp::eol;

                /*
                 * If there is an overlap between the 2 exons.
                 */

                if ((overlapMatrix[index5pT1][index5pT2] & OVERLAP) == OVERLAP) {

                  root.infoStream() << "overlap between " << t1Features[index5pT1]->getIdentifier() << " " << t2Features[index5pT2]->getIdentifier() << log4cpp::eol;
                  /*
                   * found overlap
                   */

                  if (t2Features[index5pT2]->getType() == EXON_TYPE) {

                    /*
                     * EXON <=> EXON overlap = STOP CONDITION
                     */
                    if ((overlapMatrix[index5pT1][index5pT2] & ID5P) == ID5P || bRELAX == 1) {

                      b3pSeek = false;
                      b3pExclusive = true;
                      // at which position ?
                      b3pPosition = (strand == 1) ? t2Features[index5pT2]->getStart() : t2Features[index5pT2]->getEnd();

                    } else {

                      b3pSeek = false;
                      b3pExclusive = false;

                    }
                  } else {

                    /*
                     * EXON <=> INTRON overlap = GO BACK TO MAIN LOOP TO TREAT THIS CASE
                     * we make sure that it's always an intron in the main loop
                     */
                    break;

                  }
                } else {

                  /*
                   * THERE IS NO OVERLAP BETWEEN THE 2 FEATURES
                   * WHICH MEAN THAT t2 FEATURE IS DOWNSTREAM t1 exon
                   */

                  int t1Start = computeFeatureStart(t1Features[index5pT1]->getStart(), t1Features[index5pT1]->getEnd()) ;
                  int t2End = computeFeatureEnd(t2Features[index5pT2]->getStart(), t2Features[index5pT2]->getEnd()) ;

                  if (t2Features[index5pT2]->getType() == EXON_TYPE && t2End < t1Start) {

                    mutualExonT2.push_back(t2Features[index5pT2]);
                  }

                } // no overlap (so this exon is downstream and is mutual)

                index5pT2++;

              } while (index5pT2 < ySize && b3pSeek);
            }
          } // t1 feature is an exon

      } // seek 3p

      root.infoStream() << "check criteria" << log4cpp::eol;
      // actually if one of the set is empty, it could be treated as a cassette exon event.
      if (b5pExclusive && b3pExclusive && mutualExonT1.size() > 0 && mutualExonT2.size() > 0) {

        // create the event
        SplicingEvent *event =
          new SplicingEvent(
              MXE_EVENT,
              ((strand == 1) ? b5pPosition : b3pPosition),
              ((strand == 1) ? b3pPosition : b5pPosition),
              f1->getChromosome(),
              f1->getStrand()
          );

        // store the splice sites
        list< shared_ptr< TranscriptFeature > >::const_iterator featureIterator;

        for(featureIterator=mutualExonT1.begin(); featureIterator!=mutualExonT1.end(); featureIterator++) {
            sites1.push_back(1);
            sites1.push_back((**featureIterator).getStart());
            sites1.push_back((**featureIterator).getEnd());
        }

        for(featureIterator=mutualExonT2.begin(); featureIterator!=mutualExonT2.end(); featureIterator++) {
            sites2.push_back(1);
            sites2.push_back((**featureIterator).getStart());
            sites2.push_back((**featureIterator).getEnd());
        }

        /*
         * before setting the mutual exon, we have to decide
         * in which order we want to store them: we don't want
         * ambiguity later when we'll parse the gff file
         */

        shared_ptr < TranscriptFeature > e1 = mutualExonT1.front();
        shared_ptr < TranscriptFeature > e2 = mutualExonT2.front();

        if (e1->getStart() < e2->getStart()) {

          event->setSetA(mutualExonT1); // copy
          event->setSetB(mutualExonT2); // copy

          event->setSitesA(sites1);
          event->setSitesB(sites2);

          event->addTranscriptPair(t1, t2);

        } else {

          event->setSetA(mutualExonT2); // copy
          event->setSetB(mutualExonT1); // copy

          event->setSitesA(sites2);
          event->setSitesB(sites1);

          event->addTranscriptPair(t2, t1);

        }

        event->setGene(t1->getGene());

        addNewEvent(shared_ptr<SplicingEvent>(event));

      }

     // cout << "end cassette exon event..." << endl;
  }


}
