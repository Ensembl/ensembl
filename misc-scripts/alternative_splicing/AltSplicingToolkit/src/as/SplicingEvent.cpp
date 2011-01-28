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

#include <string>
#include <algorithm>
#include <sstream>

#include "SplicingEvent.h"

using namespace as;
using namespace std;

SplicingEvent::SplicingEvent() : GeneFeature::GeneFeature()
{
}

SplicingEvent::SplicingEvent(int type, unsigned int start, unsigned int end, string chr, short int strand) : GeneFeature::GeneFeature(type, start, end, chr, strand)
{
}

SplicingEvent::~SplicingEvent()
{
}

void SplicingEvent::setSetA(const list< shared_ptr< TranscriptFeature > > &a)
{
  setA.clear(); setA.assign(a.begin(), a.end());
}
void SplicingEvent::setSetB(const list< shared_ptr< TranscriptFeature > > &b)
{
  setB.clear(); setB.assign(b.begin(), b.end());
}

void SplicingEvent::setSitesA(const list< int > &a)
{
  sitesA.clear(); sitesA.assign(a.begin(), a.end());
}
void SplicingEvent::setSitesB(const list< int > &b)
{
  sitesB.clear(); sitesB.assign(b.begin(), b.end());
}
void SplicingEvent::setConstitutiveSites(const list< int > &c)
{
  constitutiveSites.clear(); constitutiveSites.assign(c.begin(), c.end());
}

void SplicingEvent::setConstitutiveExons(const list< shared_ptr< TranscriptFeature > > &constitutives)
{
  constitutiveExons.clear(); constitutiveExons.assign(constitutives.begin(), constitutives.end());
}

/**
 * To avoid the cost of copying the list, we return
 * the list by reference rather than by value.
 */
const list< shared_ptr< TranscriptFeature > >& SplicingEvent::getSetA() const
{
  return setA;
}

/**
 * To avoid the cost of copying the list, we return
 * the list by reference rather than by value.
 */
const list< shared_ptr< TranscriptFeature > >& SplicingEvent::getSetB() const
{
  return setB;
}

const list< int >& SplicingEvent::getSitesA() const
{
  return sitesA;
}

const list< int >& SplicingEvent::getSitesB() const
{
  return sitesB;
}

const list< int >& SplicingEvent::getConstitutiveSites() const
{
  return constitutiveSites;
}

/**
 * To avoid the cost of copying the list, we return
 * the list by reference rather than by value.
 */
const list< shared_ptr< TranscriptFeature > >& SplicingEvent::getConstitutiveExons() const
{
  return constitutiveExons;
}

/**
 * Allow to compare two multisets of transcript features.
 *
 * <p>Two sets are said equals if they contains the same features with the
 * same coordinates. If one set has more features than the other set of features
 * then the coordinates of the extra features must be the same as the one in the
 * other sets. This can happen because some exons may have the same coordinates
 * but with a different identifier.</p>
 *
 * <p>Implementation: check if the coordinates of the features in set1 falls into
 * set2. Loop on set1 features and compare with the features of set2</p>
 *
 * @param set1 The first set of features.
 * @param set2 The second set of features.
 * @return true if the set contains identical exonic features, false otherwhise.
 *
 * @see SplicingEvent#equals for an usage of this method
 */
bool SplicingEvent::compareFeatures(const list< shared_ptr< TranscriptFeature > > &set1, const list< shared_ptr< TranscriptFeature > > &set2 ) const
{
  bool same = true;

  if (set1.size() > 0) {

    list< shared_ptr< TranscriptFeature > >::const_iterator mii1 = set1.begin();
    list< shared_ptr< TranscriptFeature > >::const_iterator mii2 = set2.begin();

    while (same) {

      /* don't forget the double indirection to get the actual data */
      same =  ((**mii1).getType() == (**mii2).getType() && (**mii1).getStart() == (**mii2).getStart() && (**mii1).getEnd() == (**mii2).getEnd());

      mii2++;
      mii1++;

      if (mii1 == set1.end() || mii2 == set2.end()) {
        break;
      }
    }
  }

  return same;
}

bool SplicingEvent::compareFeatureSites(const list< int > &sites1, const list< int > &sites2 ) const
{
  bool same = false;

  if (sites1.size() == sites2.size() && sites1.size() % 3 == 0) {

    same = true;

    if (sites1.size() > 0) {

      list< int >::const_iterator mii1 = sites1.begin();
      list< int >::const_iterator mii2 = sites2.begin();

      // size % 3 = 0
      while (same) {

        int sameType = ((*mii1) == (*mii2));

        mii2++;
        mii1++;

        bool sameStart = ((*mii1) == (*mii2));

        mii2++;
        mii1++;

        bool sameEnd = ((*mii1) == (*mii2));

        mii2++;
        mii1++;

        same =  sameType && sameStart && sameEnd;

        if (mii1 == sites1.end() || mii2 == sites2.end()) {
          break;
        }
      }
    }
  }

  return same;

}

string SplicingEvent::getTypeAsString() const
{
  string eventType;
  switch (type) {

    case MXE_EVENT: eventType = "MXE";  break;
    case CE_EVENT: eventType = "CE";    break;

    case EI_EVENT: eventType = "EI";    break;
    case II_EVENT: eventType = "II";    break;

    case IR_EVENT: eventType = "IR";    break;

    case AI_EVENT: eventType = "AI";    break;
    case AT_EVENT: eventType = "AT";    break;

    case AFE_EVENT: eventType = "AFE";  break;
    case ALE_EVENT: eventType = "ALE";  break;

    case CNE_EVENT: eventType = "CNE";  break;

    default: eventType = "UNDEFINED";
  }
  return eventType;
}

bool SplicingEvent::equals(const SplicingEvent &event) const
{

  // 2 events are considered equals if:
  // * they have the same gene
  // * they have the same type
  // * they share the same coordinates
  // * they have the same number of related features (with the same coordinates)
  //
  // Special cases:
  // * Alternative Initiation: same first exon 3' end (we can merge all the events with this case and change the coordinates accordingly)
  // * Alternative Termination: same last exon 5' end (we can merge all the events with this case and change the coordinates accordingly)
  // * Alternative First exon: equivalent mutually exclusive first exon (we can merge all the events with this case and change the coordinates accordingly)
  // * Alternative Last exon: equivalent mutually exclusive last exon (we can merge all the events with this case and change the coordinates accordingly)

  bool result =
    (gene->getIdentifier().compare((event.getGene())->getIdentifier()) == 0) &&
    type == event.getType() &&
    (
        (start == event.getStart() && end == event.getEnd() &&
        event.getSitesA().size() == sitesA.size() && event.getSitesB().size() == sitesB.size()) ||

        /**
         * AI: Alternative Initiation
         *    [===]- SetA
         *     [==]- SetB
         */

        (type == AI_EVENT &&
         ((strand == 1 && end == event.getEnd()) ||
         (strand == -1 && start == event.getStart()))) ||

        /**
         * AT: Alternative Termination
         *    -[===] SetA
         *    -[==] SetB
         */

        (type == AT_EVENT &&
         ((strand == 1 && start == event.getStart()) ||
         (strand == -1 && end == event.getEnd()))) ||

        /**
         * AFE: Alternative First Exon: check the constitutive second exon
         * [==]-------[===] SetA
         *      [==]--[===] SetB
         */
        (type == AFE_EVENT &&
        setA.back()->getStart() == event.getSetA().back()->getStart() &&
        setA.back()->getEnd() == event.getSetA().back()->getEnd() &&
        (strand == 1 && end == event.getEnd() && setA.front()->getEnd() == event.getSetA().front()->getEnd() && setB.front()->getEnd() == event.getSetB().front()->getEnd()) ||
        (strand == -1 && start == event.getStart() && setA.front()->getStart() == event.getSetA().front()->getStart() && setB.front()->getStart() == event.getSetB().front()->getStart())) ||

        /**
         * ALE: Alternative Last Exon
         * [===]-------[==] SetA
         * [===]--[==]      SetB
         */
        (type == ALE_EVENT &&
         setA.back()->getStart() == event.getSetA().back()->getStart() &&
         setA.back()->getEnd() == event.getSetA().back()->getEnd() &&
         (strand == -1 && end == event.getEnd() && setA.front()->getEnd() == event.getSetA().front()->getEnd() && setB.front()->getEnd() == event.getSetB().front()->getEnd()) ||
         (strand == 1 && start == event.getStart() && setA.front()->getStart() == event.getSetA().front()->getStart() && setB.front()->getStart() == event.getSetB().front()->getStart()))

    );

/*    if (type == AFE_EVENT) {
      cout << "AFE_EVENT:" << (setA.back()->getStart() == event.getSetA().back()->getStart()) << endl;
      cout << "AFE_EVENT:" <<(setA.back()->getEnd() == event.getSetA().back()->getEnd()) << endl;
      if (strand == 1) {
        cout << "AFE_EVENT 1:" << (end == event.getEnd()) << endl;
      } else {
        cout << "AFE_EVENT -1:" << (start == event.getStart()) << endl;
        cout << "AFE_EVENT -1:" << (setA.front()->getStart() == event.getSetA().front()->getStart()) << endl;
        cout << "AFE_EVENT -1:" << (setB.front()->getStart() == event.getSetB().front()->getStart()) << endl;
      }
    } else if (type == AI_EVENT) {

      cout << "AI_EVENT coord:" <<
        ((strand == 1 && end == event.getEnd()) ||
        (strand == -1 && start == event.getStart())) << endl;
      cout << "AI_EVENT gene id:" << (gene->getIdentifier().compare((event.getGene())->getIdentifier()) == 0) << endl;
      cout << "AI_EVENT types:" << (type == event.getType()) << endl;

    }*/

  //if (type == IR_EVENT) {

    //cout << "***IR " << gene->getIdentifier() << " : " << result << endl;
    //cout << start << "<>" << event.getStart() << endl;
    //cout << end << "<>" << event.getEnd() << endl;

  //}
    /**
     * Alternative Initiation, Alternative Termination,
     * Alternative First Exon, Alternative Last Exon
     * in this case, we don't compare the features.
     */
    if ((type & ALT_FIRST_LAST_EXON_EVENT) == ALT_FIRST_LAST_EXON_EVENT) {
      //cout << "ALT_FIRST_LAST_EXON_EVENT: " << result << endl;
      return result;
    }



    return result &&
      (compareFeatureSites(sitesA, event.getSitesA()) && compareFeatureSites(sitesB, event.getSitesB())) ||
      (compareFeatureSites(sitesA, event.getSitesB()) && compareFeatureSites(sitesB, event.getSitesA()));

     // (compareFeatures(setA, event.getSetA()) && compareFeatures(setB, event.getSetB())) ||
     // (compareFeatures(setA, event.getSetB()) && compareFeatures(setB, event.getSetA()));

}

void SplicingEvent::outputAsGFF(ostream &oStream, string datasource, int number) const
{
  string eventType;
  string name;

  /**
    * Add more information specific to each type of events:
    */
  std::ostringstream oss;


  switch (type) {

    case MXE_EVENT: 

			eventType = "MXE"; 
			name = "mutual exclusion"; 
			break;

    case CE_EVENT:

      eventType = "CE"; 
			name = "cassette exon";
      oss << "NbCrypticExons=" << setA.size() << "; ";
      break;

    case II_EVENT: 

			eventType = "II"; 
			name = "intron isoform"; 
			break;

    case IR_EVENT:

      eventType = "IR"; 
			name = "intron retention";
      oss << "NbIntrons=" << (setB.size()-1) << "; ";
      break;

    case EI_EVENT: 

			eventType = "EI"; 
			name = "exon isoform";

			if (strand == 1) {

				oss << "3pModification=" << std::abs((int)(**setA.begin()).getStart() - (int)(**setB.begin()).getStart()) << "bp; ";
			} else {

				oss << "3pModification=" << std::abs((int)(**setA.begin()).getEnd() - (int)(**setB.begin()).getEnd()) << "bp; ";

			}

			if (strand == -1) {

				oss << "5pModification=" << std::abs((int)(**setA.begin()).getStart() - (int)(**setB.begin()).getStart()) << "bp; ";
			} else {

				oss << "5pModification=" << std::abs((int)(**setA.begin()).getEnd() - (int)(**setB.begin()).getEnd()) << "bp; ";

			}

			break;

    case A3SS_EVENT:

      eventType = "A3SS"; name = "alternative 3' splice site";

      if (strand == 1) {
        oss << "3pModification=" << std::abs((int)(**setA.begin()).getStart() - (int)(**setB.begin()).getStart()) << "bp; ";
      } else {
        oss << "3pModification=" << std::abs((int)(**setA.begin()).getEnd() - (int)(**setB.begin()).getEnd()) << "bp; ";
      }
      break;

    case A5SS_EVENT:

      eventType = "A5SS"; name = "alternative 5' splice site";

      if (strand == -1) {
        oss << "5pModification=" << std::abs((int)(**setA.begin()).getStart() - (int)(**setB.begin()).getStart()) << "bp; ";
      } else {
        oss << "5pModification=" << std::abs((int)(**setA.begin()).getEnd() - (int)(**setB.begin()).getEnd()) << "bp; ";
      }
      break;

    case AI_EVENT: eventType = "AI"; name = "alternative initiation"; break;

    case AT_EVENT: eventType = "AT"; name = "alternative termination"; break;

    case AFE_EVENT: eventType = "AFE"; name = "alternative first exon"; break;

    case ALE_EVENT: eventType = "ALE"; name = "alternative last exon"; break;

    case CNE_EVENT: eventType = "CNE"; name = "constitutive exon"; break;

    default: eventType = "UNDEFINED"; name="undefined";
  }

  oStream << chr << "\t" << datasource << "\t" << eventType << "\t" << start << "\t" << end << "\t.\t" << getStrandAsString()
  << "\t.\tID=" << gene->getIdentifier() << "-" << eventType << "-" << number << "; Derives_from=" << gene->getIdentifier() << "; Name=" << name << "; ";
  oStream << oss.str();

  // show the mutually exclusive exons:


  list< shared_ptr< TranscriptFeature > >::const_iterator featureIterator;
  int c = 0;

  if (setA.size() > 0) {

    // create a list of all the transcripts
    list< string > listA;
    list< string > listB;
    vector< pair<boost::shared_ptr< Transcript >, boost::shared_ptr< Transcript > > >::const_iterator pii;

    for( pii=transcriptPairs.begin(); pii!=transcriptPairs.end(); pii++ )
      {
	string transcriptIdA = pii->first->getIdentifier();
	listA.push_back(transcriptIdA);
	string transcriptIdB = pii->second->getIdentifier();
	listB.push_back(transcriptIdB);
      }

    if (type == AT_EVENT || type == AI_EVENT || type == CE_EVENT) {
      listA.merge(listB);
    }

    oStream << "FeaturesA=";
    for(featureIterator=setA.begin(); featureIterator!=setA.end(); featureIterator++) {

    oStream << ((c > 0) ? "," : "") <<  (**featureIterator).getIdentifier();
    oStream << "[" << findInvolvedTranscriptIdentifier(oStream, listA, (**featureIterator));
    oStream << "]";
    c++;
    }
    oStream << "; ";
  }

  c = 0;
  if (setB.size() > 0) {

    // create a list of all the transcripts
    list< string > listB;
    
    vector< pair<boost::shared_ptr< Transcript >, boost::shared_ptr< Transcript > > >::const_iterator pii;

    for( pii=transcriptPairs.begin(); pii!=transcriptPairs.end(); pii++ )
      {
	string transcriptId = pii->second->getIdentifier();
	listB.push_back(transcriptId);
      }
    
    oStream << "FeaturesB=";
    for(featureIterator=setB.begin(); featureIterator!=setB.end(); featureIterator++) {
      oStream << ((c > 0) ? "," : "") << (**featureIterator).getIdentifier() << "[" << findInvolvedTranscriptIdentifier(oStream, listB, (**featureIterator)) << "]";
	c++;
    }
    oStream << "; ";
  }

  /* Display 5'/3' sites */
  list< int >::const_iterator siteIterator;
  c = 0;
  if (sitesA.size() > 0) {
    oStream << "SitesA=";
    for(siteIterator=sitesA.begin(); siteIterator!=sitesA.end(); siteIterator++) {

      int co = (*siteIterator);
      string connector = (co == 0) ? "i" : "e";
      //cerr << connector << endl;
      siteIterator++;
      int startSite = (*siteIterator);
      //cerr << startSite << endl;
      siteIterator++;
      int endSite = (*siteIterator);
      //cerr << endSite << endl;

      oStream << ((c > 0) ? "," : "") << connector << "(" << startSite << "-" << endSite << ")";
      c++;

    }
  oStream << "; ";
  }

  c = 0;
  if (sitesB.size() > 0) {
    oStream << "SitesB=";
  for(siteIterator=sitesB.begin(); siteIterator!=sitesB.end(); siteIterator++) {

    int co = (*siteIterator);
    string connector = (co == 0) ? "i" : "e";
    //cerr << connector << endl;
    siteIterator++;
    int startSite = (*siteIterator);
    //cerr << startSite << endl;
    siteIterator++;
    int endSite = (*siteIterator);
    //cerr << endSite << endl;

    oStream << ((c > 0) ? "," : "") << connector << "(" << startSite << "-" << endSite << ")";
    c++;

  }
  oStream << "; ";
  }

  c = 0;
  if (constitutiveExons.size() > 0) {
    oStream << "ConstitutiveExons=";
  for(featureIterator=constitutiveExons.begin(); featureIterator!=constitutiveExons.end(); featureIterator++) {
    oStream << ((c > 0) ? "," : "") << (**featureIterator).getIdentifier();
    c++;
  }
  oStream << "; ";
  }

  c = 0;
  if (constitutiveSites.size() > 0) {
    oStream << "ConstitutiveSites=";
  for(siteIterator=constitutiveSites.begin(); siteIterator!=constitutiveSites.end(); siteIterator++) {

    int co = (*siteIterator);
    string connector = (co == 0) ? "i" : "e";
    //cerr << connector << endl;
    siteIterator++;
    int startSite = (*siteIterator);
    //cerr << startSite << endl;
    siteIterator++;
    int endSite = (*siteIterator);
    //cerr << endSite << endl;

    oStream << ((c > 0) ? "," : "") << connector << "(" << startSite << "-" << endSite << ")";
    c++;

  }
  oStream << "; ";
  }

  // now display the transcript pairs
  vector< pair<boost::shared_ptr< Transcript >, boost::shared_ptr< Transcript > > >::const_iterator pii;

  for( pii=transcriptPairs.begin(); pii!=transcriptPairs.end(); pii++ )
  {
    oStream << "Pair=" << pii->first->getIdentifier() << "," << pii->second->getIdentifier() << "; ";
  }




  oStream << endl;
}

string SplicingEvent::findInvolvedTranscriptIdentifier(ostream &oStream, list<string> &transcriptList, TranscriptFeature &feature) const 
{

  list<string>::const_iterator sii;
  list<string> inserted;

  int count = 0;

  std::ostringstream oss;

  vector< shared_ptr< Transcript > >::const_iterator vii;
  for ( vii=feature.getTranscripts().begin(); vii!=feature.getTranscripts().end(); vii++ ) {

    string currentId = (*vii)->getIdentifier();

    list<string>::const_iterator sii;

    for ( sii = transcriptList.begin(); sii!=transcriptList.end(); sii++ ) {

      //oStream << currentId << " compare to " << *sii << "\n";
 
      if (currentId == (*sii)) {


	bool alreadyInserted = false;
	
	list<string>::const_iterator sii2; 
	
	for ( sii2 = inserted.begin(); sii2!=inserted.end(); sii2++ ) {

	  if (currentId == (*sii2)) {
	    alreadyInserted = true;
	    break;
	  }
	  
	}

	if (!alreadyInserted) {
	  inserted.push_back(currentId);
	
	  if (count>0)
	    oss << ":";
	  
	  oss << currentId;
	  count++;
	}
      } 

    }

  }

  return  oss.str();

}

/**
 * Find a similar splicing event from the vector.
 */
vector< shared_ptr< SplicingEvent > >::iterator SplicingEvent::find(vector< shared_ptr< SplicingEvent > > &events) const
{

  vector< shared_ptr< SplicingEvent > >::iterator eventIt;

  for(eventIt=events.begin(); eventIt!=events.end(); eventIt++) {

    if (this->equals(**eventIt)) {
      break;
    }
  }

  return eventIt;
}


bool SplicingEvent::contains(vector< shared_ptr< SplicingEvent > > &events) const
{

  vector< shared_ptr< SplicingEvent > >::iterator eventIt = find(events);

  return (eventIt != events.end());
}

const vector< pair<boost::shared_ptr< Transcript >, boost::shared_ptr< Transcript > > > &SplicingEvent::getTranscriptPairs() const
{
  return transcriptPairs;
}

void SplicingEvent::mergeTranscriptPairs(const vector< pair<boost::shared_ptr< Transcript >, boost::shared_ptr< Transcript > > > &pairs)
{
  vector < pair<boost::shared_ptr< Transcript >, boost::shared_ptr< Transcript > > > newPairs( transcriptPairs.size() + pairs.size() );

  std::merge( transcriptPairs.begin(), transcriptPairs.end(), pairs.begin(), pairs.end(), newPairs.begin() );

  transcriptPairs.clear();
  transcriptPairs.assign(newPairs.begin(), newPairs.end());

}

void SplicingEvent::addTranscriptPair(const shared_ptr< Transcript > & t1, const shared_ptr< Transcript > & t2)
{
  shared_ptr< Transcript > p1(t1);
  shared_ptr< Transcript > p2(t2);
  pair< shared_ptr< Transcript >, shared_ptr< Transcript> > pr(p1, p2);
  transcriptPairs.push_back( pr );

}

void SplicingEvent::mergeCoordinates(unsigned int start, unsigned int end)
{
  // extends the coordinates
  this->start = min(this->start, start);
  this->end = max(this->end, end);
}

void SplicingEvent::mergeSetA(const list< shared_ptr< TranscriptFeature > > &a)
{
  mergeSets(setA, a);
}

void SplicingEvent::mergeSetB(const list< shared_ptr< TranscriptFeature > > &b)
{
  mergeSets(setB, b);
}

void SplicingEvent::mergeConstitutiveExons(const list< shared_ptr< TranscriptFeature > > &constitutives)
{
  mergeSets(constitutiveExons, constitutives);
}

void SplicingEvent::mergeSets(list< shared_ptr< TranscriptFeature > > &x, const list< shared_ptr< TranscriptFeature > > &y)
{

  // to keep
  list< shared_ptr< TranscriptFeature > > keep4X;
  list< shared_ptr< TranscriptFeature > >::const_iterator featureIteratorX;
  list< shared_ptr< TranscriptFeature > >::const_iterator featureIteratorY;

  for(featureIteratorY=y.begin(); featureIteratorY!=y.end(); featureIteratorY++) {

    bool sameID = false;

    for(featureIteratorX=x.begin(); featureIteratorX!=x.end(); featureIteratorX++) {

      if ((**featureIteratorY).getIdentifier().compare((**featureIteratorX).getIdentifier()) == 0) {
        sameID = true;
        break;
      }
    }

    if (!sameID) {
      keep4X.push_back((*featureIteratorY));
    }
  }

  x.merge(keep4X);

}

AlternativeInitiation::AlternativeInitiation()
{
}

AlternativeInitiation::~AlternativeInitiation()
{
}

AlternativeTermination::AlternativeTermination()
{
}

AlternativeTermination::~AlternativeTermination()
{
}
