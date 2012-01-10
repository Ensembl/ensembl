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

#include <string>

#include "Feature.h"

namespace as {

  Feature::Feature() : Coordinates()
  {
  }

  Feature::Feature(std::string featureIdentifier) : Coordinates()
  {
    identifier = featureIdentifier;
  }

  Feature::Feature(int type, unsigned int start, unsigned int end, std::string chr, short int strand) : Coordinates(start,end), type(type), chr(chr), strand(strand)
  {
  }

  Feature::~Feature(void)
  {
    //cout << "Calling Feature destructor for type " << type << " " << identifier << endl;
  }

  std::string Feature::getIdentifier() const
  {
    return identifier;
  }

  short int Feature::getStrand() const
  {
    return strand;
  }

  std::string Feature::getStrandAsString() const
  {
    return (strand == 1) ? "+" : "-";
  }

  std::string Feature::getChromosome() const
  {
    return chr;
  }

  void Feature::setStrand(short int v)
  {
    strand = v;
  }

  void Feature::setChromosome(std::string v)
   {
     chr = v;
   }

  void Feature::setType(int featureType)
  {
    type = featureType;
  }

  void Feature::setIndex(int v) {
    index = v;
  }

  int Feature::getIndex() const {
    return index;
  }

  int Feature::getType() const {
    return type;
  }

}
// ;P
