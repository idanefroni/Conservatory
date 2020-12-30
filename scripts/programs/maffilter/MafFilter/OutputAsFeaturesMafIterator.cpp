//
// File: OutputasFeaturesMafIterator.cpp
// Authors: Julien Dutheil
// Created: Fri Mar 22 2013
//

/*
Copyright or © or Copr. Julien Y. Dutheil, (2013)

This file is part of MafFilter.

MafFilter is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MafFilter is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MafFilter.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "OutputAsFeaturesMafIterator.h"

using namespace bpp;

//From the STL:
#include <string>
#include <vector>

using namespace std;

void OutputAsFeaturesMafIterator::writeHeader(std::ostream& out) const
{
  out << "##Maf2Annotations" << endl;
  out << "Chr\tStart\tStop\tStrand\tSrcSize" << endl;
  //We probably here ultimately want to write in GFF or sthg like this?
}

void OutputAsFeaturesMafIterator::writeBlock(std::ostream& out, const MafBlock& block) const
{
  if (block.hasSequenceForSpecies(species_)) {
    vector<const MafSequence*> sequences = block.getSequencesForSpecies(species_);
    for (size_t i = 0; i < sequences.size(); ++i) {
      out << sequences[i]->getChromosome() << "\t" << sequences[i]->start() << "\t" << sequences[i]->stop() << "\t" << sequences[i]->getStrand() << "\t" << sequences[i]->getSrcSize() << endl; 
    }
  }
}

