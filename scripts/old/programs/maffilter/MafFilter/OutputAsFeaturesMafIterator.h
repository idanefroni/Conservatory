//
// File: OutputasFeaturesMafIterator.h
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

#ifndef _OUTPUTASFEATURESMAFITERATOR_H_
#define _OUTPUTASFEATURESMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/MafIterator.h>

//From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp {

/**
 * @brief This iterator print as features the sequences for a given species.
 */
class OutputAsFeaturesMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::ostream* output_;
    std::string species_;

  public:
    OutputAsFeaturesMafIterator(MafIterator* iterator, std::ostream* out, const std::string& species) :
      AbstractFilterMafIterator(iterator), output_(out), species_(species)
    {
      if (output_)
        writeHeader(*output_);
    }

  private:
    OutputAsFeaturesMafIterator(const OutputAsFeaturesMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      output_(iterator.output_),
      species_(iterator.species_)
    {}
    
    OutputAsFeaturesMafIterator& operator=(const OutputAsFeaturesMafIterator& iterator)
    {
      output_ = iterator.output_;
      species_ = iterator.species_;
      return *this;
    }


  public:
    MafBlock* analyseCurrentBlock_() {
      currentBlock_ = iterator_->nextBlock();
      if (output_ && currentBlock_)
        writeBlock(*output_, *currentBlock_);
      return currentBlock_;
    }

  private:
    void writeHeader(std::ostream& out) const;
    void writeBlock(std::ostream& out, const MafBlock& block) const;
};

} //end of namespace bpp.

#endif //_OUTPUTASFEATURESMAFITERATOR_H_

