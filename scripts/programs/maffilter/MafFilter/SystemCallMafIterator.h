//
// File: SystemCallMafIterator.h
// Authors: Julien Dutheil
// Created: Tue Sep 29 2015
//

/*
Copyright or © or Copr. Julien Y. Dutheil, (2015)

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

#ifndef _SYSTEMCALLMAFITERATOR_H_
#define _SYSTEMCALLMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/MafIterator.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/Io/OSequence.h>

//From the STL:
#include <iostream>
#include <string>
#include <memory>

namespace bpp {

/**
 * @brief This iterator calls an external program on each block.
 */
class SystemCallMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::unique_ptr<OAlignment> alnWriter_;
    std::string inputFile_;
    std::unique_ptr<IAlignment> alnReader_;
    std::string outputFile_;
    std::string call_;
    bool hotTest_;

  public:
    SystemCallMafIterator(
        MafIterator* iterator,
        OAlignment* alnWriter,
        const std::string& inputFile,
        IAlignment* alnReader,
        const std::string& outputFile,
        const std::string& callCmd,
        bool hotTest = false) :
      AbstractFilterMafIterator(iterator),
          alnWriter_(alnWriter),
          inputFile_(inputFile),
          alnReader_(alnReader),
          outputFile_(outputFile),
          call_(callCmd),
          hotTest_(hotTest)
    {}

  private:
    SystemCallMafIterator(const SystemCallMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      alnWriter_(),
      inputFile_(iterator.inputFile_),
      alnReader_(),
      outputFile_(iterator.outputFile_),
      call_(iterator.call_),
      hotTest_(iterator.hotTest_)
    {}
    
    SystemCallMafIterator& operator=(const SystemCallMafIterator& iterator)
    {
      inputFile_ = iterator.inputFile_;
      outputFile_ = iterator.outputFile_;
      call_ = iterator.call_;
      hotTest_ = iterator.hotTest_;
      return *this;
    }


  public:
    MafBlock* analyseCurrentBlock_();

};

} // namespace bpp.

#endif //_SYSTEMCALLMAFITERATOR_H_

