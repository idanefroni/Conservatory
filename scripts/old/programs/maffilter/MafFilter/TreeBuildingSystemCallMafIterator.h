//
// File: TreeBuildingSystemCallMafIterator.h
// Authors: Julien Dutheil
// Created: Sat Jun 18 2016
//

/*
Copyright or © or Copr. Julien Y. Dutheil, (2016)

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

#ifndef _TREEBUILDINGSYSTEMCALLMAFITERATOR_H_
#define _TREEBUILDINGSYSTEMCALLMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/MafIterator.h>
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Seq/Io/OSequence.h>

//From the STL:
#include <iostream>
#include <string>
#include <memory>

namespace bpp {

/**
 * @brief This iterator calls an external program on each block.
 */
class TreeBuildingSystemCallMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::unique_ptr<OAlignment> alnWriter_;
    std::string inputFile_;
    std::unique_ptr<ITree> treeReader_;
    std::string outputFile_;
    std::string call_;
    std::string propertyName_;

  public:
    TreeBuildingSystemCallMafIterator(
        MafIterator* iterator,
        OAlignment* alnWriter,
        const std::string& inputFile,
        ITree* treeReader,
        const std::string& outputFile,
        const std::string& callCmd,
        const std::string& propertyName) :
      AbstractFilterMafIterator(iterator),
          alnWriter_(alnWriter),
          inputFile_(inputFile),
          treeReader_(treeReader),
          outputFile_(outputFile),
          call_(callCmd),
          propertyName_(propertyName)
    {}

  private:
    TreeBuildingSystemCallMafIterator(const TreeBuildingSystemCallMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      alnWriter_(),
      inputFile_(iterator.inputFile_),
      treeReader_(),
      outputFile_(iterator.outputFile_),
      call_(iterator.call_),
      propertyName_(iterator.propertyName_)
    {}
    
    TreeBuildingSystemCallMafIterator& operator=(const TreeBuildingSystemCallMafIterator& iterator)
    {
      inputFile_ = iterator.inputFile_;
      outputFile_ = iterator.outputFile_;
      call_ = iterator.call_;
      propertyName_ = iterator.propertyName_;
      return *this;
    }


  public:
    MafBlock* analyseCurrentBlock_();

};

} // namespace bpp.

#endif //_TREEBUILDINGSYSTEMCALLMAFITERATOR_H_

