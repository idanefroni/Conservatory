//
// File: TreeBuildingSystemCallMafIterator.cpp
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

#include "TreeBuildingSystemCallMafIterator.h"

// From bpp-phyl
#include "Bpp/Phyl/TreeTemplate.h"

using namespace bpp;
using namespace std;

MafBlock* TreeBuildingSystemCallMafIterator::analyseCurrentBlock_() {
  currentBlock_ = iterator_->nextBlock();
  if (! currentBlock_)
    return 0;
  unique_ptr<AlignedSequenceContainer> aln(currentBlock_->getAlignment().clone());
  
  //We translate sequence names to avoid compatibility issues
  vector<string> names(aln->getNumberOfSequences());
  map<string, string> nameIndex;
  for (size_t i = 0; i < names.size(); ++i) {
    names[i] = "seq" + TextTools::toString(i);
    nameIndex[names[i]] = currentBlock_->getSequence(i).getName();
  }
  aln->setSequencesNames(names);
  
  //Write sequences to file:
  alnWriter_->writeAlignment(inputFile_, *aln, true);

  //Call the external program:
  int rc = system(call_.c_str());
  if (rc) throw Exception("TreeBuildingSystemCallMafIterator::analyseCurrentBlock_(). System call exited with non-zero status.");

  //Then read the generated tree and assign sequence names:
  unique_ptr< Tree > result(treeReader_->read(outputFile_));
  unique_ptr< TreeTemplate<Node> > tree(new TreeTemplate<Node>(*result));
  vector<Node*> leaves = tree->getLeaves();
  for (size_t i = 0; i < leaves.size(); ++i) {
    leaves[i]->setName(nameIndex[leaves[i]->getName()]);
  }
  currentBlock_->setProperty(propertyName_, tree.release());

  //Done:
  return currentBlock_;
}

