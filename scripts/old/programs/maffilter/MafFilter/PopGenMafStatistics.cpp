//
// File: PopGenMafSttistics.cpp
// Authors: Julien Dutheil
// Created: Wed Jun 26 2019
//

/*
Copyright or © or Copr. Julien Y. Dutheil, (2019)

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

#include "PopGenMafStatistics.h"

#include <Bpp/PopGen/PolymorphismSequenceContainer.h>
#include <Bpp/PopGen/SequenceStatistics.h>
#include <Bpp/App/ApplicationTools.h>

using namespace std;
using namespace bpp;

vector<string> FstMafStatistics::FstMafStatistics::getSupportedTags() const
{
  vector<string> tags;
  tags.push_back("Fst");
  if (minNbPermutations_ > 0 || maxNbPermutations_ > 0) {
    tags.push_back("Fst.p-value");
    tags.push_back("Fst.nb-perm");
  }
  return tags;
}

void FstMafStatistics::compute(const MafBlock& block)
{
  PolymorphismSequenceContainer poly(block.getAlignment());
  vector<string> names = block.getSpeciesList();
  poly.setSequencesNames(names); //Have to be unique, no paralog allowed, otherwise exception thrown.
  for (auto it = pop1_.begin(); it != pop1_.end(); ++it) {
    poly.setGroupId(*it, 1);
  }
  for (auto it = pop2_.begin(); it != pop2_.end(); ++it) {
    poly.setGroupId(*it, 2);
  }
  double fst = SequenceStatistics::fstHudson92(poly, 1, 2);
  result_.setValue("Fst", fst);
  if (minNbPermutations_ > 0 || maxNbPermutations_ > 0) {
    double nbTests = 0;
    double nbPermutations = 0;
    while(nbPermutations < minNbPermutations_ || (nbTests == 0 && nbPermutations < maxNbPermutations_)) {
      ++nbPermutations;
      if (verbose_) {
        ApplicationTools::displayGauge(nbPermutations, maxNbPermutations_, '=', "Compute Fst on permutations");
      }
      vector<string> individuals = pop1_;
      individuals.insert(individuals.end(), pop2_.begin(), pop2_.end());
      random_shuffle(individuals.begin(), individuals.end());
      for (size_t j = 0; j < pop1_.size(); ++j) {
        poly.setGroupId(individuals[j], 1);
      }
      for (size_t j = pop1_.size(); j < pop1_.size() + pop2_.size(); ++j) {
        poly.setGroupId(individuals[j], 2);
      }
      double randFst = SequenceStatistics::fstHudson92(poly, 1, 2);
      if (randFst >= fst) ++nbTests;
    }
    double pv = (nbTests + 1) / static_cast<double>(nbPermutations);
    result_.setValue("Fst.p-value", pv);
    result_.setValue("Fst.nb-perm", nbPermutations);
  }
}


