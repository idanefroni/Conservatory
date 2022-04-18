//
// File: PopGenMafSttistics.h
// Authors: Julien Dutheil
// Created: Tue Sep 29 2015
//

/*
Copyright or © or Copr. Julien Y. Dutheil, (2019)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _POPGENMAFSTATISTICS_H_
#define _POPGENMAFSTATISTICS_H_

#include <Bpp/Seq/Io/Maf/MafStatistics.h>

//From the STL:
#include <map>
#include <string>

namespace bpp {

/**
 * @brief Compute Fst value between two populations.
 *
 * Hudson's 92 estimator is used.
 */
class FstMafStatistics:
  public AbstractMafStatistics
{

  private:
    std::vector<std::string> pop1_, pop2_;
    unsigned int minNbPermutations_;
    unsigned int maxNbPermutations_;
    bool verbose_;

  public:
    FstMafStatistics(std::vector<std::string>& pop1, std::vector<std::string>& pop2,
        unsigned int minNbPermutations = 0,
        unsigned int maxNbPermutations = 0,
        bool verbose = true):
      AbstractMafStatistics(), pop1_(pop1), pop2_(pop2),
      minNbPermutations_(minNbPermutations),
      maxNbPermutations_(maxNbPermutations),
      verbose_(verbose)
    {}

    virtual ~FstMafStatistics() {}

  public:
    std::string getShortName() const { return "FstStatistics"; }
    std::string getFullName() const { return "Fst statistics."; }
    void compute(const MafBlock& block);
    std::vector<std::string> getSupportedTags() const;
};


} // end of namespace bpp

#endif //_POPGENMAFSTATISTICS_H_

