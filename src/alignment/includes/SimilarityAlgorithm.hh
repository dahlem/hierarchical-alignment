// Copyright (C) 2015 Dominik Dahlem <Dominik.Dahlem@gmail.com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

/** @file SimilarityAlgorithm.hh
 * Declaration of the interface to similarity computations
 *
 * @author Dominik Dahlem
 */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifndef __SIMILARITYALGORITH_HH__
#define __SIMILARITYALGORITH_HH__

#ifndef __STDC_CONSTANT_MACROS
# define __STDC_CONSTANT_MACROS
#endif /* __STDC_CONSTANT_MACROS */


#include "Types.hh"
#include "AbstractDistanceMeasure.hh"
#include "MemoryPool.hh"

namespace alignment
{


/** @struct alignmentResult
 * This structure keeps the result of an alignment, containing the score and a vector of alignments.
 */
struct alignmentResult {
  double score;
  std::vector<common::StringVec> alignment;
};


/** @class SimilarityAlgorithm
 * This class provides the wrapper interface to the similarity algorithms
 *
 * @author <a href="mailto:dominik.dahlem@gmail.com">Dominik Dahlem</a>
 */
class SimilarityAlgorithm
{
 public:
  SimilarityAlgorithm() {}
  virtual ~SimilarityAlgorithm() {}

  /** @fn alignmentResult align(std::vector<boost::flyweight<std::string> >, std::vector<boost::flyweight<std::string> >, AbstractDistanceMeasure &, MemoryPool &)
   *
   * This function computes similarities between two sequences and given a scoring scheme. The
   * concrete implementation needs to implement this definition in order to comply with this
   * interface.
   *
   * @param common::StringVec & the first sequence
   * @param common::StringVec & the second sequence
   * @param AbstractDistanceMeasure & the scoring scheme
   * @param MemoryPool & external memory holding the dynamic programming matrices
   */
  virtual alignmentResult align(common::StringVec & seq_a, common::StringVec & seq_b,
                                AbstractDistanceMeasure & scoring_matrix,
                                MemoryPool & mem) = 0;
};


}


#endif
