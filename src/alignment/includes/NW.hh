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

/** @file NW.hh
 * Declaration and implementation of the Needleman-Wunsch algorithm for global alignment.
 *
 * @author Dominik Dahlem
 */
#if HAVE_CONFIG_H
# include <config.h>
#endif

/**
 * taken and modified from
 * https://wiki.uni-koeln.de/biologicalphysics/index.php/Implementation_of_the_Smith-Waterman_local_alignment_algorithm
 */

#ifndef __NW_HH__
#define __NW_HH__

#ifndef __STDC_CONSTANT_MACROS
# define __STDC_CONSTANT_MACROS
#endif /* __STDC_CONSTANT_MACROS */

#ifndef NDEBUG
# include <iostream>
#endif /* NDEBUG */

#include <algorithm>
#include <string>
#include <vector>

#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>

#include "SimilarityAlgorithm.hh"
#include "Types.hh"


namespace alignment
{


class NW : public SimilarityAlgorithm
{
 private:
  bool m_justscores;

 public:
  NW(bool p_justscores) : SimilarityAlgorithm(), m_justscores(p_justscores) {}
  ~NW() {}

  alignmentResult align(
      common::StringVec & seq_a, /* sequence 1 */
      common::StringVec & seq_b, /* sequence 2 */
      AbstractDistanceMeasure & scoring_matrix,
      MemoryPool & mem) { /* scoring scheme */
#ifndef NDEBUG
    std::cout << "NW::operator()" << std::endl;
#endif /* NDEBUG */

    boost::uint32_t N_a = seq_a.size();
    boost::uint32_t N_b = seq_b.size();

    // initialize H
    /* in this case, we only initialize row 0 and col 0 */
    mem.reset(N_a + 1, N_b + 1);
    for (boost::int32_t i = 1; i < N_a+1; ++i) mem.H()[i][0] = -i * scoring_matrix.getDelta();
    for (boost::int32_t j = 1; j < N_b+1; ++j) mem.H()[0][j] = -j * scoring_matrix.getDelta();

    double temp[3];
    double *mdit;

    for (boost::uint32_t i = 1; i <= N_a; i++) {
      for (boost::uint32_t j = 1; j <= N_b; j++) {
        // calculate all possible paths to improve on prev optimal alignments
        temp[0] = mem.H()[i-1][j-1] + scoring_matrix.d(seq_a[i-1], seq_b[j-1]);
        temp[1] = mem.H()[i-1][j] - scoring_matrix.getDelta();
        temp[2] = mem.H()[i][j-1] - scoring_matrix.getDelta();

        // get max
        mdit = std::max_element(temp, temp+3);
        mem.H()[i][j] = *mdit;

        switch(std::distance(temp, mdit)) {
          case 0: // score in (i,j) stems from a match/mismatch
            mem.I_i()[i][j] = i-1;
            mem.I_j()[i][j] = j-1;
            break;
          case 1: // score in (i,j) stems from a del in sequence A
            mem.I_i()[i][j] = i-1;
            mem.I_j()[i][j] = j;
            break;
          case 2: // score in (i,j) stems from a del in sequence B
            mem.I_i()[i][j] = i;
            mem.I_j()[i][j] = j-1;
            break;
        }
      }
    }
#ifndef NDEBUG
    std::cout << "H" << std::endl;
    for (boost::uint32_t i = 0; i <= N_a; ++i) {
      for (boost::uint32_t j = 0; j <= N_b; ++j) {
        std::cout << mem.H()[i][j] << ",";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
#endif

    // store results
    alignmentResult result;
    result.score = mem.H()[N_a][N_b];
    result.alignment.resize(2);

    if (!m_justscores) {
      /* we now backtrack from the bottom right cell of H */
      boost::int32_t current_i = N_a, current_j = N_b;
      boost::int32_t next_i = mem.I_i()[current_i][current_j];
      boost::int32_t next_j = mem.I_j()[current_i][current_j];
      boost::int32_t tick = 0;

      common::StringVec consensus_a, consensus_b;
      consensus_a.resize(N_a + N_b + 2);
      consensus_b.resize(N_a + N_b + 2);

      /* we have to go from MN to 00 */
      while (current_i != 0 || current_j != 0) {
        if (next_i == current_i) { consensus_a[tick] = "-"; } // deletion in A
        else { consensus_a[tick] = scoring_matrix.match(seq_a[current_i-1], seq_b[current_j-1]); }      // match/mismatch in A

        if (next_j == current_j) { consensus_b[tick] = "-"; } // deletion in B
        else { consensus_b[tick] = scoring_matrix.match(seq_a[current_i-1], seq_b[current_j-1]); }      // match/mismatch in B

        current_i = next_i;
        current_j = next_j;
        next_i = mem.I_i()[current_i][current_j];
        next_j = mem.I_j()[current_i][current_j];
        tick++;
      }

      for(boost::int32_t i = tick-1; i >= 0; i--) { result.alignment[0].push_back(consensus_a[i]); }
      for(boost::int32_t j = tick-1; j >= 0; j--) { result.alignment[1].push_back(consensus_b[j]); }
    }

    return result;
  }
};


}


#endif
