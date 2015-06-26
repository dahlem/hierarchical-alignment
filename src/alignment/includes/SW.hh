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

/** @file SW.hh
 * Declaration and implementation of the Smith-Waterman algorithm for local alignment.
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

#ifndef __SW_HH__
#define __SW_HH__

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

#include "SimilarityAlgorithm.hh"
#include "Types.hh"


namespace alignment
{


class SW : public SimilarityAlgorithm
{
 private:
  bool m_justscores;

 public:
  SW(bool p_justscores) : SimilarityAlgorithm(), m_justscores(p_justscores) {}
  ~SW() {}

  alignmentResult align(
      common::StringVec &seq_a, /* sequence 1 */
      common::StringVec &seq_b, /* sequence 2 */
      AbstractDistanceMeasure & scoring_matrix,
      MemoryPool &mem) { /* gap penalty */

    boost::uint32_t N_a = seq_a.size();
    boost::uint32_t N_b = seq_b.size();

    // initialize H
    mem.reset(N_a + 1, N_b + 1);

    double temp[4];
    double *mdit;

    double H_max = 0.;
    boost::uint32_t i_max = 0, j_max = 0;

    for (boost::uint32_t i = 1; i <= N_a; i++) {
      for (boost::uint32_t j = 1; j <= N_b; j++) {
        // calculate all possible paths to improve on prev optimal alignments
        // we assume for 1 and 2 a linear gap scoring scheme. Hence we do not need to apply a max operator
        //  over previous indices since we know that the max value is attained at the previous index
        temp[0] = mem.H()[i-1][j-1] + scoring_matrix.d(seq_a[i-1], seq_b[j-1]);
        temp[1] = mem.H()[i-1][j] - scoring_matrix.getDelta();
        temp[2] = mem.H()[i][j-1] - scoring_matrix.getDelta();
        temp[3] = 0.0;

        // get max
        mdit = std::max_element(temp, temp+4);
        mem.H()[i][j] = *mdit;

        switch (std::distance(temp, mdit)) {
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
          case 3: // (i,j) is the beginning of a subsequence
            mem.I_i()[i][j] = i;
            mem.I_j()[i][j] = j;
            break;
        }

        // store maximum
        if(mem.H()[i][j] > H_max){
          H_max = mem.H()[i][j];
          i_max = i;
          j_max = j;
        }
      }
    }

    // store results
    alignmentResult result;
    result.score = H_max;
    result.alignment.resize(2);

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

    if (!m_justscores) {
      // Backtracking from H_max
      boost::int32_t current_i = i_max, current_j = j_max;
      boost::int32_t next_i = mem.I_i()[current_i][current_j];
      boost::int32_t next_j = mem.I_j()[current_i][current_j];
      boost::int32_t tick = 0;

      common::StringVec consensus_a, consensus_b;
      consensus_a.resize(N_a + N_b + 2);
      consensus_b.resize(N_a + N_b + 2);

      while (((current_i != next_i) || (current_j != next_j)) && (next_j >= 0) && (next_i >= 0) && (current_i > 0) && (current_j > 0)) {
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

      for (boost::int32_t i = tick-1; i >= 0; i--) { result.alignment[0].push_back(consensus_a[i]); }
      for (boost::int32_t j = tick-1; j >= 0; j--) { result.alignment[1].push_back(consensus_b[j]); }
    }

    return result;
  } // sw
};


}


#endif



