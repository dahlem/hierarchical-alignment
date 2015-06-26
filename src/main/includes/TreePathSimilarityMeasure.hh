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

/** @file TreePathSimilarityMeasure
 * Declaration of the tree-based path similarity scoring scheme for alignments.
 *
 * @author Dominik Dahlem
 */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifndef __TREEPATHSIMILARITYMEASURE_HH__
#define __TREEPATHSIMILARITYMEASURE_HH__

#ifndef __STDC_CONSTANT_MACROS
# define __STDC_CONSTANT_MACROS
#endif /* __STDC_CONSTANT_MACROS */

#ifndef NDEBUG
# include <iostream>
#endif /* NDEBUG */

#include <boost/cstdint.hpp>

#include "AbstractDistanceMeasure.hh"
#include "Types.hh"


class TreePathSimilarityMeasure : public alignment::AbstractDistanceMeasure
{
 public:
  TreePathSimilarityMeasure(double p_delta, common::IntVec &p_levels, common::StringIntMap &p_pos,
                            common::StrStrMap &p_lcas)
      : alignment::AbstractDistanceMeasure(p_delta), m_levels(p_levels), m_pos(p_pos), m_lcas(p_lcas) {}

  ~TreePathSimilarityMeasure() {}

  double d(common::Symbol &a, common::Symbol &b)
  {
    if (a == b) { return 1.0; }

    auto left = m_pos[a] < m_pos[b] ? a : b;
    auto right = m_pos[b] >= m_pos[a] ? b : a;
    auto key = std::make_tuple(left, right);

    if (m_lcas.find(key) == m_lcas.end()) { return -1.0; }

    common::Symbol lca = m_lcas[key];
    boost::uint32_t levelLCA = m_levels[m_pos[lca]];
    boost::uint32_t distLeftLCA = m_levels[m_pos[left]] - m_levels[m_pos[lca]];
    boost::uint32_t distRightLCA = m_levels[m_pos[right]] - m_levels[m_pos[lca]];

    return (1.0 + levelLCA)/(1.0 + levelLCA + distLeftLCA + distRightLCA);
  }

  common::Symbol match(common::Symbol &a, common::Symbol &b)
  {
    if (a == b) { return a; }

    auto left = m_pos[a] < m_pos[b] ? a : b;
    auto right = m_pos[b] >= m_pos[a] ? b : a;
    common::Symbol lca = m_lcas[std::make_tuple(left, right)];

#ifndef NDEBUG
    std::cout << "match: " << left << ", " << right << ", " << lca << std::endl;
#endif

    return lca;
  }

 private:
  common::IntVec &m_levels;
  common::StringIntMap &m_pos;
  common::StrStrMap &m_lcas;
};


#endif
