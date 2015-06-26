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

/** @file MemoryPool.hh
 * Declaration and implementation of a simple memory pool for the alignment algorithms.
 *
 * @author Dominik Dahlem
 */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifndef __MEMORYPOOL_HH__
#define __MEMORYPOOL_HH__

#ifndef __STDC_CONSTANT_MACROS
# define __STDC_CONSTANT_MACROS
#endif /* __STDC_CONSTANT_MACROS */

#include <algorithm>
#include <vector>

#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>


namespace alignment
{

typedef std::vector<std::vector<double> > DMatrix;
typedef std::vector<std::vector<boost::int32_t> > IMatrix;

/** @class MemoryPool
 *
 * This class is just a simple wrapper around the dynamic programming
 * matrices of the alignment algorithms. This memory pool is passed by
 * the calling function to the alignments in order to avoid
 * unnecessary memory (de)allocation in batch mode.
 */
class MemoryPool
{
 public:
  MemoryPool() {}
  ~MemoryPool() {}

  void checkDimensions(boost::uint32_t p_rows, boost::uint32_t p_cols)
  {
    if (p_rows >= m_H.size()) {
      m_H.resize(p_rows * 2);
    }
    if (p_rows >= m_I_i.size()) {
      m_I_i.resize(p_rows * 2);
    }
    if (p_rows >= m_I_j.size()) {
      m_I_j.resize(p_rows * 2);
    }
    if (p_cols >= m_H[0].size()) {
      BOOST_FOREACH(std::vector<double> & row, m_H) {
        row.resize(p_cols * 2);
      }
    }
    if (p_cols >= m_I_i[0].size()) {
      BOOST_FOREACH(std::vector<boost::int32_t> & row, m_I_i) {
        row.resize(p_cols * 2);
      }
    }
    if (p_cols >= m_I_j[0].size()) {
      BOOST_FOREACH(std::vector<boost::int32_t> & row, m_I_j) {
        row.resize(p_cols * 2);
      }
    }
  }

  void reset(boost::uint32_t p_rows, boost::uint32_t p_cols)
  {
    checkDimensions(p_rows, p_cols);
    BOOST_FOREACH(std::vector<double> & row, m_H) {
      std::fill_n(row.begin(), p_cols, 0.0);
    }
    BOOST_FOREACH(std::vector<boost::int32_t> & row, m_I_i) {
      std::fill_n(row.begin(), p_cols, 0.0);
    }
    BOOST_FOREACH(std::vector<boost::int32_t> & row, m_I_j) {
      std::fill_n(row.begin(), p_cols, 0.0);
    }
  }

  DMatrix & H()
  {
    return m_H;
  }

  IMatrix & I_i()
  {
    return m_I_i;
  }

  IMatrix & I_j()
  {
    return m_I_j;
  }

 private:
  DMatrix m_H;
  IMatrix m_I_i;
  IMatrix m_I_j;
};


}


#endif
