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

/** @file AbstractDistanceMeasure.hh
 * Declaration of the abstract distance measure to be implemented by concrete scoring schemes.
 *
 * @author Dominik Dahlem
 */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifndef __ABSTRACTDIST_HH__
#define __ABSTRACTDIST_HH__

#ifndef __STDC_CONSTANT_MACROS
# define __STDC_CONSTANT_MACROS
#endif /* __STDC_CONSTANT_MACROS */

#include "Types.hh"

namespace alignment
{

/** @class AbstractDistanceMeasure
 *
 * This class declares the abstract interface for alignment scoring schemes.
 */
class AbstractDistanceMeasure
{
 protected:
  double m_delta;

 public:
  AbstractDistanceMeasure(double p_delta) : m_delta(p_delta) {}

  virtual ~AbstractDistanceMeasure() {}

  virtual double d(common::Symbol &a, common::Symbol &b) = 0;

  virtual common::Symbol match(common::Symbol &a, common::Symbol &b) = 0;

  inline
  double getDelta()
  {
    return m_delta;
  }
};


}


#endif
