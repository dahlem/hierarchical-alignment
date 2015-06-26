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

/** @file CL.cc
 * Implementation of the command-line parsing of the main routine
 *
 * @author Dominik Dahlem
 */
#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifndef __STDC_CONSTANT_MACROS
# define __STDC_CONSTANT_MACROS
#endif /* __STDC_CONSTANT_MACROS */

#include <iostream>

#include <boost/cstdint.hpp>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;

#include "CL.hh"


CL::CL()
    : m_opt_desc(new po::options_description("Options"))
{
  // Declare the supported options.
  po::options_description opt_general("General Configuration");
  opt_general.add_options()
      (HELP.c_str(), "produce help message")
      (VERS.c_str(), "show the version")
      ;

  po::options_description opt_io("I/O Configuration");
  opt_io.add_options()
      (RESULTS_DIR.c_str(), po::value <std::string>()->default_value("./results"), "results directory.")
      (EULER_LEVELS.c_str(), po::value <std::string>()->default_value(""), "Filename of the vertex levels in the Euler Circuit.")
      (EULER_POSITIONS.c_str(), po::value <std::string>()->default_value(""), "Filename of the vertex positions in the Euler Circuit.")
      (LCA.c_str(), po::value <std::string>()->default_value(""), "Filename with the LCAs computed offline.")
      (SET_1.c_str(), po::value <std::string>()->default_value(""), "Filename of the source set.")
      (SET_2.c_str(), po::value <std::string>()->default_value(""), "Filename of the target set.")
      ;

  po::options_description opt_ha("Algorithm Configuration");
  opt_ha.add_options()
      (ALG.c_str(), po::value <boost::int32_t>()->default_value(1), "Algorithm: 1 - local alignment, 2 - global alignment.")
      (SCORES.c_str(), po::value <bool>()->default_value(0), "Compute just alignment scores, no backtracking.")
      (GAP_PENALTY.c_str(), po::value <double>()->default_value(1.33), "Gap penalty for the alignments.")
      ;

  m_opt_desc->add(opt_general);
  m_opt_desc->add(opt_io);
  m_opt_desc->add(opt_ha);
}


int CL::parse(int argc, char *argv[], args_t &p_args)
{
  po::variables_map vm;

  po::store(po::parse_command_line(argc, argv, (*m_opt_desc.get())), vm);
  po::notify(vm);

  if (vm.count(HELP)) {
    std::cout << (*m_opt_desc.get()) << std::endl;
    return EXIT_FAILURE;
  }

  if (vm.count(VERS)) {
    std::cout << PACKAGE_NAME << " " << PACKAGE_VERSION << std::endl;
    std::cout << argv[0] << std::endl;
    return EXIT_FAILURE;
  }

  if (vm.count(RESULTS_DIR.c_str())) {
    p_args.results_dir = vm[RESULTS_DIR.c_str()].as <std::string>();
  }

  if (vm.count(EULER_LEVELS.c_str())) {
    p_args.euler_levels = vm[EULER_LEVELS.c_str()].as <std::string>();
    if (!fs::exists(p_args.euler_levels)) {
      std::cerr << "The filename " << p_args.euler_levels << " with the Vertex levels in the Euler Circuit does not exist!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (vm.count(EULER_POSITIONS.c_str())) {
    p_args.euler_positions = vm[EULER_POSITIONS.c_str()].as <std::string>();
    if (!fs::exists(p_args.euler_levels)) {
      std::cerr << "The filename " << p_args.euler_positions << " with the Vertex positions in the Euler Circuit does not exist!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (vm.count(LCA.c_str())) {
    p_args.lca = vm[LCA.c_str()].as <std::string>();
    if (!fs::exists(p_args.lca)) {
      std::cerr << "The filename " << p_args.lca << " containing the LCAs does not exist!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (vm.count(SET_1.c_str())) {
    p_args.set_1 = vm[SET_1.c_str()].as <std::string>();
    if (!fs::exists(p_args.set_1)) {
      std::cerr << "The filename " << p_args.set_1 << " containing the source set does not exist!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (vm.count(SET_2.c_str())) {
    p_args.set_2 = vm[SET_2.c_str()].as <std::string>();
    if (!fs::exists(p_args.set_2)) {
      std::cerr << "The filename " << p_args.set_2 << " containing the target set does not exist!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (vm.count(ALG.c_str())) {
    p_args.alg = vm[ALG.c_str()].as <boost::int32_t>();
  }

  if (vm.count(SCORES.c_str())) {
    p_args.scores = vm[SCORES.c_str()].as <bool>();
  }

  if (vm.count(GAP_PENALTY.c_str())) {
    p_args.gap_penalty = vm[GAP_PENALTY.c_str()].as <double>();
  }

  std::cout << argv[0] << " " << PACKAGE_VERSION << std::endl;
  std::cout << PACKAGE_NAME << std::endl;
  std::cout << p_args << std::endl;

  return EXIT_SUCCESS;
}
