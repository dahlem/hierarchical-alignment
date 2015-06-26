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

#if HAVE_CONFIG_H
# include <config.h>
#endif

#ifndef __STDC_CONSTANT_MACROS
# define __STDC_CONSTANT_MACROS
#endif /* __STDC_CONSTANT_MACROS */

#ifndef NDEBUG
# include <iterator>
#endif /* NDEBUG */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include <boost/cstdint.hpp>
#include <boost/flyweight.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "CL.hh"
#include "Types.hh"
#include "TreePathSimilarityMeasure.hh"

#include "AbstractDistanceMeasure.hh"
#include "SimilarityAlgorithm.hh"
#include "NW.hh"
#include "SW.hh"


typedef boost::tokenizer <boost::escaped_list_separator <char> > Tokenizer;

static common::Symbol::initializer fw_symbol_init;


namespace std
{

std::ostream& operator<<(std::ostream &os, const common::StringIntMap::value_type &v) {
  return os << v.first << ": " << v.second;
}

std::ostream& operator<<(std::ostream &os, const common::StrStrMap::value_type &v) {
  return os << std::get<0>(v.first) << ", " << std::get<1>(v.first) << ": " << v.second;
}

std::ostream& operator<<(std::ostream &os, const common::StringVec &v) {
  std::copy(v.begin(), v.end(), std::ostream_iterator<common::Symbol>(os, "\t"));
  return os;
}

}

int main(int argc, char *argv[])
{
  args_t args;
  CL cl;

  if (cl.parse(argc, argv, args)) {
    return EXIT_SUCCESS;
  }

  // parse the euler levels
  std::ifstream eulerLevelsFile;
  eulerLevelsFile.open(args.euler_levels.c_str());
  if (!eulerLevelsFile.is_open()) {
    std::cerr << "Could not open file: " << args.euler_levels << std::endl;
    return EXIT_FAILURE;
  }

#ifndef NDEBUG
  std::cout << "Reading the Euler Levels..." << std::endl;
#endif /* NDEBUG */

  common::IntVec euler_levels;
  std::string line;
  while (!eulerLevelsFile.eof()) {
    std::getline(eulerLevelsFile, line);
#ifndef NDEBUG
    std::cout << "Read line: " << line << std::endl;
#endif /* NDEBUG */

    if (line != "") {
      boost::algorithm::trim(line);
      euler_levels.push_back(boost::lexical_cast<boost::uint32_t>(line));
    }
  }
  eulerLevelsFile.close();


  // parse the euler positions
  std::ifstream eulerPositionsFile;
  eulerPositionsFile.open(args.euler_positions.c_str());
  if (!eulerPositionsFile.is_open()) {
    std::cerr << "Could not open file: " << args.euler_positions << std::endl;
    return EXIT_FAILURE;
  }

#ifndef NDEBUG
  std::cout << "Reading the Euler Positions..." << std::endl;
#endif /* NDEBUG */

  common::StringIntMap euler_positions;
  while (!eulerPositionsFile.eof()) {
    std::getline(eulerPositionsFile, line);
#ifndef NDEBUG
    std::cout << "Read line: " << line << std::endl;
#endif /* NDEBUG */

    if (line != "") {
      common::StringVec tokens;
      boost::split(tokens, line, boost::is_any_of(","), boost::token_compress_on);
      if (tokens.size() != 2) {
        eulerPositionsFile.close();
        std::cerr << "Each line of the euler positions file can only contain two tokens: <key>,<value>!" << std::endl;
        return EXIT_FAILURE;
      } else {
        if (euler_positions.find(tokens[0]) != euler_positions.end()) {
          eulerPositionsFile.close();
          std::cerr << "The Euler Positions cannot contain duplicates: " << tokens[0] << std::endl;
          return EXIT_FAILURE;
        } else {
          euler_positions[tokens[0]] = boost::lexical_cast<boost::uint32_t>(tokens[1]);
        }
      }
    }
  }
  eulerPositionsFile.close();

#ifndef NDEBUG
  std::cout << std::endl << "Euler Levels:  ";
  std::copy(euler_levels.begin(), euler_levels.end(), std::ostream_iterator<boost::uint32_t>(std::cout, "\t"));
  std::cout << std::endl;
  std::cout << std::endl << "Euler Positions:  ";
  std::copy(euler_positions.begin(), euler_positions.end(), std::ostream_iterator<common::StringIntMap::value_type>(std::cout, "\t"));
  std::cout << std::endl;
#endif /* NDEBUG */

  std::ifstream lcaFile;
  lcaFile.open(args.lca.c_str());
  if (!lcaFile.is_open()) {
    std::cerr << "Could not open file: " << args.lca << std::endl;
    return EXIT_FAILURE;
  }

#ifndef NDEBUG
  std::cout << "Reading the lca..." << std::endl;
#endif /* NDEBUG */

  common::StrStrMap lcas;
  while (!lcaFile.eof()) {
    std::getline(lcaFile, line);
#ifndef NDEBUG
    std::cout << "Read line: " << line << std::endl;
#endif /* NDEBUG */

    if (line != "") {
      common::StringVec tokens;
      boost::split(tokens, line, boost::is_any_of(","), boost::token_compress_on);
      if (tokens.size() != 3) {
        std::cerr << "The LCA requires three tokens, got: " << tokens.size() << std::endl;
        lcaFile.close();
        return EXIT_FAILURE;
      }
      lcas[std::make_tuple(tokens[0], tokens[1])] = common::Symbol(tokens[2]);
    }
  }
  lcaFile.close();

#ifndef NDEBUG
  std::cout << std::endl << "LCAs:  ";
  std::copy(lcas.begin(), lcas.end(), std::ostream_iterator<common::StrStrMap::value_type>(std::cout, "\t"));
  std::cout << std::endl;
#endif /* NDEBUG */

  std::ifstream set1File;
  set1File.open(args.set_1.c_str());
  if (!set1File.is_open()) {
    std::cerr << "Could not open file: " << args.set_1 << std::endl;
    return EXIT_FAILURE;
  }

#ifndef NDEBUG
  std::cout << "Reading the set1..." << std::endl;
#endif /* NDEBUG */

  common::Sequences seqs_1;
  while (!set1File.eof()) {
    std::getline(set1File, line);
#ifndef NDEBUG
    std::cout << "Read line: " << line << std::endl;
#endif /* NDEBUG */

    if (line != "") {
      common::StringVec cats;
      boost::split(cats, line, boost::is_any_of(","), boost::token_compress_on);
      seqs_1.push_back(cats);
    }
  }
  set1File.close();

#ifndef NDEBUG
  std::cout << std::endl << "1. Sequences:  ";
  std::copy(seqs_1.begin(), seqs_1.end(), std::ostream_iterator<common::StringVec>(std::cout, "\n"));
  std::cout << std::endl;
#endif /* NDEBUG */

  std::ifstream set2File;
  set2File.open(args.set_2.c_str());
  if (!set2File.is_open()) {
    std::cerr << "Could not open file: " << args.set_2 << std::endl;
    return EXIT_FAILURE;
  }

#ifndef NDEBUG
  std::cout << "Reading the set2..." << std::endl;
#endif /* NDEBUG */

  common::Sequences seqs_2;
  while (!set2File.eof()) {
    std::getline(set2File, line);
#ifndef NDEBUG
    std::cout << "Read line: " << line << std::endl;
#endif /* NDEBUG */

    if (line != "") {
      common::StringVec cats;
      boost::split(cats, line, boost::is_any_of(","), boost::token_compress_on);
      seqs_2.push_back(cats);
    }
  }
  set2File.close();

#ifndef NDEBUG
  std::cout << std::endl << "2. Sequences:  ";
  std::copy(seqs_2.begin(), seqs_2.end(), std::ostream_iterator<common::StringVec>(std::cout, "\n"));
  std::cout << std::endl;
#endif /* NDEBUG */

  alignment::AbstractDistanceMeasure *scoringScheme = new TreePathSimilarityMeasure(args.gap_penalty, euler_levels, euler_positions, lcas);
  alignment::SimilarityAlgorithm *similarity;

  if (args.alg == 1) {
    similarity = new alignment::SW(args.scores);
  } else {
    similarity = new alignment::NW(args.scores);
  }

  std::string outFile = args.results_dir + "/similarity-scores.dat";
  std::ofstream out(outFile.c_str(), std::ios::out);

  for (boost::uint32_t i = 0; i < seqs_1.size(); ++i) {
    #pragma omp parallel shared(std::cout, out, seqs_1, seqs_2, i, scoringScheme, similarity) default(none)
    {
      alignment::MemoryPool mem;
      std::vector<double> scores;

      #pragma omp for
      for (boost::uint32_t j = 0; j < seqs_2.size(); ++j) {
#ifndef NDEBUG
        std::cout << "i: ";
        std::copy(seqs_1[i].begin(), seqs_1[i].end(), std::ostream_iterator<std::string>(std::cout, " "));
        std::cout << std::endl;
        std::cout << "j: ";
        std::copy(seqs_2[j].begin(), seqs_2[j].end(), std::ostream_iterator<std::string>(std::cout, " "));
        std::cout << std::endl;
#endif /* NDEBUG */

        alignment::alignmentResult res = similarity->align(seqs_1[i], seqs_2[j], *scoringScheme, mem);
        // compute a normalised similarity measure
        double numerator = res.score * res.score;
        auto minSize = std::min(seqs_1[i].size(), seqs_2[j].size());
        double denominator = boost::lexical_cast<double>(minSize * minSize);
        double score = numerator/denominator;
        scores.push_back(score);

#ifndef NDEBUG
        std::cout << "Score: " << res.score << std::endl;
        std::cout << "Alignments: " << std::endl;
        std::copy(res.alignment[0].begin(), res.alignment[0].end(), std::ostream_iterator<std::string>(std::cout, " "));
        std::cout << std::endl;
#endif /* NDEBUG */
      }

      #pragma omp critical
      {
        for (std::vector<double>::iterator s_it = scores.begin(); s_it != scores.end(); ++s_it) {
          out << *s_it << std::endl;
        }
      }
    }
  }

  delete similarity;
  delete scoringScheme;

  return EXIT_SUCCESS;
}
