//
// File: MafFilter.cpp
// Created by: Julien Dutheil
// Created on: Jul 21 2010
//

/*
Copyright or © or Copr. Julien Y. Dutheil, (2010)

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

// From the STL:
#include <iostream>
#include <iomanip>
#include <string>
#include <memory>
using namespace std;

#include "OutputAsFeaturesMafIterator.h"
#include "SystemCallMafIterator.h"
#include "TreeBuildingSystemCallMafIterator.h"
#include "PopGenMafStatistics.h"

//From boost:
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>
using namespace boost::iostreams;

// From bpp-core:
#include <Bpp/App/BppApplication.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Text/StringTokenizer.h>

// From bpp-seq:
#include <Bpp/Seq/SequenceWithQuality.h>
#include <Bpp/Seq/Io/BppOSequenceStreamReaderFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentReaderFormat.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

// From bpp-seq-omics and bpp-phyl-omics:
#include <Bpp/Seq/Io/Maf/MafParser.h>
#include <Bpp/Seq/Io/Maf/SequenceStreamToMafIterator.h>
#include <Bpp/Seq/Io/Maf/SequenceFilterMafIterator.h>
#include <Bpp/Seq/Io/Maf/OrphanSequenceFilterMafIterator.h>
#include <Bpp/Seq/Io/Maf/BlockMergerMafIterator.h>
#include <Bpp/Seq/Io/Maf/BlockLengthMafIterator.h>
#include <Bpp/Seq/Io/Maf/BlockSizeMafIterator.h>
#include <Bpp/Seq/Io/Maf/ChromosomeMafIterator.h>
#include <Bpp/Seq/Io/Maf/ConcatenateMafIterator.h>
#include <Bpp/Seq/Io/Maf/RemoveEmptySequencesMafIterator.h>
#include <Bpp/Seq/Io/Maf/DuplicateFilterMafIterator.h>
#include <Bpp/Seq/Io/Maf/OrderFilterMafIterator.h>
#include <Bpp/Seq/Io/Maf/FeatureFilterMafIterator.h>
#include <Bpp/Seq/Io/Maf/QualityFilterMafIterator.h>
#include <Bpp/Seq/Io/Maf/MaskFilterMafIterator.h>
#include <Bpp/Seq/Io/Maf/EntropyFilterMafIterator.h>
#include <Bpp/Seq/Io/Maf/OutputMafIterator.h>
#include <Bpp/Seq/Io/Maf/AlignmentFilterMafIterator.h>
#include <Bpp/Seq/Io/Maf/PlinkOutputMafIterator.h>
#include <Bpp/Seq/Io/Maf/SequenceLDhotOutputMafIterator.h>
#include <Bpp/Seq/Io/Maf/CoordinateTranslatorMafIterator.h>
#include <Bpp/Seq/Io/Maf/CoordinatesOutputMafIterator.h>
#include <Bpp/Seq/Io/Maf/FullGapFilterMafIterator.h>
#include <Bpp/Seq/Io/Maf/FilterTreeMafIterator.h>
#include <Bpp/Seq/Io/Maf/OutputTreeMafIterator.h>
#include <Bpp/Seq/Io/Maf/OutputAlignmentMafIterator.h>
#include <Bpp/Seq/Io/Maf/OutputDistanceMatrixMafIterator.h>
#include <Bpp/Seq/Io/Maf/VcfOutputMafIterator.h>
#include <Bpp/Seq/Io/Maf/MsmcOutputMafIterator.h>
#include <Bpp/Seq/Io/Maf/TableOutputMafIterator.h>
#include <Bpp/Seq/Io/Maf/SequenceStatisticsMafIterator.h>
#include <Bpp/Seq/Io/Maf/FeatureExtractorMafIterator.h>
#include <Bpp/Seq/Io/Maf/WindowSplitMafIterator.h>
#include <Bpp/Seq/Io/Maf/CountDistanceEstimationMafIterator.h>
#include <Bpp/Seq/Io/Maf/MaximumLikelihoodDistanceEstimationMafIterator.h>
#include <Bpp/Seq/Io/Maf/DistanceBasedPhylogenyReconstructionMafIterator.h>
#include <Bpp/Seq/Io/Maf/TreeManipulationMafIterators.h>
#include <Bpp/Seq/Io/Maf/MafStatistics.h>
#include <Bpp/Seq/Io/Maf/CountClustersMafStatistics.h>
#include <Bpp/Seq/Io/Maf/MaximumLikelihoodModelFitMafStatistics.h>
#include <Bpp/Seq/Io/Maf/IterationListener.h>
#include <Bpp/Seq/Feature/Gff/GffFeatureReader.h>
#include <Bpp/Seq/Feature/Gtf/GtfFeatureReader.h>
#include <Bpp/Seq/Feature/Bed/BedGraphFeatureReader.h>

// From bpp-phyl:
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/NeighborJoining.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Io/BppOSubstitutionModelFormat.h>
#include <Bpp/Phyl/Io/BppORateDistributionFormat.h>
#include <Bpp/Phyl/Io/BppOTreeReaderFormat.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "maffilter name1=value1 name2=value2").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the MafFilter Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "  Online version: https://jydu.github.io/maffilter/Manual/index.html").endLine();
  (*ApplicationTools::message << "  Or type 'info maffilter' in a terminal.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                  MAF Filter, version 1.3.1                     *" << endl;
  cout << "* Author: J. Dutheil                        Created on  10/09/10 *" << endl;
  cout << "*                                           Last Modif. 18/08/18 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    exit(0);
  }
  
  try
  {
    BppApplication maffilter(args, argv, "MafFilter");
    maffilter.startTimer();

    string inputFile = ApplicationTools::getAFilePath("input.file", maffilter.getParams(), true, true);
    string inputFormat = ApplicationTools::getStringParameter("input.format", maffilter.getParams(), "Maf", "", true, false);
    string compress = ApplicationTools::getStringParameter("input.file.compression", maffilter.getParams(), "none");
    string inputDot = ApplicationTools::getStringParameter("input.dots", maffilter.getParams(), "error", "", true, false);

    filtering_istream stream;
    if (compress == "none") {
    } else if (compress == "gzip") {
      stream.push(gzip_decompressor());
    } else if (compress == "zip") {
      stream.push(zlib_decompressor());
    } else if (compress == "bzip2") {
      stream.push(bzip2_decompressor());
    } else
      throw Exception("Bad input incompression format: " + compress);
    stream.push(file_source(inputFile));
    
    string logFile = ApplicationTools::getAFilePath("output.log", maffilter.getParams(), false, false);
    shared_ptr<StlOutputStream> log;
    if (logFile != "none") {
      log.reset(new StlOutputStream(new ofstream(logFile.c_str(), ios::out)));
    }

    MafIterator* currentIterator;

    if (inputFormat == "Maf") {
      short dotOption = MafParser::DOT_ERROR;
      if (inputDot == "as_gaps") {
        ApplicationTools::displayResult("Maf 'dotted' alignment input", string("converted to gaps"));
        dotOption = MafParser::DOT_ASGAP;
      } else if (inputDot == "as_unresolved") {
        ApplicationTools::displayResult("Maf 'dotted' alignment input", string("converted to unresolved"));
        dotOption = MafParser::DOT_ASUNRES;
      }

      bool checkSize = ApplicationTools::getBooleanParameter("input.check_sequence_size", maffilter.getParams(), true, "", true, false);
      if (!checkSize)
        ApplicationTools::displayBooleanResult("Check size of sequences", false);
      currentIterator = new MafParser(&stream, true, checkSize, dotOption);
    } else {
      if (inputDot == "as_gaps") throw Exception("'dot_as_gaps' option only available with Maf input.");
      BppOSequenceStreamReaderFormat reader;
      ISequenceStream* seqStream = reader.read(inputFormat);
      map<string, string> cmdArgs(reader.getUnparsedArguments());
      bool zeroBased = ApplicationTools::getBooleanParameter("zero_based", cmdArgs, true);
      currentIterator = new SequenceStreamToMafIterator(seqStream, &stream, false, zeroBased);
    }
    
    ApplicationTools::displayResult("Reading file", inputFile + " as " + inputFormat + (compress == "none" ? "" : "(" + compress + ")"));
    ApplicationTools::displayResult("Output log file", logFile);


    vector<string> actions = ApplicationTools::getVectorParameter<string>("maf.filter", maffilter.getParams(), ',', "", "", false, false);
    vector<MafIterator*> its;
    its.push_back(currentIterator);
    vector<filtering_ostream*> ostreams;
    for (size_t a = 0; a < actions.size(); a++) {
      string cmdName;
      map<string, string> cmdArgs;
      KeyvalTools::parseProcedure(actions[a], cmdName, cmdArgs);
      (*ApplicationTools::message << "-------------------------------------------------------------------").endLine();
      ApplicationTools::displayResult("Adding filter", cmdName);
      
      bool verbose = ApplicationTools::getBooleanParameter("verbose", cmdArgs, true, "", true, false);
      ApplicationTools::displayBooleanResult("-- Verbose", verbose);

      // +-----------------+
      // | Sequence subset |
      // +-----------------+
      if (cmdName == "Subset") {
        bool strict = ApplicationTools::getBooleanParameter("strict", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- All species should be in output blocks", strict);
        bool keep = ApplicationTools::getBooleanParameter("keep", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Sequences not in the list will be kept", keep);
        if (cmdArgs.find("rm.duplicates") != cmdArgs.end()) {
          throw Exception("rm.duplicates argument in Subset is deprecated: use remove_duplicates instead.");
        }
        bool rmdupl = ApplicationTools::getBooleanParameter("remove_duplicates", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Species should be present only once", rmdupl);
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'Subset'.");
        SequenceFilterMafIterator* iterator = new SequenceFilterMafIterator(currentIterator, species, strict, keep, rmdupl);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------------------+
      // | Sequence orphan selection |
      // +---------------------------+
      else if (cmdName == "SelectOrphans") {
        bool strict = ApplicationTools::getBooleanParameter("strict", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- All species should be in output blocks", strict);
        bool rmdupl = ApplicationTools::getBooleanParameter("remove_duplicates", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Species should be present only once", rmdupl);
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'SelectOrphans'.");
        OrphanSequenceFilterMafIterator* iterator = new OrphanSequenceFilterMafIterator(currentIterator, species, strict, rmdupl);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------+
      // | Block merging |
      // +---------------+
      else if (cmdName == "Merge") {
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'Merge'.");

        if (cmdArgs.find("dist.max") != cmdArgs.end()) {
          throw Exception("dist.max argument in Merge is deprecated: use dist_max instead.");
        }
        unsigned int distMax = ApplicationTools::getParameter<unsigned int>("dist_max", cmdArgs, 0);
        ApplicationTools::displayResult("-- Maximum distance allowed", distMax);
        bool renameChimeras = ApplicationTools::getBooleanParameter("rename_chimeric_chromosomes", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Rename chimeric chromosomes", renameChimeras);
        BlockMergerMafIterator* iterator = new BlockMergerMafIterator(currentIterator, species, distMax, renameChimeras);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        if (cmdArgs.find("ignore.chr") != cmdArgs.end()) {
          throw Exception("ignore.chr argument in Merge is deprecated: use ignore_chr instead.");
        }
        string ignoreChrList = ApplicationTools::getStringParameter("ignore_chr", cmdArgs, "none");
        if (ignoreChrList != "none") {
          if (ignoreChrList[0] == '(') {
            StringTokenizer st(ignoreChrList.substr(1, ignoreChrList.size() - 2), ",");
            while (st.hasMoreToken())
              iterator->ignoreChromosome(st.nextToken());
          } else {
            iterator->ignoreChromosome(ignoreChrList);
          }
        }
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------------+
      // | Block concatenation |
      // +---------------------+
      else if (cmdName == "Concatenate") {
        unsigned int minimumSize = ApplicationTools::getParameter<unsigned int>("minimum_size", cmdArgs, 0);
        string ref = ApplicationTools::getStringParameter("ref_species", cmdArgs, "", "", true, 2);
        ApplicationTools::displayResult("-- Minimum final block size", minimumSize);
        if (ref != "")
          ApplicationTools::displayResult("-- Reference species", ref);
        ConcatenateMafIterator* iterator = new ConcatenateMafIterator(currentIterator, minimumSize, ref);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +--------------------+
      // | Full gap filtering |
      // +--------------------+
      else if (cmdName == "XFullGap") {
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'XFullGap'.");
        FullGapFilterMafIterator* iterator = new FullGapFilterMafIterator(currentIterator, species);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------------+
      // | Alignment filtering |
      // +---------------------+
      else if (cmdName == "AlnFilter") {
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'AlnFilter'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        bool relative   = ApplicationTools::getBooleanParameter("relative", cmdArgs, false);
        unsigned int gm = 0;
        double rm = 0;
        if (relative)
          rm = ApplicationTools::getDoubleParameter("max.gap", cmdArgs, 0);
        else
          gm = ApplicationTools::getParameter<unsigned int>("max.gap", cmdArgs, 0);
        double em       = ApplicationTools::getParameter<double>("max.ent", cmdArgs, 0); //Default means no entropy threshold
        bool missingAsGap = ApplicationTools::getParameter<bool>("missing_as_gap", cmdArgs, false);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Window size", ws);
        ApplicationTools::displayResult("-- Window step", st);
        if (relative)
          ApplicationTools::displayResult("-- Max. gaps allowed in Window", TextTools::toString(rm * 100) + "%");
        else 
          ApplicationTools::displayResult("-- Max. gaps allowed in Window", gm);
        ApplicationTools::displayResult("-- Max. total entropy in Window", em);
        ApplicationTools::displayBooleanResult("-- Missing sequence replaced by gaps", missingAsGap);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
        AlignmentFilterMafIterator* iterator;
        if (relative)
          iterator = new AlignmentFilterMafIterator(currentIterator, species, ws, st, rm, em, !trash, missingAsGap);
        else
          iterator = new AlignmentFilterMafIterator(currentIterator, species, ws, st, gm, em, !trash, missingAsGap);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        its.push_back(iterator);

        if (!trash) {
          compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
          filtering_ostream* out = new filtering_ostream;
          if (compress == "none") {
          } else if (compress == "gzip") {
            out->push(gzip_compressor());
          } else if (compress == "zip") {
            out->push(zlib_compressor());
          } else if (compress == "bzip2") {
            out->push(bzip2_compressor());
          } else
            throw Exception("Bad output compression format: " + compress);
          out->push(file_sink(outputFile));
          ostreams.push_back(out);
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }


      // +-----------------------+
      // | Alignment filtering 2 |
      // +-----------------------+
      else if (cmdName == "AlnFilter2") {
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'AlnFilter2'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        bool relative   = ApplicationTools::getBooleanParameter("relative", cmdArgs, false);
        unsigned int gm = 0;
        double rm = 0;
        if (relative)
          rm = ApplicationTools::getDoubleParameter("max.gap", cmdArgs, 0);
        else
          gm = ApplicationTools::getParameter<unsigned int>("max.gap", cmdArgs, 0);
        unsigned int pm = ApplicationTools::getParameter<unsigned int>("max.pos", cmdArgs, 0);
        bool missingAsGap = ApplicationTools::getParameter<bool>("missing_as_gap", cmdArgs, false);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Window size", ws);
        ApplicationTools::displayResult("-- Window step", st);
        if (relative)
          ApplicationTools::displayResult("-- Max. gaps allowed per position", TextTools::toString(rm * 100) + "%");
        else 
          ApplicationTools::displayResult("-- Max. gaps allowed per position", gm);
        ApplicationTools::displayResult("-- Max. gap positions allowed", pm);
        ApplicationTools::displayBooleanResult("-- Missing sequence replaced by gaps", missingAsGap);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
        AlignmentFilter2MafIterator* iterator;
        if (relative)
          iterator = new AlignmentFilter2MafIterator(currentIterator, species, ws, st, rm, pm, !trash, missingAsGap);
        else
          iterator = new AlignmentFilter2MafIterator(currentIterator, species, ws, st, gm, pm, !trash, missingAsGap);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        its.push_back(iterator);

        if (!trash) {
          compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
          filtering_ostream* out = new filtering_ostream;
          if (compress == "none") {
          } else if (compress == "gzip") {
            out->push(gzip_compressor());
          } else if (compress == "zip") {
            out->push(zlib_compressor());
          } else if (compress == "bzip2") {
            out->push(bzip2_compressor());
          } else
            throw Exception("Bad output compression format: " + compress);
          out->push(file_sink(outputFile));
          ostreams.push_back(out);
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }



      // +-------------------+
      // | Entropy filtering |
      // +-------------------+
      else if (cmdName == "EntropyFilter") {
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'AlnFilter2'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        double       em = ApplicationTools::getParameter<double>      ("max.ent", cmdArgs, 0);
        unsigned int pm = ApplicationTools::getParameter<unsigned int>("max.pos", cmdArgs, 0);
        bool missingAsGap = ApplicationTools::getParameter<bool>("missing_as_gap", cmdArgs, false);
        bool ignoreGaps   = ApplicationTools::getParameter<bool>("ignore_gaps", cmdArgs, false);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Window size", ws);
        ApplicationTools::displayResult("-- Window step", st);
        ApplicationTools::displayResult("-- Max. entropy allowed per position", em);
        ApplicationTools::displayResult("-- Max. high entropy positions allowed", pm);
        ApplicationTools::displayBooleanResult("-- Missing sequence replaced by gaps", missingAsGap);
        ApplicationTools::displayBooleanResult("-- Gaps should be ignored", ignoreGaps);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
        if (ignoreGaps && missingAsGap)
          throw Exception("Error, incompatible options ingore_gaps=yes and missing_as_gap=yes.");
        EntropyFilterMafIterator* iterator = new EntropyFilterMafIterator(currentIterator, species, ws, st, em, pm, !trash, missingAsGap, ignoreGaps);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        its.push_back(iterator);

        if (!trash) {
          compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
          filtering_ostream* out = new filtering_ostream;
          if (compress == "none") {
          } else if (compress == "gzip") {
            out->push(gzip_compressor());
          } else if (compress == "zip") {
            out->push(zlib_compressor());
          } else if (compress == "bzip2") {
            out->push(bzip2_compressor());
          } else
            throw Exception("Bad output compression format: " + compress);
          out->push(file_sink(outputFile));
          ostreams.push_back(out);
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }



      // +----------------+
      // | Mask filtering |
      // +----------------+
      else if (cmdName == "MaskFilter") {
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'MaskFilter'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        unsigned int mm = ApplicationTools::getParameter<unsigned int>("max.masked", cmdArgs, 0);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Window size", ws);
        ApplicationTools::displayResult("-- Window step", st);
        ApplicationTools::displayResult("-- Max. masked sites allowed in Window", mm);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
        MaskFilterMafIterator* iterator = new MaskFilterMafIterator(currentIterator, species, ws, st, mm, !trash);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        its.push_back(iterator);

        if (!trash) {
          compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
          filtering_ostream* out = new filtering_ostream;
          if (compress == "none") {
          } else if (compress == "gzip") {
            out->push(gzip_compressor());
          } else if (compress == "zip") {
            out->push(zlib_compressor());
          } else if (compress == "bzip2") {
            out->push(bzip2_compressor());
          } else
            throw Exception("Bad output compression format: " + compress);
          out->push(file_sink(outputFile));
          ostreams.push_back(out);
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }


      // +-------------------+
      // | Quality filtering |
      // +-------------------+
      else if (cmdName == "QualFilter") {
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for command 'QualFilter'.");
        unsigned int ws = ApplicationTools::getParameter<unsigned int>("window.size", cmdArgs, 10);
        unsigned int st = ApplicationTools::getParameter<unsigned int>("window.step", cmdArgs, 5);
        double       mq = ApplicationTools::getDoubleParameter("min.qual", cmdArgs, 0);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Window size", ws);
        ApplicationTools::displayResult("-- Window step", st);
        ApplicationTools::displayResult("-- Min. average quality allowed in Window", mq);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
        QualityFilterMafIterator* iterator = new QualityFilterMafIterator(currentIterator, species, ws, st, mq, !trash);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        its.push_back(iterator);

        if (!trash) {
          compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
          filtering_ostream* out = new filtering_ostream;
          if (compress == "none") {
          } else if (compress == "gzip") {
            out->push(gzip_compressor());
          } else if (compress == "zip") {
            out->push(zlib_compressor());
          } else if (compress == "bzip2") {
            out->push(bzip2_compressor());
          } else
            throw Exception("Bad output compression format: " + compress);
          out->push(file_sink(outputFile));
          ostreams.push_back(out);
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }


      
      // +-------------------------+
      // | Feature-based filtering |
      // +-------------------------+
      else if (cmdName == "FeatureFilter") {
        string refSpecies = ApplicationTools::getStringParameter("ref_species", cmdArgs, "none");
        string featureFile = ApplicationTools::getAFilePath("feature.file", cmdArgs, false, false);
        string featureFormat = ApplicationTools::getStringParameter("feature.format", cmdArgs, "GFF");
        vector<string> featureType = ApplicationTools::getVectorParameter<string>("feature.type", cmdArgs, ',', "all");
        if (featureType.size() == 0)
          throw Exception("At least one feature should be provided for command 'FeatureFilter'.");
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        ApplicationTools::displayResult("-- Features to remove", featureFile + " (" + featureFormat + ")");
        ApplicationTools::displayResult("-- Features are for species", refSpecies);
        ApplicationTools::displayBooleanResult("-- Output removed blocks", !trash);
        compress = ApplicationTools::getStringParameter("feature.file.compression", cmdArgs, "none");
        filtering_istream featureStream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          featureStream.push(gzip_decompressor());
        } else if (compress == "zip") {
          featureStream.push(zlib_decompressor());
        } else if (compress == "bzip2") {
          featureStream.push(bzip2_decompressor());
        } else
          throw Exception("Bad input incompression format: " + compress);
        featureStream.push(file_source(featureFile));
        unique_ptr<FeatureReader> ftReader;
        SequenceFeatureSet featuresSet;
        if (featureFormat == "GFF") {
          ftReader.reset(new GffFeatureReader(featureStream));
        } else if (featureFormat == "GTF") {
          ftReader.reset(new GtfFeatureReader(featureStream));
        } else if (featureFormat == "BedGraph") {
          ftReader.reset(new BedGraphFeatureReader(featureStream));
         } else
          throw Exception("Unsupported feature format: " + featureFormat);
        if (featureType.size() == 1 && featureType[0] == "all")
          ftReader->getAllFeatures(featuresSet);
        else {
          for (size_t i = 0; i < featureType.size(); ++i) {
            ApplicationTools::displayResult("-- Filter features of type", featureType[i]);
            ftReader->getFeaturesOfType(featureType[i], featuresSet);
          }
        }
        ApplicationTools::displayResult("-- Total number of features", featuresSet.getNumberOfFeatures());
        FeatureFilterMafIterator* iterator = new FeatureFilterMafIterator(currentIterator, refSpecies, featuresSet, !trash);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        its.push_back(iterator);

        if (!trash) {
          compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
          filtering_ostream* out = new filtering_ostream;
          if (compress == "none") {
          } else if (compress == "gzip") {
            out->push(gzip_compressor());
          } else if (compress == "zip") {
            out->push(zlib_compressor());
          } else if (compress == "bzip2") {
            out->push(bzip2_compressor());
          } else
            throw Exception("Bad output compression format: " + compress);
          out->push(file_sink(outputFile));
          ostreams.push_back(out);
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }



      // +------------------------+
      // | Block length filtering |
      // +------------------------+
      else if (cmdName == "MinBlockLength") {
        if (cmdArgs.find("min.length") != cmdArgs.end()) {
          throw Exception("min.length argument in MinBlockLength is deprecated: use min_length instead.");
        }
        unsigned int minLength = ApplicationTools::getParameter<unsigned int>("min_length", cmdArgs, 0);
        ApplicationTools::displayResult("-- Minimum block length required", minLength);
        BlockLengthMafIterator* iterator = new BlockLengthMafIterator(currentIterator, minLength);
        iterator->setLogStream(log);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +----------------------+
      // | Block size filtering |
      // +----------------------+
      else if (cmdName == "MinBlockSize") {
        if (cmdArgs.find("min.size") != cmdArgs.end()) {
          throw Exception("min.size argument in MinBlockSize is deprecated: use min_size instead.");
        }
        unsigned int minSize = ApplicationTools::getParameter<unsigned int>("min_size", cmdArgs, 0);
        ApplicationTools::displayResult("-- Minimum block size required", minSize);
        if (minSize > 5)
          ApplicationTools::displayWarning("!! Warning, in previous version of maffilter BlockLength was named BlockSize... Check!");
        BlockSizeMafIterator* iterator = new BlockSizeMafIterator(currentIterator, minSize);
        iterator->setLogStream(log);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +----------------------+
      // | Chromosome filtering |
      // +----------------------+
      else if (cmdName == "SelectChr") {
        if (cmdArgs.find("reference") != cmdArgs.end()) {
          throw Exception("reference argument in SelectChr is deprecated: use ref_species instead.");
        }
        string ref = ApplicationTools::getStringParameter("ref_species", cmdArgs, "");
        ApplicationTools::displayResult("-- Reference species", ref);
        string chr = ApplicationTools::getStringParameter("chromosome", cmdArgs, "");
        ApplicationTools::displayResult("-- Chromosome", chr);
        ChromosomeMafIterator* iterator = new ChromosomeMafIterator(currentIterator, ref, chr);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------------+
      // | Duplicate filtering |
      // +---------------------+
      //Nb: this is kind of deprecated, should be done better by looking at partial overlap.
      //could be useful for debugging though. We do not report it in the documentation for now.
      else if (cmdName == "DuplicateFilter") {
        string ref = ApplicationTools::getStringParameter("reference", cmdArgs, "");
        ApplicationTools::displayResult("-- Reference species", ref);
        DuplicateFilterMafIterator* iterator = new DuplicateFilterMafIterator(currentIterator, ref);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +-----------------+
      // | Order filtering |
      // +-----------------+
      else if (cmdName == "OrderFilter") {
        string ref = ApplicationTools::getStringParameter("reference", cmdArgs, "");
        ApplicationTools::displayResult("-- Reference species", ref);

        string unsortedAction = ApplicationTools::getStringParameter("do_unsorted", cmdArgs, "exception");
        ApplicationTools::displayResult("-- Unsorted block action", unsortedAction);
        if (unsortedAction != "none" && unsortedAction != "discard" && unsortedAction != "error")
          throw Exception("'do_unsorted' must be one of 'none', 'discard' or 'error'.");

        string overlappingAction = ApplicationTools::getStringParameter("do_overlapping", cmdArgs, "exception");
        ApplicationTools::displayResult("-- Overlapping block action", overlappingAction);
        if (overlappingAction != "none" && overlappingAction != "discard" && overlappingAction != "error")
          throw Exception("'do_overlapping' must be one of 'none', 'discard' or 'error'.");

        OrderFilterMafIterator* iterator = new OrderFilterMafIterator(currentIterator, ref,
            unsortedAction == "discard", unsortedAction == "error",
            overlappingAction == "discard", overlappingAction == "error");
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      // +---------------------------+
      // | Empty sequences filtering |
      // +---------------------------+
      else if (cmdName == "RemoveEmptySequences") {
        bool unresolvedAsGaps = ApplicationTools::getBooleanParameter("unresolved_as_gaps", cmdArgs, "");
        ApplicationTools::displayBooleanResult("-- Unresolved as gaps", unresolvedAsGaps);
        RemoveEmptySequencesMafIterator* iterator = new RemoveEmptySequencesMafIterator(currentIterator, unresolvedAsGaps);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }

      
      // +----------------+
      // | Tree filtering |
      // +----------------+
      else if (cmdName == "TreeFilter") {
        string treeProperty = ApplicationTools::getStringParameter("tree", cmdArgs, "none");
        double maxBrLen = ApplicationTools::getDoubleParameter("max_brlen", cmdArgs, 0.1);
        
        ApplicationTools::displayResult("-- Max. branch length", maxBrLen);
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, false, false);
        bool trash = outputFile == "none";
        FilterTreeMafIterator* iterator = new FilterTreeMafIterator(currentIterator, treeProperty, maxBrLen, !trash);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        its.push_back(iterator);

        if (!trash) {
          compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
          filtering_ostream* out = new filtering_ostream;
          if (compress == "none") {
          } else if (compress == "gzip") {
            out->push(gzip_compressor());
          } else if (compress == "zip") {
            out->push(zlib_compressor());
          } else if (compress == "bzip2") {
            out->push(bzip2_compressor());
          } else
            throw Exception("Bad output compression format: " + compress);
          out->push(file_sink(outputFile));
          ostreams.push_back(out);
          ApplicationTools::displayResult("-- File compression for removed blocks", compress);

          //Now build an adaptor for retrieving the trashed blocks:
          TrashIteratorAdapter* trashIt = new TrashIteratorAdapter(iterator);
          //Add an output iterator:
          OutputMafIterator* outIt = new OutputMafIterator(trashIt, out);
          //And then synchronize the two iterators:
          MafIteratorSynchronizer* syncIt = new MafIteratorSynchronizer(iterator, outIt);
          //Returns last iterator:
          currentIterator = syncIt;
          //Keep track of all those iterators:
          its.push_back(trashIt);
          its.push_back(syncIt);
        } else {
          //We only get the remaining blocks here:
          currentIterator = iterator;
        }
      }





      // +---------------------+
      // | Sequence statistics |
      // +---------------------+
      else if (cmdName == "SequenceStatistics") {
        vector<string> statisticsDesc = ApplicationTools::getVectorParameter<string>("statistics", cmdArgs, ',', "", "", false, true);
        
        //Parse all statistics:
        vector<MafStatistics*> statistics;
        for (size_t i = 0; i < statisticsDesc.size(); ++i) {
          string statName;
          map<string, string> statArgs;
          KeyvalTools::parseProcedure(statisticsDesc[i], statName, statArgs);
          MafStatistics* mafStat = 0;
          string statDesc = "";
          if (statName == "BlockSize") {
            mafStat = new BlockSizeMafStatistics();
          } else if (statName == "BlockLength") {
            mafStat = new BlockLengthMafStatistics();
          } else if (statName == "SequenceLength") {
            string sp = ApplicationTools::getStringParameter("species", statArgs, "");
            mafStat = new SequenceLengthMafStatistics(sp);
          } else if (statName == "AlnScore") {
            mafStat = new AlignmentScoreMafStatistics();
          } else if (statName == "BlockCounts") {
            vector<string> species = ApplicationTools::getVectorParameter<string>("species", statArgs, ',', "", "", false, true);
            string suffix = ApplicationTools::getStringParameter("suffix", statArgs, "");
            mafStat = new CharacterCountsMafStatistics(&AlphabetTools::DNA_ALPHABET, species, suffix);
          } else if (statName == "PairwiseDivergence") {
            string sp1 = ApplicationTools::getStringParameter("species1", statArgs, "");
            string sp2 = ApplicationTools::getStringParameter("species2", statArgs, "");
            mafStat = new PairwiseDivergenceMafStatistics(sp1, sp2);
          } else if (statName == "SiteFrequencySpectrum") {
            vector<double> bounds  = ApplicationTools::getVectorParameter<double>("bounds", statArgs, ',', "", "", false, true);
            vector<string> ingroup = ApplicationTools::getVectorParameter<string>("ingroup", statArgs, ',', "", "", false, true);
            if (ingroup.size() < 2)
              throw Exception("ERROR: at least two ingroup sequences are required to compute the site frequency spectrum.");
            string outgroup        = ApplicationTools::getStringParameter("outgroup", statArgs, "", "", false, true);
            mafStat = new SiteFrequencySpectrumMafStatistics(&AlphabetTools::DNA_ALPHABET, bounds, ingroup, outgroup); 
          } else if (statName == "FourSpeciesSitePatternCounts") {
            string species1 = ApplicationTools::getStringParameter("species1", statArgs, "sp1", "", false, true);
            string species2 = ApplicationTools::getStringParameter("species2", statArgs, "sp2", "", false, true);
            string species3 = ApplicationTools::getStringParameter("species3", statArgs, "sp3", "", false, true);
            string species4 = ApplicationTools::getStringParameter("species4", statArgs, "sp4", "", false, true);
            vector<string> species;
            species.push_back(species1);
            species.push_back(species2);
            species.push_back(species3);
            species.push_back(species4);
            mafStat = new FourSpeciesPatternCountsMafStatistics(&AlphabetTools::DNA_ALPHABET, species); 
          } else if (statName == "SiteStatistics") {
            vector<string> species = ApplicationTools::getVectorParameter<string>("species", statArgs, ',', "", "", false, true);
            mafStat = new SiteMafStatistics(species); 
          } else if (statName == "PolymorphismStatistics") {
            vector<string> species1 = ApplicationTools::getVectorParameter<string>("species1", statArgs, ',', "", "", false, true);
            vector<string> species2 = ApplicationTools::getVectorParameter<string>("species2", statArgs, ',', "", "", false, true);
            vector< vector<string> > species;
            species.push_back(species1);
            species.push_back(species2);
            mafStat = new PolymorphismMafStatistics(species); 
          } else if (statName == "DiversityStatistics") {
            vector<string> species = ApplicationTools::getVectorParameter<string>("ingroup", statArgs, ',', "", "", false, true);
            if (species.size() < 2)
              throw Exception("ERROR: at least two sequences are required to compute diversity estimators.");
            mafStat = new SequenceDiversityMafStatistics(species); 
          } else if (statName == "FstStatistics") {
            vector<string> species1 = ApplicationTools::getVectorParameter<string>("species1", statArgs, ',', "", "", false, true);
            vector<string> species2 = ApplicationTools::getVectorParameter<string>("species2", statArgs, ',', "", "", false, true);
            unsigned int minNbPermutations = ApplicationTools::getParameter<unsigned int>("min_permutation_number", statArgs, 0, "", false, true);
            unsigned int maxNbPermutations = ApplicationTools::getParameter<unsigned int>("max_permutation_number", statArgs, 0, "", false, true);
            bool verboseStat = ApplicationTools::getBooleanParameter("verbose", statArgs, true, "", false, true);
            if (minNbPermutations > 0 || maxNbPermutations > 0) {
              ApplicationTools::displayResult("-- Min. Nb. permutations", minNbPermutations);
              ApplicationTools::displayResult("-- Max. Nb. permutations", maxNbPermutations);
            }
            mafStat = new FstMafStatistics(species1, species2, minNbPermutations, maxNbPermutations, verboseStat); 
          } else if (statName == "CountClusters") {
            string treeProperty = ApplicationTools::getStringParameter("tree", statArgs, "none");
            double threshold = ApplicationTools::getDoubleParameter("threshold", statArgs, 0);
            mafStat = new CountClustersMafStatistics(treeProperty, threshold);
            statDesc = " / " + treeProperty;
          } else if (statName == "ModelFit") {
            unique_ptr<SubstitutionModel> model;
            unique_ptr<SubstitutionModelSet> modelSet;
            unique_ptr<FrequenciesSet> rootFreqs;

            string modelType = ApplicationTools::getStringParameter("model_type", statArgs, "Homogeneous");
            if (modelType == "Homogeneous") {
              model.reset(PhylogeneticsApplicationTools::getSubstitutionModel(&AlphabetTools::DNA_ALPHABET, 0, 0, statArgs, "", true, true));
              ApplicationTools::displayResult("-- Substitution model", model->getName());
              string freqDescription = ApplicationTools::getStringParameter("root_freq", statArgs, "None");
              if (freqDescription != "None") {
                rootFreqs.reset(PhylogeneticsApplicationTools::getFrequenciesSet(
                    &AlphabetTools::DNA_ALPHABET, 0, freqDescription, 0, vector<double>(), 0));
                ApplicationTools::displayResult("-- Root frequencies", rootFreqs->getName());
              }
            } else if (modelType == "Nonhomogeneous") {
              modelSet.reset(PhylogeneticsApplicationTools::getSubstitutionModelSet(&AlphabetTools::DNA_ALPHABET, 0, 0, statArgs, "", true, true));
            } else {
              throw Exception("Unknown model type: " + modelType + ". Must be either Homogeneous or Nonhomogeneous.");
            }
            
            unique_ptr<DiscreteDistribution> rDist(PhylogeneticsApplicationTools::getRateDistribution(statArgs, "", true, false));
            ApplicationTools::displayResult("-- Rate distribution", rDist->getName());
            string treeProperty = ApplicationTools::getStringParameter("tree", statArgs, "none");
            vector<string> parametersOutput = ApplicationTools::getVectorParameter<string>("parameters_output", statArgs, ',', "");
            vector<string> fixedParametersNames = ApplicationTools::getVectorParameter<string>("fixed_parameters", statArgs, ',', "");
            ParameterList fixedParameters;
            if (fixedParametersNames.size() > 0) {
              ParameterList parameters;
              if (modelSet.get())
                parameters = modelSet->getParameters();
              else
                parameters = model->getParameters();
              parameters.addParameters(rDist->getParameters());
              fixedParameters = parameters.subList(fixedParametersNames);
            }
            bool reestimateBrLen = ApplicationTools::getBooleanParameter("reestimate_brlen", statArgs, true);
            ApplicationTools::displayBooleanResult("-- Reestimate branch lengths", reestimateBrLen);
            double propGapsToKeep = ApplicationTools::getDoubleParameter("max_freq_gaps", statArgs, 0.);
            ApplicationTools::displayResult("-- Max. frequency of gaps", propGapsToKeep);
            bool gapsAsUnresolved = ApplicationTools::getBooleanParameter("gaps_as_unresolved", statArgs, true);
            ApplicationTools::displayBooleanResult("-- Gaps as unresolved", gapsAsUnresolved);
            bool useClock = ApplicationTools::getBooleanParameter("global_clock", statArgs, false);
            ApplicationTools::displayBooleanResult("-- Use a global molecular clock", useClock);
            bool reparametrize = ApplicationTools::getBooleanParameter("reparametrize", statArgs, false);
            ApplicationTools::displayBooleanResult("-- Reparametrization", reparametrize);
            if (treeProperty == "none") {
              unique_ptr<Tree> tree(PhylogeneticsApplicationTools::getTree(statArgs, "", "", true, false)); 
              if (modelSet.get()) {
                mafStat = new MaximumLikelihoodModelFitMafStatistics(modelSet.release(), rDist.release(), tree.release(), parametersOutput,
                    fixedParameters, reestimateBrLen, propGapsToKeep, gapsAsUnresolved, useClock, reparametrize);
              } else {
                mafStat = new MaximumLikelihoodModelFitMafStatistics(model.release(), rDist.release(), dynamic_cast<NucleotideFrequenciesSet*>(rootFreqs.release()), tree.release(), parametersOutput,
                    fixedParameters, reestimateBrLen, propGapsToKeep, gapsAsUnresolved, useClock, reparametrize);
              }
            } else {
              mafStat = new MaximumLikelihoodModelFitMafStatistics(model.release(), rDist.release(), dynamic_cast<NucleotideFrequenciesSet*>(rootFreqs.release()), treeProperty, parametersOutput,
                  fixedParameters, reestimateBrLen, propGapsToKeep, gapsAsUnresolved, useClock, reparametrize);
            }
          } else {
            throw Exception("Unknown statistic: " + statName);
          }
          statistics.push_back(mafStat);
          ApplicationTools::displayResult("-- Adding statistic", mafStat->getFullName() + " <" + mafStat->getShortName() + ">" + statDesc);
        }

        //Get output file:
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output file", outputFile);
        filtering_ostream* out = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          out->push(gzip_compressor());
        } else if (compress == "zip") {
          out->push(zlib_compressor());
        } else if (compress == "bzip2") {
          out->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        out->push(file_sink(outputFile));
        ostreams.push_back(out);
        //ostreams.push_back(out);
        ApplicationTools::displayResult("-- File compression", compress);
        //StlOutputStream* output = new StlOutputStream(new ofstream(outputFile.c_str(), ios::out));
        StlOutputStream* output = new StlOutputStream(out);

        SequenceStatisticsMafIterator* iterator = new SequenceStatisticsMafIterator(currentIterator, statistics);
        
        if (cmdArgs.find("reference") != cmdArgs.end()) {
          throw Exception("reference argument in SequenceStatistics is deprecated: use ref_species instead.");
        }
        string ref = ApplicationTools::getStringParameter("ref_species", cmdArgs, "none");
        ApplicationTools::displayResult("-- Reference species", ref);
        CsvStatisticsOutputIterationListener* listener = new CsvStatisticsOutputIterationListener(iterator, ref, output);
        
        iterator->addIterationListener(listener);
        currentIterator = iterator;
        iterator->setVerbose(verbose);
        its.push_back(iterator);
      }

      
      // +--------------------+
      // | Feature extraction |
      // +--------------------+
      else if (cmdName == "ExtractFeature") {
        bool ignoreStrand    = ApplicationTools::getBooleanParameter("ignore_strand", cmdArgs, false);
        bool completeOnly    = ApplicationTools::getBooleanParameter("complete", cmdArgs, false);
        string refSpecies    = ApplicationTools::getStringParameter("ref_species", cmdArgs, "none");
        string featureFile   = ApplicationTools::getAFilePath("feature.file", cmdArgs, false, false);
        string featureFormat = ApplicationTools::getStringParameter("feature.format", cmdArgs, "GFF");
        vector<string> featureType = ApplicationTools::getVectorParameter<string>("feature.type", cmdArgs, ',', "all");
        if (featureType.size() == 0)
          throw Exception("At least one feature should be provided for command 'ExtractFeature'.");
        ApplicationTools::displayResult("-- Features to extract", featureFile + " (" + featureFormat + ")");
        ApplicationTools::displayResult("-- Features are for species", refSpecies);
        ApplicationTools::displayBooleanResult("-- Features are strand-aware", !ignoreStrand);
        ApplicationTools::displayBooleanResult("-- Extract incomplete features", !completeOnly);
        compress = ApplicationTools::getStringParameter("feature.file.compression", cmdArgs, "none");
        filtering_istream featureStream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          featureStream.push(gzip_decompressor());
        } else if (compress == "zip") {
          featureStream.push(zlib_decompressor());
        } else if (compress == "bzip2") {
          featureStream.push(bzip2_decompressor());
        } else
          throw Exception("Bad input incompression format: " + compress);
        featureStream.push(file_source(featureFile));
        unique_ptr<FeatureReader> ftReader;
        SequenceFeatureSet featuresSet;
        if (featureFormat == "GFF") {
          ftReader.reset(new GffFeatureReader(featureStream));
        } else if (featureFormat == "GTF") {
          ftReader.reset(new GtfFeatureReader(featureStream));
        } else if (featureFormat == "BedGraph") {
          ftReader.reset(new BedGraphFeatureReader(featureStream));
         } else
          throw Exception("Unsupported feature format: " + featureFormat);
        if (featureType.size() == 1 && featureType[0] == "all")
          ftReader->getAllFeatures(featuresSet);
        else {
          for (size_t i = 0; i < featureType.size(); ++i) {
            ApplicationTools::displayResult("-- Extract features of type", featureType[i]);
            ftReader->getFeaturesOfType(featureType[i], featuresSet);
          }
        }
        ApplicationTools::displayResult("-- Total number of features", featuresSet.getNumberOfFeatures());
        FeatureExtractorMafIterator* iterator = new FeatureExtractorMafIterator(currentIterator, refSpecies, featuresSet, completeOnly, ignoreStrand);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        its.push_back(iterator);

        currentIterator = iterator;
      }



      // +------------------+
      // | Window splitting |
      // +------------------+
      else if (cmdName == "WindowSplit") {
        if (cmdArgs.find("preferred.size") != cmdArgs.end()) {
          throw Exception("preferred.size argument in WindowSplit is deprecated: use preferred_size instead.");
        }
        unsigned int preferredSize = ApplicationTools::getParameter<unsigned int>("preferred_size", cmdArgs, 0);
        ApplicationTools::displayResult("-- Preferred size", preferredSize);
        string splitOptionStr = ApplicationTools::getStringParameter("align", cmdArgs, "center");
        short splitOption;
        if (splitOptionStr == "ragged_left")
          splitOption = WindowSplitMafIterator::RAGGED_LEFT;
        else if (splitOptionStr == "ragged_right")
          splitOption = WindowSplitMafIterator::RAGGED_RIGHT;
        else if (splitOptionStr == "center")
          splitOption = WindowSplitMafIterator::CENTER;
        else if (splitOptionStr == "adjust")
          splitOption = WindowSplitMafIterator::ADJUST;
        else throw Exception("Unvalid alignment option for WindowSplit: " + splitOptionStr);
        ApplicationTools::displayResult("-- Alignment option", splitOptionStr);
        bool keepSmallBlocks = ApplicationTools::getBooleanParameter("keep_small_blocks", cmdArgs, false);
        if (splitOptionStr == "adjust")
          ApplicationTools::displayBooleanResult("-- Keep small blocks", keepSmallBlocks);

        WindowSplitMafIterator* iterator = new WindowSplitMafIterator(currentIterator, preferredSize, splitOption, keepSmallBlocks);
        iterator->setLogStream(log);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +---------------------+
      // | Distance estimation |
      // +---------------------+
      else if (cmdName == "DistanceEstimation") {
        string distMethod = ApplicationTools::getStringParameter("method", cmdArgs, "count");
        ApplicationTools::displayResult("-- Method", distMethod);
        if (distMethod == "count") {
          string gapOption = ApplicationTools::getStringParameter("gap_option", cmdArgs, "no_gap");
          if (gapOption == "all") {
            gapOption = SiteContainerTools::SIMILARITY_ALL;
          } else if (gapOption == "no_gap") {
            gapOption = SiteContainerTools::SIMILARITY_NOGAP;
          } else if (gapOption == "no_full_gap") {
            gapOption = SiteContainerTools::SIMILARITY_NOFULLGAP;
          } else if (gapOption == "no_double_gap") {
            gapOption = SiteContainerTools::SIMILARITY_NODOUBLEGAP;
          } else {
            throw Exception("Unrecognized gap option, should be either 'all', 'no_full_gap', 'no_double_gap' or 'no_gap'.");
          }
          ApplicationTools::displayResult("-- Gap option", gapOption);
          bool unresolvedAsGap = ApplicationTools::getBooleanParameter("unresolved_as_gap", cmdArgs, "no");
          ApplicationTools::displayBooleanResult("-- Unresolved as gaps", unresolvedAsGap);
          bool extendedSeqNames = ApplicationTools::getBooleanParameter("extended_names", cmdArgs, true);
          ApplicationTools::displayBooleanResult("-- Use extended names in matrix", extendedSeqNames);

          CountDistanceEstimationMafIterator* iterator = new CountDistanceEstimationMafIterator(currentIterator, gapOption, unresolvedAsGap, extendedSeqNames);
          ApplicationTools::displayResult("-- Block-wise matrices are registered as", iterator->getPropertyName());
          iterator->setLogStream(log);
          currentIterator = iterator;
          its.push_back(iterator);
        } else if (distMethod == "ml") {
          string modelDesc = ApplicationTools::getStringParameter("model", cmdArgs, "JC()");
          string rdistDesc = ApplicationTools::getStringParameter("rate", cmdArgs, "Constant()");
          string paramOpt  = ApplicationTools::getStringParameter("parameter_estimation", cmdArgs, "initial");
          if (paramOpt == "initial") {
            paramOpt = OptimizationTools::DISTANCEMETHOD_INIT;
          } else if (paramOpt == "pairwise") {
            paramOpt = OptimizationTools::DISTANCEMETHOD_PAIRWISE;
          } else {
            throw Exception("Unrecognized parameter option, should be either 'initial', 'pairwise'.");
          }
          string prPath = ApplicationTools::getAFilePath("profiler", cmdArgs, false, false);
          string mhPath = ApplicationTools::getAFilePath("message_handler", cmdArgs, false, false);
          double propGapsToKeep = ApplicationTools::getDoubleParameter("max_freq_gaps", cmdArgs, 0.);
          bool gapsAsUnresolved = ApplicationTools::getBooleanParameter("gaps_as_unresolved", cmdArgs, true);
          
          ApplicationTools::displayResult("-- Max. frequency of gaps", propGapsToKeep);
          ApplicationTools::displayBooleanResult("-- Gaps as unresolved", gapsAsUnresolved);
          
          bool extendedSeqNames = ApplicationTools::getBooleanParameter("extended_names", cmdArgs, true);
          ApplicationTools::displayBooleanResult("-- Use extended names in matrix", extendedSeqNames);
          
          BppOSubstitutionModelFormat modelReader(BppOSubstitutionModelFormat::DNA, false, false, true, true, 1);
          unique_ptr<SubstitutionModel> model(modelReader.read(&AlphabetTools::DNA_ALPHABET, modelDesc, 0, true));
          BppORateDistributionFormat rdistReader(true);
          unique_ptr<DiscreteDistribution> rdist(rdistReader.read(rdistDesc, true)); 
          unique_ptr<DistanceEstimation> distEst(new DistanceEstimation(model.release(), rdist.release()));
          
          OutputStream* profiler =
            (prPath == "none") ? 0 :
              (prPath == "std") ? ApplicationTools::message.get() :
              new StlOutputStream(new ofstream(prPath.c_str(), ios::out));
          if (profiler)
            profiler->setPrecision(20);
          if (verbose)
            ApplicationTools::displayResult("-- Optimization profile in", prPath);
          distEst->getOptimizer()->setProfiler(profiler);

          OutputStream* messenger =
            (mhPath == "none") ? 0 :
              (mhPath == "std") ? ApplicationTools::message.get() :
              new StlOutputStream(new ofstream(mhPath.c_str(), ios::out));
          if (messenger)
            messenger->setPrecision(20);
          if (verbose)
            ApplicationTools::displayResult("-- Optimization messages in", mhPath);
          distEst->getOptimizer()->setMessageHandler(messenger);

          MaximumLikelihoodDistanceEstimationMafIterator* iterator = new MaximumLikelihoodDistanceEstimationMafIterator(currentIterator,
              distEst.release(), propGapsToKeep, gapsAsUnresolved, paramOpt, extendedSeqNames);
          ApplicationTools::displayResult("-- Block-wise matrices are registered as", iterator->getPropertyName());
          iterator->setLogStream(log);
          iterator->setVerbose(verbose);
          currentIterator = iterator;
          its.push_back(iterator);
        } else {
          throw Exception("Unknown distance method: " + distMethod);
        }
      }



      // +--------------------------+
      // | Phylogeny reconstruction |
      // +--------------------------+
      else if (cmdName == "DistanceBasedPhylogeny") {
        string distMethodName = ApplicationTools::getStringParameter("method", cmdArgs, "bionj");
        string distProperty = ApplicationTools::getStringParameter("dist_mat", cmdArgs, "none");
        DistanceMethod* distMethod = 0;
        if (distMethodName == "upgma") {
          distMethod = new PGMA(false);
        } else if (distMethodName == "wpgma") {
          distMethod = new PGMA(true);
        } else if (distMethodName == "nj") {
          distMethod = new NeighborJoining(false, false);
        } else if (distMethodName == "bionj") {
          distMethod = new BioNJ(false, false);
        } else {
          throw Exception("Unknown distance-based phylogenetic method: " + distMethodName); 
        }
        distMethod->setVerbose(false);
        ApplicationTools::displayResult("-- Reading distance matrix from", distProperty);
        ApplicationTools::displayResult("-- Build distance tree using", distMethodName);

        DistanceBasedPhylogenyReconstructionMafIterator* iterator = new DistanceBasedPhylogenyReconstructionMafIterator(currentIterator, distMethod, distProperty);
        ApplicationTools::displayResult("-- Writing block-wise trees to", iterator->getPropertyName());
        iterator->setLogStream(log);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +-----------------------------------+
      // | External phylogeny reconstruction |
      // +-----------------------------------+
      else if (cmdName == "ExternalTreeBuilding") {
        string name = ApplicationTools::getStringParameter("name", cmdArgs, "external");
        
        string programInputFile = ApplicationTools::getAFilePath("input.file", cmdArgs, true, false);
        string programInputFormat = ApplicationTools::getStringParameter("input.format", cmdArgs, "Fasta");
        BppOAlignmentWriterFormat bppoWriter(1);
        OAlignment* alnWriter(bppoWriter.read(programInputFormat));

        string programOutputFile = ApplicationTools::getAFilePath("output.file", cmdArgs, true, false);
        string programOutputFormat = ApplicationTools::getStringParameter("output.format", cmdArgs, "Newick");
        BppOTreeReaderFormat bppoReader(1);
        ITree* treeReader(bppoReader.read(programOutputFormat));

        string propertyName = ApplicationTools::getStringParameter("property_name", cmdArgs, "ExternalTree");
        ApplicationTools::displayResult("-- Registering block-wise trees to", propertyName);

        string command = ApplicationTools::getStringParameter("call", cmdArgs, "echo \"TODO: implement wrapper!\"");
        
        ApplicationTools::displayResult("-- External call (tree building)", name);
        ApplicationTools::displayResult("   Command", command);

        TreeBuildingSystemCallMafIterator* iterator = new TreeBuildingSystemCallMafIterator(currentIterator, alnWriter, programInputFile, treeReader, programOutputFile, command, propertyName);

        iterator->setLogStream(log);
        currentIterator = iterator;
        its.push_back(iterator);
      }




      // +-------------------+
      // | Phylogeny rooting |
      // +-------------------+
      else if (cmdName == "NewOutgroup") {
        string treePropertyInput = ApplicationTools::getStringParameter("tree_input", cmdArgs, "none");
        string treePropertyOutput = ApplicationTools::getStringParameter("tree_output", cmdArgs, "none");
        string outgroup = ApplicationTools::getStringParameter("outgroup", cmdArgs, "none");
        ApplicationTools::displayResult("-- Reading tree from", treePropertyInput);
        ApplicationTools::displayResult("-- Rerooting according to species", outgroup);
        NewOutgroupMafIterator* iterator = new NewOutgroupMafIterator(currentIterator, treePropertyInput, treePropertyOutput, outgroup);
        ApplicationTools::displayResult("-- Writing tree to", treePropertyOutput);
        iterator->setLogStream(log);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +------------------------+
      // | Phylogeny drop species |
      // +------------------------+
      else if (cmdName == "DropSpecies") {
        string treePropertyInput = ApplicationTools::getStringParameter("tree_input", cmdArgs, "none");
        string treePropertyOutput = ApplicationTools::getStringParameter("tree_output", cmdArgs, "none");
        string species = ApplicationTools::getStringParameter("species", cmdArgs, "none");
        ApplicationTools::displayResult("-- Reading tree from", treePropertyInput);
        ApplicationTools::displayResult("-- Removing leaves from species", species);
        DropSpeciesMafIterator* iterator = new DropSpeciesMafIterator(currentIterator, treePropertyInput, treePropertyOutput, species);
        ApplicationTools::displayResult("-- Writing tree to", treePropertyOutput);
        iterator->setLogStream(log);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +--------+
      // | Output |
      // +--------+
      else if (cmdName == "Output") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output file", outputFile);
        filtering_ostream* out = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          out->push(gzip_compressor());
        } else if (compress == "zip") {
          out->push(zlib_compressor());
        } else if (compress == "bzip2") {
          out->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        out->push(file_sink(outputFile));
        ostreams.push_back(out);
        ApplicationTools::displayResult("-- File compression", compress);
        bool mask = ApplicationTools::getBooleanParameter("mask", cmdArgs, true);
        ApplicationTools::displayBooleanResult("-- Output mask", mask);
        OutputMafIterator* iterator = new OutputMafIterator(currentIterator, out, mask);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +-------------------+
      // | Output alignments |
      // +-------------------+
      else if (cmdName == "OutputAlignments") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        bool multipleFiles = (outputFile.find("%i") != string::npos);
        ApplicationTools::displayResult("-- Output alignment file" + string(multipleFiles ? "s" : ""), outputFile);
        bool mask = ApplicationTools::getBooleanParameter("mask", cmdArgs, true);
        ApplicationTools::displayBooleanResult("-- Output mask", mask);
        bool coords = ApplicationTools::getBooleanParameter("coordinates", cmdArgs, true);
        ApplicationTools::displayBooleanResult("-- Output coordinates", coords);
        bool header = ApplicationTools::getBooleanParameter("ldhat_header", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Output header line", header);
         string reference = ApplicationTools::getStringParameter("reference", cmdArgs, "", "", true, 1);
        if (reference != "")
          ApplicationTools::displayResult("-- Reference species", reference);
        
        OutputAlignmentMafIterator* iterator; 
        BppOAlignmentWriterFormat bppoWriter(1);
        string description = ApplicationTools::getStringParameter("format", cmdArgs, "Clustal");
        OAlignment* oAln = bppoWriter.read(description);
        if (multipleFiles) {
          iterator = new OutputAlignmentMafIterator(currentIterator, outputFile, oAln, mask, coords, header, reference);
        } else {
          compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
          filtering_ostream* out = new filtering_ostream;
          if (compress == "none") {
          } else if (compress == "gzip") {
            out->push(gzip_compressor());
          } else if (compress == "zip") {
            out->push(zlib_compressor());
          } else if (compress == "bzip2") {
            out->push(bzip2_compressor());
          } else
            throw Exception("Bad output compression format: " + compress);
          out->push(file_sink(outputFile));
          ostreams.push_back(out);
          ApplicationTools::displayResult("-- File compression", compress);
          iterator = new OutputAlignmentMafIterator(currentIterator, out, oAln, mask, coords, header, reference);
        }
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +--------------------+
      // | Output as features |
      // +--------------------+
      else if (cmdName == "OutputAsFeatures") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output feature file", outputFile);
        filtering_ostream* out = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          out->push(gzip_compressor());
        } else if (compress == "zip") {
          out->push(zlib_compressor());
        } else if (compress == "bzip2") {
          out->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        out->push(file_sink(outputFile));
        ostreams.push_back(out);
        ApplicationTools::displayResult("-- File compression", compress);
        string species = ApplicationTools::getStringParameter("species", cmdArgs, "");
        if (species == "")
          throw Exception("A species name should be provided for command 'OutputAsFeatures'.");
        ApplicationTools::displayResult("-- Species to use", species);

        OutputAsFeaturesMafIterator* iterator = new OutputAsFeaturesMafIterator(currentIterator, out, species);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +-----------------+
      // | Output as table |
      // +-----------------+
      else if (cmdName == "OutputAsTable") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output table file", outputFile);
        filtering_ostream* out = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          out->push(gzip_compressor());
        } else if (compress == "zip") {
          out->push(zlib_compressor());
        } else if (compress == "bzip2") {
          out->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        out->push(file_sink(outputFile));
        ostreams.push_back(out);
        ApplicationTools::displayResult("-- File compression", compress);
        string reference = ApplicationTools::getStringParameter("reference", cmdArgs, "");
        if (reference != "")
        ApplicationTools::displayResult("-- Reference sequence", reference);
 
        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "", "", false, false);

        TableOutputMafIterator* iterator = new TableOutputMafIterator(currentIterator, out, species, reference);
        currentIterator = iterator;
        its.push_back(iterator);
      }




      // +------------+
      // | VCF output |
      // +------------+
      else if (cmdName == "VcfOutput") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output file", outputFile);
        filtering_ostream* out = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          out->push(gzip_compressor());
        } else if (compress == "zip") {
          out->push(zlib_compressor());
        } else if (compress == "bzip2") {
          out->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        out->push(file_sink(outputFile));
        ostreams.push_back(out);
        ApplicationTools::displayResult("-- File compression", compress);

        string reference = ApplicationTools::getStringParameter("reference", cmdArgs, "");
        if (reference == "")
          throw Exception("A reference sequence should be provided for filter 'VcfOutput'.");
        ApplicationTools::displayResult("-- Reference sequence", reference);
        
        vector<string> genotypes = ApplicationTools::getVectorParameter<string>("genotypes", cmdArgs, ',', "");
        for (size_t i = 0; i < genotypes.size(); ++i) {
          ApplicationTools::displayResult("-- Adding genotype info for", genotypes[i]);
        }
        
        bool outputAll = ApplicationTools::getBooleanParameter("all", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Output non-variable positions", outputAll);

        bool outputDiploids = ApplicationTools::getBooleanParameter("diploids", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Output (homozygous) diploids", outputDiploids);

        VcfOutputMafIterator* iterator = new VcfOutputMafIterator(currentIterator, out, reference, genotypes, outputAll, outputDiploids);

        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +-------------+
      // | MSMC output |
      // +-------------+
      else if (cmdName == "MsmcOutput") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output file", outputFile);
        filtering_ostream* out = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          out->push(gzip_compressor());
        } else if (compress == "zip") {
          out->push(zlib_compressor());
        } else if (compress == "bzip2") {
          out->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        out->push(file_sink(outputFile));
        ostreams.push_back(out);
        ApplicationTools::displayResult("-- File compression", compress);

        string reference = ApplicationTools::getStringParameter("reference", cmdArgs, "");
        if (reference == "")
          throw Exception("A reference sequence should be provided for filter 'MsmcOutput'.");
        ApplicationTools::displayResult("-- Reference sequence", reference);
        
        vector<string> species = ApplicationTools::getVectorParameter<string>("genotypes", cmdArgs, ',', "");
        if (species.size() < 2)
          throw Exception("MsmcOutput: at least two genomes are necessary to call SNPs.");
        MsmcOutputMafIterator* iterator = new MsmcOutputMafIterator(currentIterator, out, species, reference);

        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +--------------+
      // | PLINK output |
      // +--------------+
      else if (cmdName == "PlinkOutput") {
        string outputPedFile = ApplicationTools::getAFilePath("ped_file", cmdArgs, true, false);
        string outputMapFile = ApplicationTools::getAFilePath("map_file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output Ped file", outputPedFile);
        ApplicationTools::displayResult("-- Output Map file", outputMapFile);
        filtering_ostream* outPed = new filtering_ostream;
        filtering_ostream* outMap = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          outPed->push(gzip_compressor());
          outMap->push(gzip_compressor());
        } else if (compress == "zip") {
          outPed->push(zlib_compressor());
          outMap->push(zlib_compressor());
        } else if (compress == "bzip2") {
          outPed->push(bzip2_compressor());
          outMap->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        outPed->push(file_sink(outputPedFile));
        outMap->push(file_sink(outputMapFile));
        ostreams.push_back(outPed);
        ostreams.push_back(outMap);
        ApplicationTools::displayResult("-- File compression", compress);

        string reference = ApplicationTools::getStringParameter("reference", cmdArgs, "");
        if (reference == "")
          throw Exception("A reference sequence should be provided for filter 'PlinkOutput'.");
        ApplicationTools::displayResult("-- Reference sequence", reference);
        
        bool map3 = ApplicationTools::getBooleanParameter("map3", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Output map3 file", map3);

        bool recodeChr = ApplicationTools::getBooleanParameter("recode_chr", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Recode chromosomes", recodeChr);

        vector<string> species = ApplicationTools::getVectorParameter<string>("genotypes", cmdArgs, ',', "");
        if (species.size() < 2)
          throw Exception("PlinkOutput: at least two genomes are necessary to call SNPs.");

        PlinkOutputMafIterator* iterator = new PlinkOutputMafIterator(currentIterator, outPed, outMap, species, reference, map3, recodeChr);

        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +----------------------+
      // | SequenceLDhot output |
      // +----------------------+
      else if (cmdName == "SequenceLDhotOutput") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        ApplicationTools::displayResult("-- Output file", outputFile);

        string reference = ApplicationTools::getStringParameter("reference", cmdArgs, "");
        if (!TextTools::isEmpty(reference))
          ApplicationTools::displayResult("-- Reference sequence", reference);
        else
          reference = "";
        
        bool completeOnly = ApplicationTools::getBooleanParameter("complete_only", cmdArgs, true);
        ApplicationTools::displayBooleanResult("-- Use only complete sites", completeOnly);

        SequenceLDhotOutputMafIterator* iterator = new SequenceLDhotOutputMafIterator(currentIterator, outputFile, completeOnly, reference);

        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }



      // +--------------------+
      // | Coordinates output |
      // +--------------------+
      else if (cmdName == "OutputCoordinates") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output file", outputFile);
        filtering_ostream* out = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          out->push(gzip_compressor());
        } else if (compress == "zip") {
          out->push(zlib_compressor());
        } else if (compress == "bzip2") {
          out->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        out->push(file_sink(outputFile));
        ostreams.push_back(out);
        ApplicationTools::displayResult("-- File compression", compress);

        vector<string> species = ApplicationTools::getVectorParameter<string>("species", cmdArgs, ',', "");
        if (species.size() == 0)
          throw Exception("At least one species should be provided for filter 'OutputCoordinates'.");
        ApplicationTools::displayResult("-- Output coordinates for", TextTools::toString(species.size()) + " species");
        
        bool includeSrcSize = ApplicationTools::getBooleanParameter("output_src_size", cmdArgs, true);
        ApplicationTools::displayBooleanResult("-- Output src size", includeSrcSize);
        
        for (size_t i = 0; i < species.size(); ++i) {
          ApplicationTools::displayResult("-- Output coordinates for species", species[i]);
        }
        CoordinatesOutputMafIterator* iterator = new CoordinatesOutputMafIterator(currentIterator, out, species, includeSrcSize);

        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }




      // +------------------------+
      // | Coordinates conversion |
      // +------------------------+
      else if (cmdName == "LiftOver") {
        //Input
        string refSpecies    = ApplicationTools::getStringParameter("ref_species", cmdArgs, "none");
        string targetSpecies = ApplicationTools::getStringParameter("target_species", cmdArgs, "none");
        string featureFile   = ApplicationTools::getAFilePath("feature.file", cmdArgs, false, false);
        string featureFormat = ApplicationTools::getStringParameter("feature.format", cmdArgs, "GFF");
        bool outputClosest   = ApplicationTools::getBooleanParameter("target_closest_position", cmdArgs, true);
        ApplicationTools::displayResult("-- Features to lift over", featureFile + " (" + featureFormat + ")");
        ApplicationTools::displayResult("-- from species", refSpecies);
        ApplicationTools::displayResult("-- to species", targetSpecies);
        ApplicationTools::displayBooleanResult("-- closest position if gap", outputClosest);
        compress = ApplicationTools::getStringParameter("feature.file.compression", cmdArgs, "none");
        filtering_istream featureStream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          featureStream.push(gzip_decompressor());
        } else if (compress == "zip") {
          featureStream.push(zlib_decompressor());
        } else if (compress == "bzip2") {
          featureStream.push(bzip2_decompressor());
        } else
          throw Exception("Bad input incompression format: " + compress);
        featureStream.push(file_source(featureFile));
        unique_ptr<FeatureReader> ftReader;
        SequenceFeatureSet featuresSet;
        if (featureFormat == "GFF") {
          ftReader.reset(new GffFeatureReader(featureStream));
        } else if (featureFormat == "GTF") {
          ftReader.reset(new GtfFeatureReader(featureStream));
        } else if (featureFormat == "BedGraph") {
          ftReader.reset(new BedGraphFeatureReader(featureStream));
        } else
          throw Exception("Unsupported feature format: " + featureFormat);
        ftReader->getAllFeatures(featuresSet);
        ApplicationTools::displayResult("-- Total number of features", featuresSet.getNumberOfFeatures());
        
        //Output
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output file", outputFile);
        filtering_ostream* out = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          out->push(gzip_compressor());
        } else if (compress == "zip") {
          out->push(zlib_compressor());
        } else if (compress == "bzip2") {
          out->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        out->push(file_sink(outputFile));
        ostreams.push_back(out);
        ApplicationTools::displayResult("-- File compression", compress);
        
        //Iterator initialization:
        CoordinateTranslatorMafIterator* iterator = new CoordinateTranslatorMafIterator(currentIterator, refSpecies, targetSpecies, featuresSet, *out, outputClosest);
        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        its.push_back(iterator);

        currentIterator = iterator;
      }



      // +--------------+
      // | Output trees |
      // +--------------+
      else if (cmdName == "OutputTrees") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output tree file", outputFile);
        filtering_ostream* out = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          out->push(gzip_compressor());
        } else if (compress == "zip") {
          out->push(zlib_compressor());
        } else if (compress == "bzip2") {
          out->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        out->push(file_sink(outputFile));
        ostreams.push_back(out);
        ApplicationTools::displayResult("-- File compression", compress);
        string treeProperty = ApplicationTools::getStringParameter("tree", cmdArgs, "none");
        ApplicationTools::displayResult("-- Tree to write", treeProperty);
        bool stripNames = ApplicationTools::getBooleanParameter("strip_names", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Strip names", stripNames);

        OutputTreeMafIterator* iterator = new OutputTreeMafIterator(currentIterator, out, treeProperty, !stripNames);
        currentIterator = iterator;
        its.push_back(iterator);
      }
    



      // +--------------------------+
      // | Output Distance matrices |
      // +--------------------------+
      else if (cmdName == "OutputDistanceMatrices") {
        string outputFile = ApplicationTools::getAFilePath("file", cmdArgs, true, false);
        compress = ApplicationTools::getStringParameter("compression", cmdArgs, "none");
        ApplicationTools::displayResult("-- Output matrix file", outputFile);
        filtering_ostream* out = new filtering_ostream;
        if (compress == "none") {
        } else if (compress == "gzip") {
          out->push(gzip_compressor());
        } else if (compress == "zip") {
          out->push(zlib_compressor());
        } else if (compress == "bzip2") {
          out->push(bzip2_compressor());
        } else
          throw Exception("Bad output compression format: " + compress);
        out->push(file_sink(outputFile));
        ostreams.push_back(out);
        ApplicationTools::displayResult("-- File compression", compress);
        string distProperty = ApplicationTools::getStringParameter("distance", cmdArgs, "none");
        ApplicationTools::displayResult("-- Matrix to write", distProperty);
        bool stripNames = ApplicationTools::getBooleanParameter("strip_names", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Strip names", stripNames);

        OutputDistanceMatrixMafIterator* iterator = new OutputDistanceMatrixMafIterator(currentIterator, out, distProperty, !stripNames);
        currentIterator = iterator;
        its.push_back(iterator);
      }
    



      // +--------------------------+
      // | External program wrapper |
      // +--------------------------+
      else if (cmdName == "SystemCall") {
        string name = ApplicationTools::getStringParameter("name", cmdArgs, "external");

        string programInputFile = ApplicationTools::getAFilePath("input.file", cmdArgs, true, false);
        string programInputFormat = ApplicationTools::getStringParameter("input.format", cmdArgs, "Fasta");
        BppOAlignmentWriterFormat bppoWriter(1);
        OAlignment* alnWriter(bppoWriter.read(programInputFormat));

        string programOutputFile = ApplicationTools::getAFilePath("output.file", cmdArgs, true, false);
        string programOutputFormat = ApplicationTools::getStringParameter("output.format", cmdArgs, "Fasta");
        BppOAlignmentReaderFormat bppoReader(1);
        IAlignment* alnReader(bppoReader.read(programOutputFormat));
        
        bool hotTest = ApplicationTools::getBooleanParameter("hot", cmdArgs, false);
        ApplicationTools::displayBooleanResult("-- Compute HoT score", hotTest);

        string command = ApplicationTools::getStringParameter("call", cmdArgs, "echo \"TODO: implement wrapper!\"");
        
        ApplicationTools::displayResult("-- External call", name);
        ApplicationTools::displayResult("   Command", command);

        SystemCallMafIterator* iterator = new SystemCallMafIterator(currentIterator, 
            alnWriter, programInputFile, alnReader, programOutputFile, command, hotTest);

        iterator->setLogStream(log);
        iterator->setVerbose(verbose);
        currentIterator = iterator;
        its.push_back(iterator);
      }


      else 
        throw Exception("Unknown filter: " + cmdName);
    }

    //Now loop over the last iterator and that's it!
    size_t blockCounter = 0;
    size_t alnSize = 0;
    cout << "Parsing..." << endl;
    while (MafBlock* block = currentIterator->nextBlock())
    {
      alnSize += block->getNumberOfSites();
      cout << '\r' << ++blockCounter << " blocks kept, totalizing " << alnSize << "bp.";
      cout.flush();
      //ApplicationTools::displayUnlimitedGauge(blockCounter++, "Parsing...");
      delete block;
    }
    ApplicationTools::message->endLine();

    //Flush all streams:
    for (size_t i = 0; i < ostreams.size(); ++i) {
      close(*ostreams[i]);
    }

    //Clean memory:
    for (size_t i = 0; i < its.size(); ++i) {
      delete its[i];
    }

    maffilter.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}

