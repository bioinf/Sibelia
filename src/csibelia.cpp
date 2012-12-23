//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include <tclap/CmdLine.h>
#include "util.h"
#include "variantcalling/variantcaller.h"

int main(int argc, char * argv[])
{	/*
	signal(SIGINT, SignalHandler);
	signal(SIGABRT, SignalHandler);	
	signal(SIGTERM, SignalHandler);

	std::stringstream parsets;		
	const std::string parameterSetNameArray[] = {"loose", "fine"};
	std::vector<std::string> parameterSetName(parameterSetNameArray, parameterSetNameArray + sizeof(parameterSetNameArray) / sizeof(parameterSetNameArray[0]));
	std::map<std::string, std::vector<std::pair<int, int> > > defaultParameters;
	defaultParameters["loose"] = LooseStageFile();
	defaultParameters["fine"] = FineStageFile();

	try
	{  
		TCLAP::CmdLine cmd("Program for finding syteny blocks in closely related genomes", ' ', "2.5.0");
		TCLAP::ValueArg<unsigned int> maxIterations("i",
			"maxiterations",
			"Maximum number of iterations during a stage of simplification, default = 4.",
			false,
			4,
			"integer",
			cmd);

		TCLAP::ValueArg<std::string> stageFile("k",
			"stagefile",
			"File that contains manually chosen simplifications parameters. See USAGE file for more information.",
			false,
			"",
			"file name");

		std::string description = std::string("Parameters set, used for the simplification. ") +
			std::string("Option \"loose\" produces fewer blocks, but they are larger (\"fine\" is opposite).");
		TCLAP::ValuesConstraint<std::string> allowedParametersVals(parameterSetName);
		TCLAP::ValueArg<std::string> parameters("s",
			"parameters",
			description,
			false,
			parameterSetName[0],
			&allowedParametersVals);

		TCLAP::ValueArg<unsigned int> minBlockSize("m",
			"minblocksize",
			"Minimum size of a synteny block, default value = 5000 BP.",
			false,
			5000,
			"integer",
			cmd);

		TCLAP::ValueArg<std::string> outFileDir("o",
			"outdir",
			"Directory where output files are written",
			false,
			".",
			"dir name",
			cmd);
		
		TCLAP::ValueArg<std::string> referenceFile("r",
			"reference",
			"FASTA file with the reference genome.",
			true,
			".",
			"file name",
			cmd);

		TCLAP::ValueArg<std::string> assemblyFile("a",
			"assembly",
			"FASTA file with the assembly.",
			true,
			".",
			"file name",
			cmd);

		cmd.xorAdd(parameters, stageFile);
		cmd.parse(argc, argv);

		std::vector<std::pair<int, int> > stage;
		if(parameters.isSet())
		{
			stage = defaultParameters[parameters.getValue()];
		}
		else
		{
			stage = ReadStageFile(stageFile.getValue());
		}
		
		size_t totalSize = 0;
		std::vector<SyntenyFinder::FASTARecord> chrList;
		std::string fileName[] = {referenceFile.getValue(), assemblyFile.getValue()};		
		for(size_t file = 0; file < 2; file++)
		{
			SyntenyFinder::FASTAReader reader(fileName[file]);
			if(!reader.IsOk())
			{
				throw std::runtime_error(("Cannot open file " + fileName[file]).c_str());
			}

			size_t seqCount = reader.GetSequences(chrList);			
			if(file == 0 && seqCount != 1)
			{
				throw std::runtime_error(("File with the reference must contain exactly one sequence " + fileName[file]).c_str());
			}
		}
		
		for(size_t i = 0; i < chrList.size(); i++)
		{
			totalSize += chrList[i].GetSequence().size();
		}

		if(totalSize > SyntenyFinder::MAX_INPUT_SIZE)
		{
			throw std::runtime_error("Input is larger than 1 GB, can't proceed");
		}

		int trimK = INT_MAX;
		SyntenyFinder::BlockFinder finder(chrList);
		for(size_t i = 0; i < stage.size(); i++)
		{
			trimK = std::min(trimK, stage[i].first);
			std::cout << "Simplification stage " << i + 1 << " of " << stage.size() << std::endl;
			std::cout << "Enumerating vertices of the graph, then performing bulge removal..." << std::endl;
			finder.PerformGraphSimplifications(stage[i].first, stage[i].second, maxIterations.getValue(), PutProgressChr);
		}
		
		std::vector<SyntenyFinder::BlockInstance> blockList;
		std::cout << "Finding variants and generating the output..." << std::endl;
		size_t lastK = std::min(stage.back().first, static_cast<int>(minBlockSize.getValue()));
		trimK = std::min(trimK, static_cast<int>(minBlockSize.getValue()));
		finder.GenerateSyntenyBlocks(lastK, trimK, minBlockSize.getValue(), blockList, false, PutProgressChr);
		size_t refSeqId = chrList[0].GetId();
		
		std::vector<SyntenyFinder::Variant> variant;
		SyntenyFinder::VariantCaller caller(refSeqId, blockList, trimK);
	//	caller.CallVariants(variant);

		SyntenyFinder::OutputGenerator generator(chrList, blockList);
		SyntenyFinder::CreateDirectory(outFileDir.getValue());
		const std::string defaultCoordsFile = outFileDir.getValue() + "/block_coords.txt";
		const std::string defaultPermutationsFile = outFileDir.getValue() + "/genomes_permutations.txt";
		const std::string defaultCoverageReportFile = outFileDir.getValue() + "/coverage_report.txt";
		const std::string defaultSequencesFile = outFileDir.getValue() + "/blocks_sequences.fasta";
		const std::string defaultGraphFile = outFileDir.getValue() + "/de_bruijn_graph.dot";

		const std::string defaultCircosDir = outFileDir.getValue() + "/circos";
		const std::string defaultCircosFile = defaultCircosDir + "/circos.conf";
		const std::string defaultD3File = outFileDir.getValue() + "/d3_blocks_diagram.html";
		const std::string defaultVariantFile = outFileDir.getValue() + "/variant.txt";


		std::string buf;
		std::ofstream variantFile(defaultVariantFile.c_str());
		for(size_t i = 0; i < variant.size(); i++)
		{
			variant[i].ToString(buf);
			variantFile << buf << std::endl;
		}

		boost::function<void()> outFunction[] = 
		{
			boost::bind(&SyntenyFinder::OutputGenerator::ListChromosomesAsPermutations, boost::cref(generator), defaultPermutationsFile),
			boost::bind(&SyntenyFinder::OutputGenerator::GenerateReport, boost::cref(generator), defaultCoverageReportFile),
			boost::bind(&SyntenyFinder::OutputGenerator::ListBlocksIndices, boost::cref(generator), defaultCoordsFile),
			boost::bind(&SyntenyFinder::OutputGenerator::ListBlocksSequences, boost::cref(generator), defaultSequencesFile),		
			boost::bind(&SyntenyFinder::OutputGenerator::GenerateCircosOutput, boost::cref(generator), defaultCircosFile, defaultCircosDir),
			boost::bind(&SyntenyFinder::OutputGenerator::GenerateD3Output, boost::cref(generator), defaultD3File)
		};
		
		size_t length = sizeof(outFunction) / sizeof(outFunction[0]);
		for(size_t i = 0; i < length; i++)
		{
			outFunction[i]();
		}

		std::stringstream buffer;
		finder.SerializeCondensedGraph(lastK, buffer, PutProgressChr);
		generator.OutputBuffer(defaultGraphFile, buffer.str());
	} 
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
	}
	catch (std::runtime_error & e)
	{
		std::cerr << "error: " << e.what() << std::endl;
	}
	catch(...)
	{
		SyntenyFinder::TempFile::Cleanup();
	}
	*/
	return 0;
}
