//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include <tclap/CmdLine.h>
#include "util.h"
#include "variantcalling/variantcaller.h"
#include "variantcalling/simplerearrangements.h"

int main(int argc, char * argv[])
{	
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
		TCLAP::CmdLine cmd("Program for finding syteny blocks in closely related genomes", ' ', "3.0.0");
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
			"Minimum size of a synteny block, default value = 500 BP.",
			false,
			500,
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

		TCLAP::ValueArg<std::string> variantsFile("v",
			"variant",
			"VCF file with found variants.",
			false,
			"variant.vcf",
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
		std::set<size_t> referenceSequenceId;
		std::vector<SyntenyFinder::FASTARecord> chrList;
		std::string fileName[] = {referenceFile.getValue(), assemblyFile.getValue()};		
		for(size_t file = 0; file < 2; file++)
		{
			SyntenyFinder::FASTAReader reader(fileName[file]);
			if(!reader.IsOk())
			{
				throw std::runtime_error(("Cannot open file " + fileName[file]).c_str());
			}
			
			reader.GetSequences(chrList);
			if(file == 0)
			{
				for(size_t i = 0; i < chrList.size(); i++)
				{
					referenceSequenceId.insert(chrList[i].GetId());
				}
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
		std::vector<std::vector<SyntenyFinder::BlockInstance> > history(stage.size() + 1);
		finder.GenerateSyntenyBlocks(stage[0].first, trimK, stage[0].first, history[0], false, PutProgressChr);
		for(size_t i = 0; i < stage.size(); i++)
		{
			trimK = std::min(trimK, stage[i].first);
			std::cout << "Simplification stage " << i + 1 << " of " << stage.size() << std::endl;
			std::cout << "Enumerating vertices of the graph, then performing bulge removal..." << std::endl;
			finder.PerformGraphSimplifications(stage[i].first, stage[i].second, maxIterations.getValue(), PutProgressChr);
			if(i < stage.size() - 1)
			{
				finder.GenerateSyntenyBlocks(stage[i].first, trimK, stage[i].first, history[i + 1], false, PutProgressChr);
			}
		}
				
		std::cout << "Finding variants and generating the output..." << std::endl;
		size_t lastK = std::min(stage.back().first, static_cast<int>(minBlockSize.getValue()));
		trimK = std::min(trimK, static_cast<int>(minBlockSize.getValue()));
		finder.GenerateSyntenyBlocks(lastK, trimK, minBlockSize.getValue(), history.back(), false, PutProgressChr);
		size_t refSeqId = chrList[0].GetId();				
		std::vector<SyntenyFinder::Variant> variant;
		std::vector<SyntenyFinder::Reversal> reversal;
		std::vector<SyntenyFinder::Translocation> translocation;		
		SyntenyFinder::VariantCaller caller(chrList, referenceSequenceId, history, trimK, minBlockSize.getValue());
		caller.CallVariants(variant);
		caller.GetHistory(history);
		std::vector<SyntenyFinder::BlockInstance> blockList = history.back();

		SyntenyFinder::OutputGenerator generator(chrList);
		SyntenyFinder::CreateDirectory(outFileDir.getValue());
		const std::string defaultCoordsFile = outFileDir.getValue() + "/blocks_coords.txt";
		const std::string defaultPermutationsFile = outFileDir.getValue() + "/genomes_permutations.txt";
		const std::string defaultCoverageReportFile = outFileDir.getValue() + "/coverage_report.txt";
		const std::string defaultSequencesFile = outFileDir.getValue() + "/blocks_sequences.fasta";
		const std::string defaultGraphFile = outFileDir.getValue() + "/de_bruijn_graph.dot";
		const std::string defaultCircosDir = outFileDir.getValue() + "/circos";
		const std::string defaultCircosFile = defaultCircosDir + "/circos.conf";
		const std::string defaultD3File = outFileDir.getValue() + "/d3_blocks_diagram.html";		
		const std::string defaultPlainVariantFile = outFileDir.getValue() + "/variant.txt";
		std::string defaultVariantFile = variantsFile.getValue();
		if(outFileDir.isSet())
		{
			defaultVariantFile = outFileDir.getValue() + "/" + defaultVariantFile;
		}
				
		generator.GenerateVariantOutput(variant, referenceSequenceId, fileName[1], defaultVariantFile);
		std::ofstream plainVariantStream(defaultPlainVariantFile.c_str());
		std::copy(variant.begin(), variant.end(), std::ostream_iterator<SyntenyFinder::Variant>(plainVariantStream, "\n"));
		if(outFileDir.isSet())
		{			
			generator.ListChromosomesAsPermutations(blockList, defaultPermutationsFile);
			generator.GenerateReport(blockList, defaultCoverageReportFile);
			generator.ListBlocksIndices(blockList, defaultCoordsFile);
			generator.ListBlocksSequences(blockList, defaultSequencesFile);
			generator.GenerateD3Output(blockList, defaultD3File);		
			generator.GenerateCircosOutput(blockList, defaultCircosFile, defaultCircosDir);
		}
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
	
	return 0;
}
