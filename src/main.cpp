//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include <tclap/CmdLine.h>
#include "outputgenerator.h"
#include "unrolledlisttest.h"


std::string IntToStr(size_t value)
{
	std::stringstream ss;
	ss << value;
	return ss.str();
}

std::vector<std::pair<int, int> > ReadStageFile(const std::string & fileName)
{
	int count = 0;
	std::ifstream stageFile(fileName.c_str());
	if(!stageFile)
	{
		throw std::runtime_error("cannot open stage file");
	}

	if(!(stageFile >> count))
	{
		throw std::runtime_error("cannot read stage file");
	}

	if(count < 0)
	{
		throw std::runtime_error("number of stages must be nonnegative");
	}

	std::vector<std::pair<int, int> > ret(count);
	for(int i = 0; i < count; i++)
	{
		if(!(stageFile >> ret[i].first >> ret[i].second))
		{
			throw std::runtime_error("too few records in the stage file");
		}

		if(ret[i].first < 2)
		{
			throw std::runtime_error("vertex size in stage record must be at least 2");
		}

		if(ret[i].second < 0)
		{
			throw std::runtime_error("minimum branch size in stage record must be nonnegative");
		}
	}

	return ret;
}

std::vector<std::pair<int, int> > LooseStageFile()
{
	std::pair<int, int> stage[] = 
	{
		std::make_pair(30, 150),
		std::make_pair(100, 1000),
		std::make_pair(1000, 5000),		
		std::make_pair(5000, 15000)
	};

	return std::vector<std::pair<int, int> >(stage, stage + sizeof(stage) / sizeof(stage[0]));
}

std::vector<std::pair<int, int> > FineStageFile()
{
	std::pair<int, int> stage[] = 
	{
		std::make_pair(30, 150),
		std::make_pair(100, 1000),
		std::make_pair(1000, 2500),		
	};

	return std::vector<std::pair<int, int> >(stage, stage + sizeof(stage) / sizeof(stage[0]));
}

void PutProgressChr(size_t progress, SyntenyFinder::BlockFinder::State state)
{
	static size_t prev = 0;
	while(prev < progress)
	{
		prev++;
		std::cout << '.';
		std::cout.flush();
	}

	switch(state)
	{
	case SyntenyFinder::BlockFinder::start:
		prev = 0;
		std::cout << '[';
		break;
	case SyntenyFinder::BlockFinder::end:
		std::cout << "]" << std::endl;
		break;
	case SyntenyFinder::BlockFinder::run:
		;
	}
}

const std::string DELIMITER(80, '-');


int main(int argc, char * argv[])
{	
	testUnrolledList();
	return 0;

	std::stringstream parsets;		
	const std::string parameterSetNameArray[] = {"loose", "fine"};
	std::vector<std::string> parameterSetName(parameterSetNameArray, parameterSetNameArray + sizeof(parameterSetNameArray) / sizeof(parameterSetNameArray[0]));
	std::map<std::string, std::vector<std::pair<int, int> > > defaultParameters;
	defaultParameters["loose"] = LooseStageFile();
	defaultParameters["fine"] = FineStageFile();

	try
	{  
		TCLAP::CmdLine cmd("Program for finding syteny blocks in closely related genomes", ' ', "0.7071");
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

		TCLAP::ValueArg<std::string> graphFile("g",
			"graphfile",
			"File for resulting condensed de Bruijn graph (in dot format), default = not set.",
			false,
			"",
			"file name",
			cmd);

		TCLAP::ValueArg<std::string> sequencesFile("q",
			"sequncesfile",
			"File for sequences of synteny blocks (FASTA format), default = not set.",
			false,
			"",
			"file name",
			cmd);

		TCLAP::ValueArg<std::string> reportFile("r",
			"reportfile",
			"File for coverage report, default = \"coverage_report.txt\".",
			false,
			"coverage_report.txt",
			"file name",
			cmd);

		TCLAP::ValueArg<std::string> coordsFile("c",
			"coordsfile",
			"File for listing coordinates of synteny blocks in simple human-readable format, default = \"block_coords.txt\".",
			false,
			"block_coords.txt",
			"file name",
			cmd);

		TCLAP::ValueArg<std::string> chrFile("p",
			"permfile",
			"File for listing genomes represented as signed permutations of synteny blocks, default = \"genomes_permutations.txt\".",
			false,
			"genomes_permutations.txt",
			"file name",
			cmd);
		
		TCLAP::ValueArg<std::string> circosDir("d",
			"circosdir",
			"Directory for circos files, default = not set",
			false,
			"",
			"dir name",
			cmd);

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

		TCLAP::SwitchArg sharedOnly("a",
			"sharedonly",
			"Output only blocks that occur exactly once in each input sequence.",			
			cmd,
			false);

		TCLAP::UnlabeledMultiArg<std::string> fileName("filenames",
			"FASTA file(s) with nucleotide sequences",
			true,
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

		std::vector<SyntenyFinder::FASTARecord> chrList;
		for(std::vector<std::string>::const_iterator it = fileName.begin(); it != fileName.end(); it++)
		{
			SyntenyFinder::FASTAReader reader(*it);
			if(!reader.IsOk())
			{
				throw std::runtime_error(("Cannot open file " + *it).c_str());
			}

			reader.GetSequences(chrList);
		}

		SyntenyFinder::BlockFinder finder(chrList);
		for(size_t i = 0; i < stage.size(); i++)
		{
			std::cout << "Simplification stage " << i + 1 << " of " << stage.size() << std::endl;
			std::cout << "Enumerating vertices of the graph, then performing bulge removal..." << std::endl;
			finder.PerformGraphSimplifications(stage[i].first, stage[i].second, maxIterations.getValue(), PutProgressChr);
		}
		
		std::vector<SyntenyFinder::BlockInstance> blockList;
		std::cout << "Finding synteny blocks and generating the output..." << std::endl;
		size_t lastK = std::min(stage.back().first, static_cast<int>(minBlockSize.getValue()));
		finder.GenerateSyntenyBlocks(lastK, minBlockSize.getValue(), blockList, sharedOnly.getValue(), PutProgressChr);
		SyntenyFinder::OutputGenerator generator(chrList, blockList);

	//	const std::string templateCircosConf = "circos.template.conf";
		const std::string defaultCircosOutFile = "circos.conf";
		
		bool doOutput[] = {true, true, true, sequencesFile.isSet(), circosDir.isSet()};
		boost::function<void()> outFunction[] = 
		{
			boost::bind(&SyntenyFinder::OutputGenerator::ListChromosomesAsPermutations, boost::cref(generator), chrFile.getValue()),
			boost::bind(&SyntenyFinder::OutputGenerator::GenerateReport, boost::cref(generator), reportFile.getValue()),
			boost::bind(&SyntenyFinder::OutputGenerator::ListBlocksIndices, boost::cref(generator), coordsFile.getValue()),
			boost::bind(&SyntenyFinder::OutputGenerator::ListBlocksSequences, boost::cref(generator), sequencesFile.getValue()),		
			boost::bind(&SyntenyFinder::OutputGenerator::GenerateCircosOutput, boost::cref(generator), circosDir.getValue() + "/" + defaultCircosOutFile, circosDir.getValue())
		};

		size_t length = sizeof(doOutput) / sizeof(doOutput[0]);
		for(size_t i = 0; i < length; i++)
		{
			if(doOutput[i])
			{
				outFunction[i]();
			}
		}

		if(graphFile.isSet())
		{
			std::stringstream buffer;
			finder.SerializeCondensedGraph(lastK, buffer, PutProgressChr);
			generator.OutputBuffer(graphFile.getValue(), buffer.str());
		}

		std::cout.setf(std::cout.fixed);
		std::cout.precision(2);
		std::cout << "Time elapsed: " << double(clock()) / CLOCKS_PER_SEC << " seconds" << std::endl;
	} 
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
	}
	catch (std::runtime_error & e)
	{
		std::cerr << "error: " << e.what() << std::endl;
	}

	return 0;
}
