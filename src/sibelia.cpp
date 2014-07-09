//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


#include <bitset>
#include <tclap/CmdLine.h>
#include "postprocessor.h"
#include "util.h"

const std::string VERSION("3.0.5");

class GreaterIntegerConstraint: public TCLAP::Constraint<int>
{
public:
	GreaterIntegerConstraint(int value): threshold_(value)
	{
		std::stringstream ss;
		ss << "integer > " << value;
		typeDesc_ = ss.str();
	}

	bool check(const int & value) const
	{
		return value > threshold_;
	}

	std::string shortID() const
	{
		return typeDesc_;
	}

	std::string description() const
	{
		return typeDesc_;   
	}

private:
	int threshold_;
	std::string typeDesc_;
};

size_t GenerateMask(size_t k, size_t weight)
{
	std::vector<size_t> a(k, 1);
	while(weight < k)
	{
		size_t pos = (rand() % (k - 2)) + 1;
		if(a[pos] == 1)
		{
			a[pos] = 0;
			++weight;
		}
	}

	size_t mask = 0;
	for(size_t i = 0; i < k; i++)
	{
		if(a[i] == 1)
		{
			mask = mask | (3 << (i * 2));
		}
	}
	
	return mask;
}

int main(int argc, char * argv[])
{
	srand(time(0));
	signal(SIGINT, SignalHandler);
	signal(SIGABRT, SignalHandler);	
	signal(SIGTERM, SignalHandler);

	std::stringstream parsets;		
	const std::string parameterSetNameArray[] = {"loose", "fine"};
	std::vector<std::string> parameterSetName(parameterSetNameArray, parameterSetNameArray + sizeof(parameterSetNameArray) / sizeof(parameterSetNameArray[0]));
	std::map<std::string, std::vector<ParameterSet> > defaultParameters;
	defaultParameters["loose"] = LooseStageFile();
	defaultParameters["fine"] = FineStageFile();
	GreaterIntegerConstraint greaterThanOne(1);
	GreaterIntegerConstraint greaterThanZero(0);
	try
	{  
		TCLAP::CmdLine cmd("Program for finding syteny blocks in closely related genomes", ' ', VERSION);
		TCLAP::ValueArg<int> maxIterations("i",
			"maxiterations",
			"Maximum number of iterations during a stage of simplification, default = 4.",
			false,
			4,
			&greaterThanZero,
			cmd);

		TCLAP::ValueArg<int> correctBoundariesFlag("",
			"correctboundaries",
			"Correction boundaries of unique synteny blocks, number of iterations.",
			false,
			0,
			"int",
			cmd);

		TCLAP::SwitchArg ("",
			"correctboundaries",
			"",
			cmd,
			false);

		TCLAP::SwitchArg noPostProcessingFlag("",
			"nopostprocess",
			"Do not perform postprocessing (stripe gluing).",
			cmd,
			false);

		TCLAP::SwitchArg GFFFormatFlag("",
			"gff",
			"Use GFF format for reporting blocks coordinates",
			cmd,
			false);		

		TCLAP::SwitchArg allStagesFlag("",
			"allstages",
			"Output coordinates of synteny blocks from all stages",
			cmd,
			false);

		TCLAP::ValueArg<int> lastKValue("",
			"lastk",
			"Value of K used for the synteny blocks inferring.",
			false,
			5000,
			&greaterThanOne,
			cmd);

		TCLAP::ValueArg<std::string> tempFileDir("t",
			"tempdir",
			"Directory where temporary files are stored.",
			false,
			".",
			"dir name",
			cmd);

		TCLAP::ValueArg<std::string> stageFile("k",
			"stagefile",
			"File that contains manually chosen simplifications parameters. See USAGE file for more information.",
			false,
			"",
			"file name");

		TCLAP::SwitchArg hierarchyPicture("v",
			"visualize",
			"Draw circos diagram with blocks at different stages.",
			cmd,
			false);

		TCLAP::SwitchArg graphFile("g",
			"graphfile",
			"Output resulting condensed de Bruijn graph (in dot format).",
			cmd,
			false);

		TCLAP::SwitchArg sequencesFile("q",
			"sequencesfile",
			"Output sequences of synteny blocks (FASTA format).",
			cmd,
			false);		

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

		TCLAP::SwitchArg inRAM("r",
			"inram",
			"Perform all computations in RAM, don't create temp files.",
			cmd,
			false);

		TCLAP::UnlabeledMultiArg<std::string> fileName("filenames",
			"FASTA file(s) with nucleotide sequences.",
			true,
			"fasta files with genomes",
			cmd);

		TCLAP::ValueArg<std::string> outFileDir("o",
			"outdir",
			"Directory where output files are written",
			false,
			".",
			"dir name",
			cmd);		

		cmd.xorAdd(parameters, stageFile);
		cmd.parse(argc, argv);
		std::vector<ParameterSet> stage;
		if(parameters.isSet())
		{
			stage = defaultParameters[parameters.getValue()];
		}
		else
		{
			stage = ReadStageFile(stageFile.getValue());
		}
		
		int trimK = INT_MAX;
		size_t totalSize = 0;
		std::set<size_t> referenceChrId;
		bool allStages = allStagesFlag.isSet();		
		bool hierarchy = hierarchyPicture.isSet();
		bool noPostProcessing = noPostProcessingFlag.isSet();
		int correctBoundaries = correctBoundariesFlag.getValue();

		std::vector<SyntenyFinder::FASTARecord> chrList;
		for(std::vector<std::string>::const_iterator it = fileName.begin(); it != fileName.end(); it++)
		{
			SyntenyFinder::FASTAReader reader(*it);
			if(!reader.IsOk())
			{
				throw std::runtime_error(("Cannot open file " + *it).c_str());
			}

			reader.GetSequences(chrList);
			if(it == fileName.begin())
			{
				for(size_t i = 0; i < chrList.size(); i++)
				{
					referenceChrId.insert(chrList[i].GetId());
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
		
		std::vector<std::vector<SyntenyFinder::BlockInstance> > history(stage.size() + 1);
		std::string tempDir = tempFileDir.isSet() ? tempFileDir.getValue() : outFileDir.getValue();		
		std::auto_ptr<SyntenyFinder::BlockFinder> finder(inRAM.isSet() ? new SyntenyFinder::BlockFinder(chrList) : new SyntenyFinder::BlockFinder(chrList, tempDir));
		SyntenyFinder::Postprocessor processor(chrList, minBlockSize.getValue());

		size_t totalBulges = 0;
		for(size_t i = 0; i < stage.size(); i++)
		{
			trimK = std::min(trimK, stage[i].k);
			if(hierarchy || allStages)
			{
				finder->GenerateSyntenyBlocks(stage[i].k, trimK, stage[i].k, history[i], sharedOnly.getValue());
				if(!noPostProcessing)
				{
					processor.GlueStripes(history[i]);
				}
			}

			size_t bulges;
			std::cout << "Simplification stage " << i + 1 << " of " << stage.size() << std::endl;
			std::cout << "Enumerating vertices of the graph, then performing bulge removal..." << std::endl;			
			bulges = finder->PerformGraphSimplifications(stage[i].k, stage[i].maxBranchSize, maxIterations.getValue(), PutProgressChr);

			std::cout << "Bulges = " << bulges << std::endl;
			std::cerr << std::endl;
		}

		std::cout << "Finding synteny blocks and generating the output..." << std::endl;
		trimK = std::min(trimK, static_cast<int>(minBlockSize.getValue()));
		size_t lastK = lastKValue.isSet() ? lastKValue.getValue() : std::min(stage.back().k, static_cast<int>(minBlockSize.getValue()));
		finder->GenerateSyntenyBlocks(lastK, trimK, minBlockSize.getValue(), history.back(), sharedOnly.getValue(), PutProgressChr);
		if(!noPostProcessing)
		{
			processor.GlueStripes(history.back());
		}

		for(size_t i = 0; i < correctBoundaries && processor.ImproveBlockBoundaries(history.back()); i++);

		bool oldFormat = !GFFFormatFlag.isSet();
		SyntenyFinder::OutputGenerator generator(chrList);
		SyntenyFinder::CreateOutDirectory(outFileDir.getValue());
		boost::function<void(const std::vector<SyntenyFinder::BlockInstance>&, const std::string&)> coordsWriter = 
			oldFormat ? boost::bind(&SyntenyFinder::OutputGenerator::ListBlocksIndices, boost::cref(generator), _1, _2)
					  : boost::bind(&SyntenyFinder::OutputGenerator::ListBlocksIndicesGFF, boost::cref(generator), _1, _2);
		const std::string defaultCoordsFile = outFileDir.getValue() + "/blocks_coords" + (oldFormat ? ".txt" : ".gff");
		const std::string defaultPermutationsFile = outFileDir.getValue() + "/genomes_permutations.txt";
		const std::string defaultCoverageReportFile = outFileDir.getValue() + "/coverage_report.txt";
		const std::string defaultSequencesFile = outFileDir.getValue() + "/blocks_sequences.fasta";
		const std::string defaultGraphFile = outFileDir.getValue() + "/de_bruijn_graph.dot";
		const std::string defaultCircosDir = outFileDir.getValue() + "/circos";
		const std::string defaultCircosFile = defaultCircosDir + "/circos.conf";
		const std::string defaultD3File = outFileDir.getValue() + "/d3_blocks_diagram.html";      		
		if(allStages)
		{			
			for(size_t i = 0; i < history.size(); i++)
			{
				std::stringstream file;
				file << outFileDir.getValue() << "/blocks_coords" << i << (oldFormat ? ".txt" : ".gff");
				coordsWriter(history[i], file.str());
			}
		}
		else
		{
			coordsWriter(history.back(), defaultCoordsFile);
		}

		generator.ListChromosomesAsPermutations(history.back(), defaultPermutationsFile);
		generator.GenerateReport(history.back(), defaultCoverageReportFile);		
		generator.GenerateD3Output(history.back(), defaultD3File);
		if(sequencesFile.isSet())
		{
			generator.ListBlocksSequences(history.back(), defaultSequencesFile);
		}

		if(!hierarchy)
		{
			generator.GenerateCircosOutput(history.back(), defaultCircosFile, defaultCircosDir);
		}
		else
		{
			generator.GenerateHierarchyCircosOutput(history, defaultCircosFile, defaultCircosDir);
		}		

		if(graphFile.isSet())
		{
			std::stringstream buffer;
			finder->SerializeCondensedGraph(lastK, buffer, PutProgressChr);
			generator.OutputBuffer(defaultGraphFile, buffer.str());
		}

		std::cout.setf(std::cout.fixed);
		std::cout.precision(2);
	//	std::cout << "Time elapsed: " << double(clock()) / CLOCKS_PER_SEC << " seconds" << std::endl;
	} 
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	catch (std::runtime_error & e)
	{
		std::cerr << "error: " << e.what() << std::endl;
		return 1;
	}
	catch(...)
	{
		SyntenyFinder::TempFile::Cleanup();
		return 1;
	}

	return 0;
}
