/***********************************************************************
 * Authors: Josh Marshall
 *
 * This code is intended to replicate MUMmer 3.x's functionality created
 * by Dr. Stefan Kurtz.  The reason for recoding the project is to
 * critically improve code quality, reduce LOC, and reduce memory usage.
 *
 * This code is distributed under the GPLv3
 *
 *    This file is part of MUMmer 4.0 .
 *
 *    MUMmer 4.0 is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published
 * 		by the Free Software Foundation, either version 3 of the License,
 * 		or (at your option) any later version.
 *
 *    MUMmer 4.0 is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ********************************************************************/

/***********************************************************************
 * DOCUMENTATION:
 *
 * This file contains main() and is the entry point of the program.  It
 * is called and used as per the usage statement.
 *
 * USAGE:
 *
 * mummer [options] <reference file> <query file1> . . . [query file32]
 *
 * There must be exactly one reference file and at least one query file.
 * Both the reference and query files should be in multi-FastA format
 * and may contain any set of upper and lowercase characters, thus DNA
 * and protein sequences are both allowed and matching is case
 * insensitive. The maximum number of query files is 32, but there is no
 * limit on how many sequences each reference or query file may contain.
 *
 * Program options
 * -mum
 * 		Compute MUMs, i.e. matches that are unique in both the reference
 * 		and query
 *
 * -mumreference
 * 		Compute MUM-candidates, i.e. matches that are unique in the
 * 		reference but not necessarily in the query
 *
 * -maxmatch
 * 		Compute all maximal matches regardless of their uniqueness
 *
 * -n
 * 		Only match the characters a, c, g, or t (case insensitive)
 *
 * -l int
 * 		Minimum match length (default 20)
 *
 * -b
 * 		Compute both forward and reverse complement matches
 *
 * -r
 * 		Only compute reverse complement matches
 *
 * -s
 * 		Show the matching substring in the output
 *
 * -c
 * 		Report the query position of a reverse complement match relative
 * 		to the forward strand of the query sequence
 *
 * -F
 * 		Force 4 column output format that prepends every match line with
 * 		the reference sequence identifier
 *
 * -L
 * 		Show the length of the query sequence on the header line
 *
 * -help
 * 		Show the possible options and exit
 *
 * ********************************************************************/

/***********************************************************************
 * INCLUDES*************************************************************
 * ********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cerrno>
#include <libfasta/fasta.h>
#include <suffixarray.h>

#include "strings.h"
#include "flags.h"
#include "defines.h"

using namespace std;

/***********************************************************************
 * STRUCTURES***********************************************************
 * ********************************************************************/


typedef struct globalArgs{
	int minimumMatchLength;
	unsigned long flags;
	suffixArrayContainer *sequences;
	size_t seqCount;
}globalArgs;



/***********************************************************************
 * HELPER FUNCTIONS*****************************************************
 * ********************************************************************/
void addRecord(globalArgs *storage, suffixArrayContainer toAdd){
		storage->seqCount++;
		storage->sequences = (suffixArrayContainer*) realloc(storage->sequences, 
											sizeof(suffixArrayContainer) * storage->seqCount);
		memcpy(&storage->sequences[storage->seqCount - 1], &toAdd, sizeof(suffixArrayContainer));
}


void verifyInput(int argc, char** argv){
	int i, toIgnore;
	unsigned long flags;
	bool dataToProcess;

	dataToProcess = false;
	flags = 0;
	i = 0;
	
	while (++i < argc){
		if(!strcmp(argv[i], "-mum")){										flags |= COMPUTE_MUMS;
			continue;
		}else if(!strcmp(argv[i], "-mumreference")){		flags |= MUM_REFERENCE;
			continue;
		}else if(!strcmp(argv[i], "-maxmatch")){				flags |= MAX_MATCH;
			continue;
		}else if(!strcmp(argv[i], "-n")){								flags |= MATCH_N_CHARS;
			continue;
		}else if(!strcmp(argv[i], "-l")){								flags |= MATCH_MINIMUM_LENGTH;
			if(EOF == sscanf(argv[++i], "%d", &toIgnore))	flags |= HELP_STATEMENT;
			if(toIgnore < 1)															flags |= HELP_STATEMENT;
			continue;
		}else if(!strcmp(argv[i], "-b")){								flags |= MATCH_BIDIRECTIONAL_COMPLEMENT_MATCHES;
			continue;
		}else if(!strcmp(argv[i], "-r")){								flags |= MATCH_REVERSE_COMPLEMENT_MATCHES_ONLY;
			continue;
		}else if(!strcmp(argv[i], "-s")){								flags |= SHOW_MATCHES_IN_OUTPUT;
			continue;
		}else if(!strcmp(argv[i], "-c")){								flags |= REPORT_MATCH_OFFSET;
			continue;
		}else if(!strcmp(argv[i], "-F")){								flags |= FOUR_COUMN_OUTPUT;
			continue;
		}else if(!strcmp(argv[i], "-L")){								flags |= SHOW_QUERY_LENGTH_ON_SEPERATE_HEADER_LINE;
			continue;
		}else if(!strcmp(argv[i], "-help")){						flags |= HELP_STATEMENT;
			continue;
		}else{//check if the file can be found and has a valid format.
			FASTA *file = fasta_open(argv[i], 0, NULL);
			if(file == NULL){
				flags |= HELP_STATEMENT;
				continue;
			}

			FASTA_rec_t *fasta_record = fasta_read(file, NULL, 0, NULL);
			if(fasta_record == NULL){
				fasta_close(file);
				flags |= HELP_STATEMENT;
				continue;
			}

			dataToProcess = true;
			fasta_rec_free(fasta_record);
			fasta_close(file);
		}
	}

	if(flags & HELP_STATEMENT || !dataToProcess){
		fprintf(stderr, USAGE_STATEMENT);
		exit(-errno);
	}
}


globalArgs parseArgs(int argc, char** argv){
	globalArgs toReturn;
	int i;

	toReturn.sequences = NULL;
	toReturn.seqCount = 0;
	toReturn.minimumMatchLength = DEFAULT_MINIMUM_MATCH_LENGTH;
	i = 0;

	while (++i < argc){
		if(!strcmp(argv[i], "-mum")){
			toReturn.flags |= COMPUTE_MUMS;
			continue;
		}else if(!strcmp(argv[i], "-mumreference")){
			toReturn.flags |= MUM_REFERENCE;
			continue;
		}else if(!strcmp(argv[i], "-maxmatch")){
			toReturn.flags |= MAX_MATCH;
			continue;
		}else if(!strcmp(argv[i], "-n")){
			toReturn.flags |= MATCH_N_CHARS;
			continue;
		}else if(!strcmp(argv[i], "-l")){
			toReturn.flags |= MATCH_MINIMUM_LENGTH;
			sscanf(argv[++i], "%d", &toReturn.minimumMatchLength);
			continue;
		}else if(!strcmp(argv[i], "-b")){
			toReturn.flags |= MATCH_BIDIRECTIONAL_COMPLEMENT_MATCHES;
			continue;
		}else if(!strcmp(argv[i], "-r")){
			toReturn.flags |= MATCH_REVERSE_COMPLEMENT_MATCHES_ONLY;
			continue;
		}else if(!strcmp(argv[i], "-s")){
			toReturn.flags |= SHOW_MATCHES_IN_OUTPUT;
			continue;
		}else if(!strcmp(argv[i], "-c")){
			toReturn.flags |= REPORT_MATCH_OFFSET;
			continue;
		}else if(!strcmp(argv[i], "-F")){
			toReturn.flags |= FOUR_COUMN_OUTPUT;
			continue;
		}else if(!strcmp(argv[i], "-L")){
			toReturn.flags |= SHOW_QUERY_LENGTH_ON_SEPERATE_HEADER_LINE;
			continue;
		}else{//Open files
			//This is not *totally* correct, as error checking should be done
			//again, but we're not going to check again.  If there's problems
			//in this segment of code, it's because of extreme instability on
			//the host system or the libfasta library.
			FASTA *inFile = fasta_open(argv[i], 0, NULL);
			FASTA_rec_t *tmpFastaRecord;
			while((tmpFastaRecord = fasta_read(inFile, NULL, 0, NULL)) != NULL){
				suffixArrayContainer tmp = copySequenceToLocal(makeSuffixArray(tmpFastaRecord->seq_mem, tmpFastaRecord->seq_len));
				addRecord(&toReturn, tmp);
				fasta_rec_free(tmpFastaRecord);
			}
			fasta_close(inFile);
		}
	}
	return toReturn;
}


void freeGlobalArgs(globalArgs toFree){
	size_t endIndex = toFree.seqCount;
	for(size_t i = 0; i < endIndex; i++){
		freeSuffixArray(&(toFree.sequences[i]));
	}
}


int main(int argc, char** argv){
	//DECLARATIONS////////////////////////////////////////////////////////
	globalArgs processedArgs;

	//INITIALIZATION AND INPUT PARSING////////////////////////////////////

	//If verifyInput() fails, it handles forced program termination and
	//error feedback to the user on stderr.
	verifyInput(argc, argv);
	processedArgs = parseArgs(argc, argv);

	//MAIN LOOP///////////////////////////////////////////////////////////

	for(size_t i = 0; i < processedArgs.seqCount; i++){
		for(size_t j = 0; j < processedArgs.sequences[i].length; j++){
			printf("%c", processedArgs.sequences[i].sequence[j]);
		}
		printf("\n\n");
	}

	//CLEANUP/////////////////////////////////////////////////////////////
	freeGlobalArgs(processedArgs);

	return 0;
}
