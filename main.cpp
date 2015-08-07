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
 *     by the Free Software Foundation, either version 3 of the License,
 *     or (at your option) any later version.
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
 *     Compute MUMs, i.e. matches that are unique in both the reference
 *     and query
 *
 * -mumreference
 *     Compute MUM-candidates, i.e. matches that are unique in the
 *     reference but not necessarily in the query
 *
 * -maxmatch
 *     Compute all maximal matches regardless of their uniqueness
 *
 * -n
 *     Only match the characters a, c, g, or t (case insensitive)
 *
 * -l int
 *     Minimum match length (default 20)
 *
 * -b
 *     Compute both forward and reverse complement matches
 *
 * -r
 *     Only compute reverse complement matches
 *
 * -s
 *     Show the matching substring in the output
 *
 * -c
 *     Report the query position of a reverse complement match relative
 *     to the forward strand of the query sequence
 *
 * -F
 *     Force 4 column output format that prepends every match line with
 *     the reference sequence identifier
 *
 * -L
 *     Show the length of the query sequence on the header line
 *
 * -help
 *     Show the possible options and exit
 *
 * ********************************************************************/

/***********************************************************************
 * INCLUDES*************************************************************
 * ********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <libfasta/fasta.h>
#include <suffixarray.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

#include "strings.h"
#include "flags.h"
#include "defines.h"

//using namespace std;

/***********************************************************************
 * STRUCTURES***********************************************************
 * ********************************************************************/


typedef struct globalArgs{
  int minimumMatchLength;
  unsigned long flags;
  suffixArray *sequences;
  size_t seqCount;
}globalArgs;



/***********************************************************************
 * HELPER FUNCTIONS*****************************************************
 * ********************************************************************/
void addRecord(globalArgs *storage, suffixArray toAdd){
    storage->seqCount++;
    storage->sequences = (suffixArray*) realloc(storage->sequences, 
                      sizeof(suffixArray) * storage->seqCount);
    memcpy(&storage->sequences[storage->seqCount - 1], &toAdd, sizeof(suffixArray));
}


void *loadFASTAAsESA(void *fileNameIn){
  char* fileName = (char*) fileNameIn;
  FASTA *inFile;
  FASTA_rec_t *tmpFastaRecord;
  suffixArrayCaster **toReturn;
  int ESAindex;
  
#ifdef DEBUG
  time_t timeStat;
  time(&timeStat);
	
	char* str0 = "Starting file open\n";
	write(fileno(stdout), str0, strlen(str0));
#endif

  ESAindex = 0;
  fileName = (char*) fileNameIn;
  inFile = fasta_open(fileName, 0, NULL);
  toReturn = (suffixArrayCaster**) malloc(sizeof(suffixArray) * (fasta_count(inFile) + 1));

#ifdef DEBUG
  char *str1 = "File load took %.f seconds\n";
	size_t strLength = 1024;
	char *toPrint = (char*) calloc(sizeof(char), strLength);
	time_t currTime;
	time(&currTime);
	sprintf(toPrint, str1, difftime(timeStat, currTime));
	write(fileno(stdout), toPrint, strlen(toPrint));
	
	time(&timeStat);
#endif
  
  while((tmpFastaRecord = fasta_read(inFile, NULL, 0, NULL)) != NULL){
    toReturn[ESAindex] = (suffixArrayCaster*) malloc(sizeof(suffixArray)); 

#ifdef DEBUG
		char *str2 = "FASTA read took %.f seconds\n";
		memset(toPrint, 0, strLength);
		time(&currTime);
		sprintf(toPrint, str2, difftime(timeStat, currTime));
		write(fileno(stdout), toPrint, strlen(toPrint));
	
		time(&timeStat);
#endif

    suffixArray tmp = copySequenceToLocal(makeSuffixArray(tmpFastaRecord->seq_mem, tmpFastaRecord->seq_len));
    
    memcpy(toReturn[ESAindex++], &tmp, sizeof(suffixArray));

#ifdef DEBUG
		char *str3 = "suffixArray construction took %.f seconds\n";
		memset(toPrint, 0, strLength);
		time(&currTime);
		sprintf(toPrint, str3, difftime(timeStat, currTime));
		write(fileno(stdout), toPrint, strlen(toPrint));
	
		time(&timeStat);
#endif
    
    fasta_rec_free(tmpFastaRecord);
		
#ifdef DEBUG
		char *str4 = "FASTA record free took %.f seconds\n";
		memset(toPrint, 0, strLength);
		time(&currTime);
		sprintf(toPrint, str4, difftime(timeStat, currTime));
		write(fileno(stdout), toPrint, strlen(toPrint));
	
		time(&timeStat);
#endif
  }
  fasta_close(inFile);
		
#ifdef DEBUG
		char *str5 = "File close took %.f seconds\n";
		memset(toPrint, 0, strLength);
		time(&currTime);
		sprintf(toPrint, str5, difftime(timeStat, currTime));
		write(fileno(stdout), toPrint, strlen(toPrint));
	
		time(&timeStat);
#endif
  
  toReturn[ESAindex] = NULL;
  return (void*) toReturn;
}


globalArgs parseArgs(int argc, char** argv){
  globalArgs toReturn;
  int i, *fileArgs, numFileArgs;
  pthread_t *loaders;
  

  toReturn.sequences = NULL;
  toReturn.seqCount = 0;
  toReturn.minimumMatchLength = DEFAULT_MINIMUM_MATCH_LENGTH;
  i = 0;
  fileArgs = NULL;
  numFileArgs = 0;
  loaders = NULL;

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
      if(EOF == sscanf(argv[++i], "%d", &toReturn.minimumMatchLength)){
        fprintf(stderr, USAGE_STATEMENT);
        exit(-errno);
      }
      if(toReturn.minimumMatchLength < 1){
        fprintf(stderr, USAGE_STATEMENT);
        exit(-errno);
      }
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
    }else if(!strcmp(argv[i], "-help")){
      fprintf(stderr, USAGE_STATEMENT);
      exit(-errno);
    }else{//Open files
      //This is not *totally* correct, as error checking should be done
      //again, but we're not going to check again.  If there's problems
      //in this segment of code, it's because of extreme instability on
      //the host system or the libfasta library.
      
      numFileArgs++;
      fileArgs = (int*) realloc(fileArgs, sizeof(int) * numFileArgs);
      fileArgs[numFileArgs -1] = i;
    }
  }
  
  if(numFileArgs == 0){
    fprintf(stderr, USAGE_STATEMENT);
    exit(-errno);
  }
    
  loaders = (pthread_t*) malloc(sizeof(pthread_t) * numFileArgs);
  for(i = 0; i < numFileArgs; i++)
    pthread_create(&loaders[i], NULL, loadFASTAAsESA, argv[fileArgs[i]]);
  
  for(i = 0; i < numFileArgs; i++){
    suffixArray **returned;
    int index = 0;
    pthread_join(loaders[i], (void**) &returned);
    while(returned[index] != NULL){
      toReturn.sequences = (suffixArray*) realloc(toReturn.sequences, sizeof(suffixArray) * toReturn.seqCount + 1);
      memcpy(&toReturn.sequences[toReturn.seqCount++], returned[index++], sizeof(suffixArray));
      free(returned[index - 1]);  
    }
    free(returned);
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
