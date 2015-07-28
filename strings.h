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
 * This file contains the strings used in MUMmer 4.0 .
 * 
 * ********************************************************************/

#ifndef STRINGS_H
#define STRINGS_H

#define USAGE_STATEMENT 																								\
"USAGE:																																	\n\
																																				\n\
mummer [options] <reference file> <query file1> . . . [query file32]		\n\
																																				\n\
There must be exactly one reference file and at least one query file.		\n\
Both the reference and query files should be in multi-FastA format 			\n\
and may contain any set of upper and lowercase characters, thus DNA 		\n\
and protein sequences are both allowed and matching is case 						\n\
insensitive. 																														\n\
																																				\n\
Program options																													\n\
-mum 																																		\n\
		Compute MUMs, i.e. matches that are unique in both the reference 		\n\
		and query																														\n\
																																				\n\
-mumreference 																													\n\
		Compute MUM-candidates, i.e. matches that are unique in the 				\n\
		reference but not necessarily in the query													\n\
																																				\n\
-maxmatch 																															\n\
		Compute all maximal matches regardless of their uniqueness					\n\
																																				\n\
-n 																																			\n\
		Only match the characters a, c, g, or t (case insensitive)					\n\
																																				\n\
-l int 																																	\n\
		Minimum match length (default 20)																		\n\
																																				\n\
-b 																																			\n\
		Compute both forward and reverse complement matches									\n\
																																				\n\
-r 																																			\n\
		Only compute reverse complement matches															\n\
																																				\n\
-s 																																			\n\
		Show the matching substring in the output														\n\
																																				\n\
-c 																																			\n\
		Report the query position of a reverse complement match relative 		\n\
		to the forward strand of the query sequence													\n\
																																				\n\
-F 																																			\n\
		Force 4 column output format that prepends every match line with		\n\
		the reference sequence identifier																		\n\
																																				\n\
-L 																																			\n\
		Show the length of the query sequence on the header line						\n\
																																				\n\
-help 																																	\n\
		Show the possible options and exit\n\n"

#endif