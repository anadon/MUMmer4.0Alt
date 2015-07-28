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
 * This file contains the flags definitions used in MUMmer 4.0 .
 * 
 * ********************************************************************/

#ifndef FLAGS_H
#define FLAGS_H

#define COMPUTE_MUMS 															(((unsigned long)1) << 0x0)
#define MUM_REFERENCE 														(((unsigned long)1) << 0x1)
#define MAX_MATCH																	(((unsigned long)1) << 0x2)
#define MATCH_N_CHARS															(((unsigned long)1) << 0x3)
#define MATCH_MINIMUM_LENGTH 											(((unsigned long)1) << 0x4)
#define MATCH_BIDIRECTIONAL_COMPLEMENT_MATCHES		(((unsigned long)1) << 0x5)
#define MATCH_REVERSE_COMPLEMENT_MATCHES_ONLY			(((unsigned long)1) << 0x6)
#define SHOW_MATCHES_IN_OUTPUT										(((unsigned long)1) << 0x7)
#define REPORT_MATCH_OFFSET												(((unsigned long)1) << 0x8)
#define FOUR_COUMN_OUTPUT													(((unsigned long)1) << 0x9)
#define SHOW_QUERY_LENGTH_ON_SEPERATE_HEADER_LINE	(((unsigned long)1) << 0xa)
#define HELP_STATEMENT														(((unsigned long)1) << 0xb)

#endif