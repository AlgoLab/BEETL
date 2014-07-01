/**
 ** Copyright (c) 2011 Illumina, Inc.
 **
 **
 ** This software is covered by the "Illumina Non-Commercial Use Software
 ** and Source Code License Agreement" and any user of this software or
 ** source file is bound by the terms therein (see accompanying file
 ** Illumina_Non-Commercial_Use_Software_and_Source_Code_License_Agreement.pdf)
 **
 ** This file is part of the BEETL software package.
 **
 ** Citation: Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 ** Lightweight BWT Construction for Very Large String Collections.
 ** Proceedings of CPM 2011, pp.219-231
 **
 **/

#include "EndPosFile.hh"

#include <cassert>
#include <fstream>
#include <stdexcept>

#include "Types.hh"

using namespace std;


EndPosFile::EndPosFile( const string &bwtFilenamePrefix )
    : file_( bwtFilenamePrefix + "-end-pos" )
    , sequenceGroupCount_( 0 )
    , sequenceCountInGroup_( 0 )
    , hasRevComp_( 0 )
{
    file_.read( reinterpret_cast< char * >( &sequenceGroupCount_ ), sizeof( SequenceNumber ) );
    file_.read( reinterpret_cast< char * >( &sequenceCountInGroup_ ), sizeof( uint8_t ) );
    file_.read( reinterpret_cast< char * >( &hasRevComp_ ), sizeof( uint8_t ) );
    //    assert( file_.good() );
    dollarSignCount_ = sequenceGroupCount_ * sequenceCountInGroup_ * ( hasRevComp_ ? 2 : 1 );

#ifdef SLIM_STRUCTURES
	 if (sequenceCountInGroup_ > 1 || hasRevComp_) {
		throw std::logic_error("'end-pos' files with paired-end reads or reverse and complement reads are not supported!");
	 }
#endif

}

SequenceNumber EndPosFile::convertDollarNumToSequenceNum( const SequenceNumber dollarNum )
//SequenceNumber EndPosFile_convertDollarNumToSequenceNum( const SequenceNumber dollarNum )
{
    assert( file_.good() && "Error: -end-pos file not readable" );

    assert( dollarNum < dollarSignCount_ );
    /*
      if ( dollarPos >= numDollarEntries )
      {
      cout << "Warning: dollarPos " << dollarPos << " >= numDollarEntries " << numDollarEntries << endl;
      //                    continue;
      dollarPos %= numDollarEntries;
      }
    */
#ifdef SLIM_STRUCTURES
    file_.seekg( sizeof( SequenceNumber ) + 2 * sizeof( uint8_t ) + ( dollarNum ) * sizeof( SequenceNumber ) );
#else
    file_.seekg( sizeof( SequenceNumber ) + 2 * sizeof( uint8_t ) + ( dollarNum ) * ( sizeof( SequenceNumber ) + sizeof( uint8_t ) ) );
#endif

	 SequenceNumber sequenceGroupNum;
    file_.read( reinterpret_cast< char * >( &sequenceGroupNum ), sizeof( SequenceNumber ) );
    uint8_t positionInGroup= (uint8_t)1;
#ifndef SLIM_STRUCTURES
    file_.read( reinterpret_cast< char * >( &positionInGroup ), sizeof( uint8_t ) );
#endif

    SequenceNumber sequenceNum = sequenceGroupNum + positionInGroup * sequenceGroupCount_;

    return sequenceNum;
}
