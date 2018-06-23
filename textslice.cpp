/*
    Copyright 2005-2012 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

//! Holds a slice of text.

#include "textslice.h"


/** Instances *must* be allocated/freed using methods herein, because the C++ declaration
    represents only the header of a much larger object in memory. */
//! Free a TextSlice object 
void TextSlice::free() {
  tbb::tbb_allocator<char>().deallocate((char*)this,sizeof(TextSlice)+(physical_end-begin())+1);
} 
//! Pointer to beginning of sequence
char* TextSlice::begin() {return (char*)(this+1);}
//! Pointer to one past last character in sequence
char* TextSlice::end() {return logical_end;}
//! Length of sequence
size_t TextSlice::size() const {return logical_end-(char*)(this+1);}
//! Maximum number of characters that can be appended to sequence
size_t TextSlice::avail() const {return physical_end-logical_end;}
//! Append sequence [first,last) to this sequence.
void TextSlice::append( char* first, char* last ) {
  memcpy( logical_end, first, last-first );
  logical_end += last-first;
}
//! Set end() to given value.
void TextSlice::set_end( char* p ) {logical_end=p;}

