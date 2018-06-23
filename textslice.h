#ifndef __TEXTSLICE_H_INCLUDED__
#define __TEXTSLICE_H_INCLUDED__ 

#include <cstring>
#include "tbb/tbb.h"

class TextSlice {
    //! Pointer to one past last character in sequence
    char* logical_end;
    //! Pointer to one past last available byte in sequence.
    char* physical_end;
public:
    //! Allocate a TextSlice object that can hold up to max_size characters.
    static TextSlice* allocate( size_t max_size ) {
      // +1 leaves room for a terminating null character.
      TextSlice* t = (TextSlice*)tbb::tbb_allocator<char>().allocate( sizeof(TextSlice)+max_size+1 );
      t->logical_end = t->begin();
      t->physical_end = t->begin()+max_size;
      return t;
    }

    //! Free a TextSlice object 
    void free();
    //! Pointer to beginning of sequence
    char* begin();
    //! Pointer to one past last character in sequence
    char* end();
    //! Length of sequence
    size_t size() const;
    //! Maximum number of characters that can be appended to sequence
    size_t avail() const;
    //! Append sequence [first,last) to this sequence.
    void append( char* first, char* last );
    //! Set end() to given value.
    void set_end( char* p );
};

#endif


