
/******************************************************* Common subroutines. */

#ifndef REDUCE2_COMMON_HPP
#define REDUCE2_COMMON_HPP

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <iostream>
#include <vector>

  

const char *  common_hostname         (void);

void          enable_timing           (void);
void          disable_timing          (void);
void          push_time               (void);
double        pop_time                (void);
void          pop_print_time          (const char *legend);
void print_int_array(FILE *out, int l, const int *a);

namespace reduce {
  template<typename T>
  void printArray(std::ostream& out, const std::vector<T>& array) {
    for (size_t cursor = 0; cursor < array.size(); ++cursor) {
      size_t lookahead = cursor + 1;
      while (lookahead < array.size() &&
             array[lookahead-1]+1 == array[lookahead])
        ++lookahead;

      if (lookahead - cursor > 5) {
        out << (cursor == 0 ? "" : " ")
            << array[cursor] + 1 << " "
            << array[cursor+1] + 1 << " .. "
            << array[lookahead-1] + 1;
        cursor = lookahead - 1;
      }
      else {
        out << (cursor == 0 ? "" : " ")
            <<  array[cursor] + 1;
      }
    }
  }
}

#endif // REDUCE2_COMMON_HPP
