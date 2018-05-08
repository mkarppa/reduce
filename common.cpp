
/******************************************************* Common subroutines. */

// #include <stdio.h>
// #include <stdlib.h>
// #include <assert.h>



#include "common.hpp"
#include <sys/utsname.h>
#include <cstring>
#include <ctime>
#include <stdexcept>

/********************************************************** Error reporting. */

/******************************************************** Get the host name. */

#define MAX_HOSTNAME 256

const char *common_hostname(void)
{
    static char hn[MAX_HOSTNAME];

    struct utsname undata;
    uname(&undata);
    strcpy(hn, undata.nodename);
    return hn;
}

/****************************************************************** Timings. */

#define TIME_STACK_CAPACITY 256

clock_t time_stack[TIME_STACK_CAPACITY];
int     time_stack_top = -1;
int     do_time        = 0;

void enable_timing(void)
{
    do_time = 1;
}

void disable_timing(void)
{
    do_time = 0;
}

void push_time(void) 
{
    if (do_time) {
        if(time_stack_top + 1 > TIME_STACK_CAPACITY)
          throw std::runtime_error("timing stack out of capacity");
        time_stack[++time_stack_top] = clock();
    }
}

double pop_time(void)
{
    if(do_time) {
        clock_t stop = clock();
        if(time_stack_top < 0)
          throw std::runtime_error("pop on an empty timing stack");
        clock_t start = time_stack[time_stack_top--];
        return (double) (1000.0*((double) (stop-start))/CLOCKS_PER_SEC);    
    } else {
        return -1.0;
    }
}

void pop_print_time(const char *legend)
{
    if(do_time) {
        fprintf(stderr, " {%s: %.2fms}", legend, pop_time());
        fflush(stderr);
    }
}

/***************************************************************** Printing. */

void print_int_array(FILE *out, int l, const int *a)
{   
    for(int cursor = 0; cursor < l; cursor++) {
        int lookahead = cursor + 1;
        for(; 
	    lookahead < l && a[lookahead-1]+1 == a[lookahead] ; 
            lookahead++)
            ;
        if(lookahead - cursor > 5) {
            fprintf(out, 
                    "%s%d %d ... %d",
                    cursor == 0 ? "" : " ",
                    a[cursor] + 1, a[cursor+1] + 1, a[lookahead-1] + 1);
            cursor = lookahead - 1;
        } else {
            fprintf(out,
                    "%s%d",
                    cursor == 0 ? "" : " ",
                    a[cursor] + 1);
        }
    }
}

