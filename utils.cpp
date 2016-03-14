#include "utils.h"

namespace fctrl {

int whitespace_counter = 0;

std::string align()
{ return std::string(whitespace_counter, ' '); }

void incr_align()
{ whitespace_counter++; }

void decr_align()
{ whitespace_counter--; }

double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

} //namespace fctrl
