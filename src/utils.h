#include "Rcpp.h"
#ifndef UTILS_H
#define UTILS_H

bool check_logical_scalar(Rcpp::RObject x, const char* thing);

int check_integer_scalar(Rcpp::RObject x, const char* thing);

double check_numeric_scalar(Rcpp::RObject x, const char* thing);

Rcpp::String check_string(Rcpp::RObject x, const char* thing);

#endif
