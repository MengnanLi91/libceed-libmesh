#pragma once
#include <libmesh/fe_type.h>

#include "structs.h"

void print_FEproblemData(const FEproblemData_ &data);

void printCeedVector(CeedVector vec);

// General function to print the contents of any struct using Boost.PFR
void verifyQFunctionContext(CeedQFunctionContext build_ctx);
