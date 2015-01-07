#ifndef SIMIT_EXPR_QUERIES_H
#define SIMIT_EXPR_QUERIES_H

#include "ir.h"
#include "indexvar.h"

namespace simit {
namespace ir {

std::vector<IndexVar> getFreeVars(Expr expr);
std::vector<IndexVar> getReductionVars(Expr expr);

bool containsFreeVar(Expr expr);
bool containsFreeVar(Stmt stmt);

bool containsReductionVar(Expr expr);
bool containsReductionVar(Stmt stmt);

/// Returns true if the statement has been flattened (only contains one index
/// expression), and false otherwise.
bool isFlattened(Stmt stmt);


/// Returns true if it is an assignment, tensor write or field write, whose
/// rhs is a blocked tensor
bool isBlocked(Stmt stmt);


std::vector<Func> getCallTree(Func func);

}}

#endif
