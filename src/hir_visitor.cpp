#include "hir_visitor.h"
#include "hir.h"

namespace simit {
namespace hir {

void HIRVisitor::visit(Program::Ptr program) {
  for (auto elem : program->elems) {
    elem->accept(this);
  }
}

void HIRVisitor::visit(StmtBlock::Ptr stmtBlock) {
  for (auto stmt : stmtBlock->stmts) {
    stmt->accept(this);
  }
}

void HIRVisitor::visit(SetType::Ptr type) {
  type->element->accept(this);
  for (auto endpoint : type->endpoints) {
    endpoint->accept(this);
  }
}

void HIRVisitor::visit(TupleType::Ptr type) {
  type->element->accept(this);
}

void HIRVisitor::visit(NonScalarTensorType::Ptr type) {
  for (auto indexSet : type->indexSets) {
    indexSet->accept(this);
  }
  type->blockType->accept(this);
}

void HIRVisitor::visit(Field::Ptr field) {
  field->type->accept(this);
}

void HIRVisitor::visit(ElementTypeDecl::Ptr decl) {
  for (auto field : decl->fields) {
    field->accept(this);
  }
}

void HIRVisitor::visit(IdentDecl::Ptr decl) {
  decl->type->accept(this);
}

void HIRVisitor::visit(Argument::Ptr arg) {
  visit(static_cast<IdentDecl::Ptr>(arg));
}

void HIRVisitor::visit(ExternDecl::Ptr decl) {
  decl->var->accept(this);
}

void HIRVisitor::visit(FuncDecl::Ptr decl) {
  for (auto arg : decl->args) {
    arg->accept(this);
  }
  for (auto result : decl->results) {
    result->accept(this);
  }
  decl->body->accept(this);
}

void HIRVisitor::visit(ProcDecl::Ptr decl) {
  visit(static_cast<FuncDecl::Ptr>(decl)); 
}

void HIRVisitor::visit(VarDecl::Ptr decl) {
  decl->var->accept(this);
  if (decl->initVal) {
    decl->initVal->accept(this);
  }
}

void HIRVisitor::visit(ConstDecl::Ptr decl) {
  visit(static_cast<VarDecl::Ptr>(decl));
}

void HIRVisitor::visit(WhileStmt::Ptr stmt) {
  stmt->cond->accept(this);
  stmt->body->accept(this);
}

void HIRVisitor::visit(DoWhileStmt::Ptr stmt) {
  visit(static_cast<WhileStmt::Ptr>(stmt));
}

void HIRVisitor::visit(IfStmt::Ptr stmt) {
  stmt->cond->accept(this);
  stmt->ifBody->accept(this);
  if (stmt->elseBody) {
    stmt->elseBody->accept(this);
  }
}

void HIRVisitor::visit(IndexSetDomain::Ptr domain) {
  domain->domain->accept(this);
}

void HIRVisitor::visit(RangeDomain::Ptr domain) {
  domain->lower->accept(this);
  domain->upper->accept(this);
}

void HIRVisitor::visit(ForStmt::Ptr stmt) {
  stmt->domain->accept(this);
  stmt->body->accept(this);
}

void HIRVisitor::visit(PrintStmt::Ptr stmt) {
  stmt->expr->accept(this);
}

void HIRVisitor::visit(ExprStmt::Ptr stmt) {
  stmt->expr->accept(this);
}

void HIRVisitor::visit(AssignStmt::Ptr stmt) {
  for (auto lhs : stmt->lhs) {
    lhs->accept(this);
  }
  visit(static_cast<ExprStmt::Ptr>(stmt));
}

void HIRVisitor::visit(ExprParam::Ptr param) {
  param->expr->accept(this);
}

void HIRVisitor::visit(MapExpr::Ptr expr) {
  for (auto param : expr->partialActuals) {
    param->accept(this);
  }
}

void HIRVisitor::visit(UnaryExpr::Ptr expr) {
  expr->operand->accept(this);
}

void HIRVisitor::visit(BinaryExpr::Ptr expr) {
  expr->lhs->accept(this);
  expr->rhs->accept(this);
}

void HIRVisitor::visit(NaryExpr::Ptr expr) {
  for (auto operand : expr->operands) {
    operand->accept(this);
  }
}

void HIRVisitor::visit(OrExpr::Ptr expr) {
  visit(static_cast<BinaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(AndExpr::Ptr expr) {
  visit(static_cast<BinaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(XorExpr::Ptr expr) {
  visit(static_cast<BinaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(EqExpr::Ptr expr) {
  visit(static_cast<NaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(NotExpr::Ptr expr) {
  visit(static_cast<UnaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(AddExpr::Ptr expr) {
  visit(static_cast<BinaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(SubExpr::Ptr expr) {
  visit(static_cast<BinaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(MulExpr::Ptr expr) {
  visit(static_cast<BinaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(DivExpr::Ptr expr) {
  visit(static_cast<BinaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(ElwiseMulExpr::Ptr expr) {
  visit(static_cast<BinaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(ElwiseDivExpr::Ptr expr) {
  visit(static_cast<BinaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(NegExpr::Ptr expr) {
  visit(static_cast<UnaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(ExpExpr::Ptr expr) {
  visit(static_cast<BinaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(TransposeExpr::Ptr expr) {
  visit(static_cast<UnaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(CallExpr::Ptr expr) {
  visit(static_cast<NaryExpr::Ptr>(expr));
}

void HIRVisitor::visit(TensorReadExpr::Ptr expr) {
  expr->tensor->accept(this);
  for (auto param : expr->indices) {
    param->accept(this);
  }
}

void HIRVisitor::visit(FieldReadExpr::Ptr expr) {
  expr->setOrElem->accept(this);
}

void HIRVisitor::visit(DenseNDTensorLiteral::Ptr tensor) {
  for (auto elem : tensor->elems) {
    elem->accept(this);
  }
}

void HIRVisitor::visit(Test::Ptr test) {
  for (auto arg : test->args) {
    arg->accept(this);
  }
  test->expected->accept(this);
}

}
}

