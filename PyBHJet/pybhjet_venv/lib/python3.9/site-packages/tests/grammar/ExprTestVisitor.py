from collections import namedtuple

from gsl.antlr import ParseTreeVisitor
from .ExprTestParser import ExprTestParser




class ExprTestVisitor(ParseTreeVisitor):
    def visitExpr(self, ctx: ExprTestParser.ExprContext):
        return self.visitNode(self.get_child(ctx, ExprTestParser.Expr0Context))

    def visitExpr0(self, ctx: ExprTestParser.Expr0Context):
        return self.visitNodes(self.get_children(ctx))

    def visitSummand(self, ctx: ExprTestParser.SummandContext):
        return self.visitNodes(self.get_children(ctx))

    def visitNumber(self, ctx: ExprTestParser.NumberContext):
        return self.visitNode(self.get_child(ctx))

    def visitSubexpr(self, ctx: ExprTestParser.SubexprContext):
        return self.visitNode(self.get_child(ctx, ExprTestParser.Expr0Context))

