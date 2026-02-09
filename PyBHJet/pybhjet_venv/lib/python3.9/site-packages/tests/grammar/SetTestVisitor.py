from collections import namedtuple

from gsl.antlr import ParseTreeVisitor
from .SetTestParser import SetTestParser




class SetTestVisitor(ParseTreeVisitor):
    def visitExpr(self, ctx: SetTestParser.ExprContext):
        return self.visitNode(self.get_child(ctx, SetTestParser.Set0Context))

    def visitSet0(self, ctx: SetTestParser.Set0Context):
        return self.visitNodes(self.get_children(ctx, SetTestParser.ElementContext))

    def visitElement(self, ctx: SetTestParser.ElementContext):
        return self.visitNode(self.get_child(ctx))

    def visitIntElement(self, ctx: SetTestParser.IntElementContext):
        return self.visitNode(self.get_child(ctx))

