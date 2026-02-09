from collections import namedtuple

from gsl.antlr import ParseTreeVisitor
from .HedgehogTestParser import HedgehogTestParser


Message = namedtuple('Message', ('messageType', 'discriminator', 'label', 'docstring', 'fields', 'messageClasses',))
Field = namedtuple('Field', ('repeated', 'fieldType', 'name', 'label', 'languageFieldSpecs',))
Oneof = namedtuple('Oneof', ('name', 'fields',))


class HedgehogTestVisitor(ParseTreeVisitor):
    def visitExpr(self, ctx: HedgehogTestParser.ExprContext):
        return self.visitNodes(self.get_children(ctx, HedgehogTestParser.MessageContext))

    def visitMessage(self, ctx: HedgehogTestParser.MessageContext):
        return Message(
            self.visitNode(ctx.messageType),
            self.visitNode(ctx.discriminator),
            self.visitNode(ctx.label),
            self.visitNode(ctx.DOCSTRING()) if ctx.DOCSTRING() else None,
            self.visitNodes(self.get_children(ctx, HedgehogTestParser.FieldContext, HedgehogTestParser.OneofContext)),
            self.visitNodes(self.get_children(ctx, HedgehogTestParser.MessageClassContext)),
        )

    def visitField(self, ctx: HedgehogTestParser.FieldContext):
        return Field(
            bool(ctx.rep),
            self.visitNode(ctx.fieldType),
            self.visitNode(ctx.name),
            self.visitNode(ctx.label),
            self.visitNodes(self.get_children(ctx, HedgehogTestParser.LanguageFieldSpecContext)),
        )

    def visitOneof(self, ctx: HedgehogTestParser.OneofContext):
        return Oneof(
            self.visitNode(ctx.name),
            self.visitNodes(self.get_children(ctx, HedgehogTestParser.FieldContext)),
        )

    def visitQualifiedIdentifier(self, ctx: HedgehogTestParser.QualifiedIdentifierContext):
        return self.visitNodes(self.get_children(ctx, HedgehogTestParser.IdentifierContext))

    def visitIdentifier(self, ctx: HedgehogTestParser.IdentifierContext):
        return self.visitNode(self.get_child(ctx))

    def visitNumber(self, ctx: HedgehogTestParser.NumberContext):
        return self.visitNode(self.get_child(ctx))

