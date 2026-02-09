# Generated from tests/grammar/ExprTest.g4 by ANTLR 4.7
# encoding: utf-8
from antlr4 import *
from io import StringIO
from typing.io import TextIO
import sys

def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\3\n")
        buf.write("%\4\2\t\2\4\3\t\3\4\4\t\4\4\5\t\5\3\2\3\2\3\2\3\3\3\3")
        buf.write("\3\3\7\3\21\n\3\f\3\16\3\24\13\3\3\4\3\4\3\4\7\4\31\n")
        buf.write("\4\f\4\16\4\34\13\4\3\5\3\5\3\5\3\5\3\5\5\5#\n\5\3\5\2")
        buf.write("\2\6\2\4\6\b\2\4\3\2\3\4\3\2\5\6\2#\2\n\3\2\2\2\4\r\3")
        buf.write("\2\2\2\6\25\3\2\2\2\b\"\3\2\2\2\n\13\5\4\3\2\13\f\7\2")
        buf.write("\2\3\f\3\3\2\2\2\r\22\5\6\4\2\16\17\t\2\2\2\17\21\5\6")
        buf.write("\4\2\20\16\3\2\2\2\21\24\3\2\2\2\22\20\3\2\2\2\22\23\3")
        buf.write("\2\2\2\23\5\3\2\2\2\24\22\3\2\2\2\25\32\5\b\5\2\26\27")
        buf.write("\t\3\2\2\27\31\5\b\5\2\30\26\3\2\2\2\31\34\3\2\2\2\32")
        buf.write("\30\3\2\2\2\32\33\3\2\2\2\33\7\3\2\2\2\34\32\3\2\2\2\35")
        buf.write("#\7\t\2\2\36\37\7\7\2\2\37 \5\4\3\2 !\7\b\2\2!#\3\2\2")
        buf.write("\2\"\35\3\2\2\2\"\36\3\2\2\2#\t\3\2\2\2\5\22\32\"")
        return buf.getvalue()


class ExprTestParser ( Parser ):

    grammarFileName = "ExprTest.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "'+'", "'-'", "'*'", "'/'", "'('", "')'" ]

    symbolicNames = [ "<INVALID>", "ADD", "SUB", "MUL", "DIV", "LPAR", "RPAR", 
                      "NUMBER", "WS" ]

    RULE_expr = 0
    RULE_expr0 = 1
    RULE_summand = 2
    RULE_factor = 3

    ruleNames =  [ "expr", "expr0", "summand", "factor" ]

    EOF = Token.EOF
    ADD=1
    SUB=2
    MUL=3
    DIV=4
    LPAR=5
    RPAR=6
    NUMBER=7
    WS=8

    def __init__(self, input:TokenStream, output:TextIO = sys.stdout):
        super().__init__(input, output)
        self.checkVersion("4.7")
        self._interp = ParserATNSimulator(self, self.atn, self.decisionsToDFA, self.sharedContextCache)
        self._predicates = None



    class ExprContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def expr0(self):
            return self.getTypedRuleContext(ExprTestParser.Expr0Context,0)


        def EOF(self):
            return self.getToken(ExprTestParser.EOF, 0)

        def getRuleIndex(self):
            return ExprTestParser.RULE_expr

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitExpr" ):
                return visitor.visitExpr(self)
            else:
                return visitor.visitChildren(self)




    def expr(self):

        localctx = ExprTestParser.ExprContext(self, self._ctx, self.state)
        self.enterRule(localctx, 0, self.RULE_expr)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 8
            self.expr0()
            self.state = 9
            self.match(ExprTestParser.EOF)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class Expr0Context(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def summand(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(ExprTestParser.SummandContext)
            else:
                return self.getTypedRuleContext(ExprTestParser.SummandContext,i)


        def ADD(self, i:int=None):
            if i is None:
                return self.getTokens(ExprTestParser.ADD)
            else:
                return self.getToken(ExprTestParser.ADD, i)

        def SUB(self, i:int=None):
            if i is None:
                return self.getTokens(ExprTestParser.SUB)
            else:
                return self.getToken(ExprTestParser.SUB, i)

        def getRuleIndex(self):
            return ExprTestParser.RULE_expr0

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitExpr0" ):
                return visitor.visitExpr0(self)
            else:
                return visitor.visitChildren(self)




    def expr0(self):

        localctx = ExprTestParser.Expr0Context(self, self._ctx, self.state)
        self.enterRule(localctx, 2, self.RULE_expr0)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 11
            self.summand()
            self.state = 16
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while _la==ExprTestParser.ADD or _la==ExprTestParser.SUB:
                self.state = 12
                _la = self._input.LA(1)
                if not(_la==ExprTestParser.ADD or _la==ExprTestParser.SUB):
                    self._errHandler.recoverInline(self)
                else:
                    self._errHandler.reportMatch(self)
                    self.consume()
                self.state = 13
                self.summand()
                self.state = 18
                self._errHandler.sync(self)
                _la = self._input.LA(1)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class SummandContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def factor(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(ExprTestParser.FactorContext)
            else:
                return self.getTypedRuleContext(ExprTestParser.FactorContext,i)


        def MUL(self, i:int=None):
            if i is None:
                return self.getTokens(ExprTestParser.MUL)
            else:
                return self.getToken(ExprTestParser.MUL, i)

        def DIV(self, i:int=None):
            if i is None:
                return self.getTokens(ExprTestParser.DIV)
            else:
                return self.getToken(ExprTestParser.DIV, i)

        def getRuleIndex(self):
            return ExprTestParser.RULE_summand

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitSummand" ):
                return visitor.visitSummand(self)
            else:
                return visitor.visitChildren(self)




    def summand(self):

        localctx = ExprTestParser.SummandContext(self, self._ctx, self.state)
        self.enterRule(localctx, 4, self.RULE_summand)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 19
            self.factor()
            self.state = 24
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while _la==ExprTestParser.MUL or _la==ExprTestParser.DIV:
                self.state = 20
                _la = self._input.LA(1)
                if not(_la==ExprTestParser.MUL or _la==ExprTestParser.DIV):
                    self._errHandler.recoverInline(self)
                else:
                    self._errHandler.reportMatch(self)
                    self.consume()
                self.state = 21
                self.factor()
                self.state = 26
                self._errHandler.sync(self)
                _la = self._input.LA(1)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class FactorContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser


        def getRuleIndex(self):
            return ExprTestParser.RULE_factor

     
        def copyFrom(self, ctx:ParserRuleContext):
            super().copyFrom(ctx)



    class NumberContext(FactorContext):

        def __init__(self, parser, ctx:ParserRuleContext): # actually a ExprTestParser.FactorContext
            super().__init__(parser)
            self.copyFrom(ctx)

        def NUMBER(self):
            return self.getToken(ExprTestParser.NUMBER, 0)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitNumber" ):
                return visitor.visitNumber(self)
            else:
                return visitor.visitChildren(self)


    class SubexprContext(FactorContext):

        def __init__(self, parser, ctx:ParserRuleContext): # actually a ExprTestParser.FactorContext
            super().__init__(parser)
            self.copyFrom(ctx)

        def LPAR(self):
            return self.getToken(ExprTestParser.LPAR, 0)
        def expr0(self):
            return self.getTypedRuleContext(ExprTestParser.Expr0Context,0)

        def RPAR(self):
            return self.getToken(ExprTestParser.RPAR, 0)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitSubexpr" ):
                return visitor.visitSubexpr(self)
            else:
                return visitor.visitChildren(self)



    def factor(self):

        localctx = ExprTestParser.FactorContext(self, self._ctx, self.state)
        self.enterRule(localctx, 6, self.RULE_factor)
        try:
            self.state = 32
            self._errHandler.sync(self)
            token = self._input.LA(1)
            if token in [ExprTestParser.NUMBER]:
                localctx = ExprTestParser.NumberContext(self, localctx)
                self.enterOuterAlt(localctx, 1)
                self.state = 27
                self.match(ExprTestParser.NUMBER)
                pass
            elif token in [ExprTestParser.LPAR]:
                localctx = ExprTestParser.SubexprContext(self, localctx)
                self.enterOuterAlt(localctx, 2)
                self.state = 28
                self.match(ExprTestParser.LPAR)
                self.state = 29
                self.expr0()
                self.state = 30
                self.match(ExprTestParser.RPAR)
                pass
            else:
                raise NoViableAltException(self)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx





