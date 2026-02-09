# Generated from tests/grammar/SetTest.g4 by ANTLR 4.7
# encoding: utf-8
from antlr4 import *
from io import StringIO
from typing.io import TextIO
import sys

def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\3\7")
        buf.write("!\4\2\t\2\4\3\t\3\4\4\t\4\4\5\t\5\3\2\3\2\3\2\3\3\3\3")
        buf.write("\3\3\3\3\7\3\22\n\3\f\3\16\3\25\13\3\5\3\27\n\3\3\3\3")
        buf.write("\3\3\4\3\4\5\4\35\n\4\3\5\3\5\3\5\2\2\6\2\4\6\b\2\2\2")
        buf.write("\37\2\n\3\2\2\2\4\r\3\2\2\2\6\34\3\2\2\2\b\36\3\2\2\2")
        buf.write("\n\13\5\4\3\2\13\f\7\2\2\3\f\3\3\2\2\2\r\26\7\4\2\2\16")
        buf.write("\23\5\6\4\2\17\20\7\3\2\2\20\22\5\6\4\2\21\17\3\2\2\2")
        buf.write("\22\25\3\2\2\2\23\21\3\2\2\2\23\24\3\2\2\2\24\27\3\2\2")
        buf.write("\2\25\23\3\2\2\2\26\16\3\2\2\2\26\27\3\2\2\2\27\30\3\2")
        buf.write("\2\2\30\31\7\5\2\2\31\5\3\2\2\2\32\35\5\b\5\2\33\35\5")
        buf.write("\4\3\2\34\32\3\2\2\2\34\33\3\2\2\2\35\7\3\2\2\2\36\37")
        buf.write("\7\6\2\2\37\t\3\2\2\2\5\23\26\34")
        return buf.getvalue()


class SetTestParser ( Parser ):

    grammarFileName = "SetTest.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "','", "'{'", "'}'" ]

    symbolicNames = [ "<INVALID>", "COMMA", "LBRACE", "RBRACE", "NUMBER", 
                      "WS" ]

    RULE_expr = 0
    RULE_set0 = 1
    RULE_element = 2
    RULE_simpleElement = 3

    ruleNames =  [ "expr", "set0", "element", "simpleElement" ]

    EOF = Token.EOF
    COMMA=1
    LBRACE=2
    RBRACE=3
    NUMBER=4
    WS=5

    def __init__(self, input:TokenStream, output:TextIO = sys.stdout):
        super().__init__(input, output)
        self.checkVersion("4.7")
        self._interp = ParserATNSimulator(self, self.atn, self.decisionsToDFA, self.sharedContextCache)
        self._predicates = None



    class ExprContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def set0(self):
            return self.getTypedRuleContext(SetTestParser.Set0Context,0)


        def EOF(self):
            return self.getToken(SetTestParser.EOF, 0)

        def getRuleIndex(self):
            return SetTestParser.RULE_expr

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitExpr" ):
                return visitor.visitExpr(self)
            else:
                return visitor.visitChildren(self)




    def expr(self):

        localctx = SetTestParser.ExprContext(self, self._ctx, self.state)
        self.enterRule(localctx, 0, self.RULE_expr)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 8
            self.set0()
            self.state = 9
            self.match(SetTestParser.EOF)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class Set0Context(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def LBRACE(self):
            return self.getToken(SetTestParser.LBRACE, 0)

        def RBRACE(self):
            return self.getToken(SetTestParser.RBRACE, 0)

        def element(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(SetTestParser.ElementContext)
            else:
                return self.getTypedRuleContext(SetTestParser.ElementContext,i)


        def COMMA(self, i:int=None):
            if i is None:
                return self.getTokens(SetTestParser.COMMA)
            else:
                return self.getToken(SetTestParser.COMMA, i)

        def getRuleIndex(self):
            return SetTestParser.RULE_set0

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitSet0" ):
                return visitor.visitSet0(self)
            else:
                return visitor.visitChildren(self)




    def set0(self):

        localctx = SetTestParser.Set0Context(self, self._ctx, self.state)
        self.enterRule(localctx, 2, self.RULE_set0)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 11
            self.match(SetTestParser.LBRACE)
            self.state = 20
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            if _la==SetTestParser.LBRACE or _la==SetTestParser.NUMBER:
                self.state = 12
                self.element()
                self.state = 17
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                while _la==SetTestParser.COMMA:
                    self.state = 13
                    self.match(SetTestParser.COMMA)
                    self.state = 14
                    self.element()
                    self.state = 19
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)



            self.state = 22
            self.match(SetTestParser.RBRACE)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class ElementContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def simpleElement(self):
            return self.getTypedRuleContext(SetTestParser.SimpleElementContext,0)


        def set0(self):
            return self.getTypedRuleContext(SetTestParser.Set0Context,0)


        def getRuleIndex(self):
            return SetTestParser.RULE_element

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitElement" ):
                return visitor.visitElement(self)
            else:
                return visitor.visitChildren(self)




    def element(self):

        localctx = SetTestParser.ElementContext(self, self._ctx, self.state)
        self.enterRule(localctx, 4, self.RULE_element)
        try:
            self.state = 26
            self._errHandler.sync(self)
            token = self._input.LA(1)
            if token in [SetTestParser.NUMBER]:
                self.enterOuterAlt(localctx, 1)
                self.state = 24
                self.simpleElement()
                pass
            elif token in [SetTestParser.LBRACE]:
                self.enterOuterAlt(localctx, 2)
                self.state = 25
                self.set0()
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

    class SimpleElementContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser


        def getRuleIndex(self):
            return SetTestParser.RULE_simpleElement

     
        def copyFrom(self, ctx:ParserRuleContext):
            super().copyFrom(ctx)



    class IntElementContext(SimpleElementContext):

        def __init__(self, parser, ctx:ParserRuleContext): # actually a SetTestParser.SimpleElementContext
            super().__init__(parser)
            self.copyFrom(ctx)

        def NUMBER(self):
            return self.getToken(SetTestParser.NUMBER, 0)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitIntElement" ):
                return visitor.visitIntElement(self)
            else:
                return visitor.visitChildren(self)



    def simpleElement(self):

        localctx = SetTestParser.SimpleElementContext(self, self._ctx, self.state)
        self.enterRule(localctx, 6, self.RULE_simpleElement)
        try:
            localctx = SetTestParser.IntElementContext(self, localctx)
            self.enterOuterAlt(localctx, 1)
            self.state = 28
            self.match(SetTestParser.NUMBER)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx





