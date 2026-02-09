# Generated from tests/grammar/SetTest.g4 by ANTLR 4.7
from antlr4 import *
from io import StringIO
from typing.io import TextIO
import sys


def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\2\7")
        buf.write("\37\b\1\4\2\t\2\4\3\t\3\4\4\t\4\4\5\t\5\4\6\t\6\3\2\3")
        buf.write("\2\3\3\3\3\3\4\3\4\3\5\6\5\25\n\5\r\5\16\5\26\3\6\6\6")
        buf.write("\32\n\6\r\6\16\6\33\3\6\3\6\2\2\7\3\3\5\4\7\5\t\6\13\7")
        buf.write("\3\2\4\3\2\62;\4\2\13\13\"\"\2 \2\3\3\2\2\2\2\5\3\2\2")
        buf.write("\2\2\7\3\2\2\2\2\t\3\2\2\2\2\13\3\2\2\2\3\r\3\2\2\2\5")
        buf.write("\17\3\2\2\2\7\21\3\2\2\2\t\24\3\2\2\2\13\31\3\2\2\2\r")
        buf.write("\16\7.\2\2\16\4\3\2\2\2\17\20\7}\2\2\20\6\3\2\2\2\21\22")
        buf.write("\7\177\2\2\22\b\3\2\2\2\23\25\t\2\2\2\24\23\3\2\2\2\25")
        buf.write("\26\3\2\2\2\26\24\3\2\2\2\26\27\3\2\2\2\27\n\3\2\2\2\30")
        buf.write("\32\t\3\2\2\31\30\3\2\2\2\32\33\3\2\2\2\33\31\3\2\2\2")
        buf.write("\33\34\3\2\2\2\34\35\3\2\2\2\35\36\b\6\2\2\36\f\3\2\2")
        buf.write("\2\5\2\26\33\3\2\3\2")
        return buf.getvalue()


class SetTestLexer(Lexer):

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    COMMA = 1
    LBRACE = 2
    RBRACE = 3
    NUMBER = 4
    WS = 5

    channelNames = [ u"DEFAULT_TOKEN_CHANNEL", u"HIDDEN" ]

    modeNames = [ "DEFAULT_MODE" ]

    literalNames = [ "<INVALID>",
            "','", "'{'", "'}'" ]

    symbolicNames = [ "<INVALID>",
            "COMMA", "LBRACE", "RBRACE", "NUMBER", "WS" ]

    ruleNames = [ "COMMA", "LBRACE", "RBRACE", "NUMBER", "WS" ]

    grammarFileName = "SetTest.g4"

    def __init__(self, input=None, output:TextIO = sys.stdout):
        super().__init__(input, output)
        self.checkVersion("4.7")
        self._interp = LexerATNSimulator(self, self.atn, self.decisionsToDFA, PredictionContextCache())
        self._actions = None
        self._predicates = None


