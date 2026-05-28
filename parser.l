%option prefix="gmp_expr"
%option noyywrap
%option batch

%{
    #include <stdlib.h>
    #include "bison.gmp_expr.tab.h" 

int gmp_exprlex_destroy(void);

static char *gmp_expr_lex_data = 0;
static bool gmp_expr_valid_data = false;

void gmp_expr_lex_init(char * str)
{
    gmp_expr_valid_data = true;
    gmp_expr_lex_data = str;
}

bool gmp_expr_lex_terminate(void)
{
    gmp_exprlex_destroy();
    return gmp_expr_valid_data;
}

int myinput (char *buf, int buflen) 
{
    strncpy(buf, gmp_expr_lex_data, buflen - 1);
    buf[buflen-1] = 0;
    int i = strlen(buf);
    gmp_expr_lex_data += i;
    return i;
}

#define YY_INPUT(buf, result, buflen) (result = myinput(buf, buflen))

%}
%%
"+" { return ADD; }
"-" { return SUB; }
"*" { return MUL; }
"/" { return DIV; }
"^" { return EXP; }
"(" { return OPEN_B; }
")" { return CLOSE_B; }
(("0X"|"0x")[0-9a-fA-F]+)|("0")|([1-9][0-9]*)|("0"[0-7]+)|(("0B"|"0b")[0-1]+)  { gmp_exprlval.str = yytext; return NUMBER; /* hexadecimal, decimal, octal and binary */ }
[ \t\f\n\r]   { /* ignore whitespace */ }
.   { printf("invalid character %c\n", *yytext); gmp_expr_valid_data = false; return INVALID; }
%%


