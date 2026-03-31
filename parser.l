%option prefix="gmp_expr"
%option noyywrap

%{

    #include <stdlib.h>
    #include "bison.gmp_expr.tab.h" 

int gmp_exprlex_destroy(void);

static char *gmp_expr_lex_data = 0;
void gmp_expr_lex_init(char * str)
{
    gmp_expr_lex_data = str;
}

void gmp_expr_lex_terminate(void)
{
gmp_exprlex_destroy();
}

int myinput (char *buf, int buflen) {
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
("0x"[0-9a-fA-F]+)|([0-9]+)  { gmp_exprlval.str = yytext; return NUMBER; }
[ \t\f]   { /* ignore whitespace */ }
.   { printf("invalid character %c\n", *yytext); }
%%


