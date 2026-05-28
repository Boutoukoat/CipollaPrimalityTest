%{
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "gmp.h"
#include "bison.gmp_expr.h"

void gmp_expr_lex_init(char * str);
bool gmp_expr_lex_terminate(void);

int yylex();
int yyerror(const char *s);

#define STACK_SIZE 20
static mpz_t stack[STACK_SIZE]; // this static data is not thread-safe
static mpz_ptr stack_ptr = stack[0];
static bool gmp_expr_valid = false;

%}

%define api.prefix {gmp_expr}
%file-prefix "bison.gmp_expr"

/* token type */
%union {
char *str;
int lval;
}

/* token declaration */
%token <str> NUMBER
%token OPEN_B CLOSE_B ADD SUB MUL DIV EXP INVALID

%%

calclist:
 sum { if (mpz_cmp_ui(stack_ptr, 0) < 0) yyerror("negative expression"); 
       mpz_set(stack[0], stack_ptr); }
 ;
sum: factor
 | sum ADD factor { stack_ptr--; mpz_add(stack_ptr, stack_ptr, stack_ptr+1); }
 | sum SUB factor { stack_ptr--; mpz_sub(stack_ptr, stack_ptr, stack_ptr+1); }
 ;
factor: exp
 | factor MUL exp { stack_ptr--; mpz_mul(stack_ptr, stack_ptr, stack_ptr+1); }
 | factor DIV exp { if (mpz_cmp_ui(stack_ptr, 0) == 0) yyerror("division by 0"); 
                     stack_ptr--; mpz_fdiv_q(stack_ptr, stack_ptr, stack_ptr+1); }
 ;
exp: sterm 
 | sterm EXP exp { stack_ptr--;
                  if (mpz_cmp_si(stack_ptr, 2) >= 0 || mpz_cmp_si(stack_ptr, -2) <= 0)
                  {
                      if (mpz_cmp_ui(stack_ptr+1, 0) < 0) yyerror("negative exponent"); 
                  }
                  else
                  {
			mpz_set_si(stack_ptr+1, mpz_cmp_ui(stack_ptr+1, 0));
                  }
                  if (mpz_sizeinbase(stack_ptr+1, 2) >= 64) yyerror("exponentiation overflow"); 
                  mpz_pow_ui(stack_ptr, stack_ptr, mpz_get_ui(stack_ptr+1)); }
 ;
sterm: term
 | ADD term { }
 | SUB term { mpz_neg(stack_ptr, stack_ptr); }
 ;
term: INVALID { yyerror("invalid input"); }
 | NUMBER { if (stack_ptr - stack[0] >= STACK_SIZE - 1) yyerror("expression overflow"); 
               stack_ptr++; int v = mpz_set_str(stack_ptr, $1, 0); if (v) yyerror("invalid number representation"); }
 | OPEN_B sum CLOSE_B
 ;

%%

#line 1 "number.c"

static int number_parse(const char *s)
{
	const char *YYMARKER = 0;

#line 9 "<stdout>"
{
	char yych;
	yych = *s;
	switch (yych) {
		case '0': goto yy3;
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9': goto yy4;
		default: goto yy1;
	}
yy1:
	++s;
yy2:
#line 11 "number.c"
	{ return 1; }
#line 31 "<stdout>"
yy3:
	yych = *(YYMARKER = ++s);
	switch (yych) {
		case 0x00: goto yy5;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7': goto yy6;
		case 'B':
		case 'b': goto yy8;
		case 'X':
		case 'x': goto yy9;
		default: goto yy2;
	}
yy4:
	yych = *(YYMARKER = ++s);
	switch (yych) {
		case 0x00: goto yy5;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9': goto yy10;
		default: goto yy2;
	}
yy5:
	++s;
#line 10 "number.c"
	{ return 0; }
#line 70 "<stdout>"
yy6:
	yych = *++s;
	switch (yych) {
		case 0x00: goto yy5;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7': goto yy6;
		default: goto yy7;
	}
yy7:
	s = YYMARKER;
	goto yy2;
yy8:
	yych = *++s;
	if (yych <= 0x00) goto yy7;
	goto yy12;
yy9:
	yych = *++s;
	if (yych <= 0x00) goto yy7;
	goto yy14;
yy10:
	yych = *++s;
	switch (yych) {
		case 0x00: goto yy5;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9': goto yy10;
		default: goto yy7;
	}
yy11:
	yych = *++s;
yy12:
	switch (yych) {
		case 0x00: goto yy5;
		case '0':
		case '1': goto yy11;
		default: goto yy7;
	}
yy13:
	yych = *++s;
yy14:
	switch (yych) {
		case 0x00: goto yy5;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
		case 'A':
		case 'B':
		case 'C':
		case 'D':
		case 'E':
		case 'F':
		case 'a':
		case 'b':
		case 'c':
		case 'd':
		case 'e':
		case 'f': goto yy13;
		default: goto yy7;
	}
}
#line 12 "number.c"


    return 2;
}

bool mpz_expression_parse(mpz_t n, char *str)
{
 gmp_expr_valid = true;
// check a valid single number
 int v = number_parse(str);
 if (!v) v = mpz_set_str(n, str, 0); 
if (v) 
{
// check an expression
 for (unsigned i = 0; i < STACK_SIZE; i++) { mpz_init(stack[i]); }
 stack_ptr = stack[0];
gmp_expr_lex_init(str);
 yyparse();
 bool valid = gmp_expr_lex_terminate();
 mpz_set(n , stack[0]);
 for (unsigned i = 0; i < STACK_SIZE; i++) { mpz_clear (stack[i]); }
 if (!valid) yyerror("expression is not valid"); 
}
return gmp_expr_valid;
}

int yyerror(const char *s)
{
 fprintf(stderr,"Parsing error: %s\n",s);
 gmp_expr_valid = false;
 return 1;
}


