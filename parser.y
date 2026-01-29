%{
#include<stdio.h>
#include "gmp.h"
#include "bison.gmp_expr.h"

int yylex();
int yyerror(const char *s);

#define STACK_SIZE 20
static mpz_t stack[STACK_SIZE]; // this static data is not thread-safe
static mpz_ptr stack_ptr = stack[0];

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
%token OPEN_B CLOSE_B ADD SUB MUL DIV EXP

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
term: NUMBER { if (stack_ptr - stack[0] >= STACK_SIZE - 1) yyerror("expression overflow"); 
               stack_ptr++; mpz_set_str(stack_ptr, $1, 0); }
 | OPEN_B sum CLOSE_B
 ;

%%

extern char* gmp_expr_lex_data;

void mpz_expression_parse(mpz_t n, char *str)
{
 for (unsigned i = 0; i < STACK_SIZE; i++) mpz_init (stack[i]);
 stack_ptr = stack[0];
 gmp_expr_lex_data = str;
 yyparse();
 mpz_set(n , stack[0]);
 for (unsigned i = 0; i < STACK_SIZE; i++) mpz_clear (stack[i]);
}

int yyerror(const char *s)
{
 fprintf(stderr,"error: %s\n",s);
 return 0;
}

