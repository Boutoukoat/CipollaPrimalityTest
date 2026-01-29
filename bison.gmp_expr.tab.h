/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_GMP_EXPR_BISON_GMP_EXPR_TAB_H_INCLUDED
# define YY_GMP_EXPR_BISON_GMP_EXPR_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef GMP_EXPRDEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define GMP_EXPRDEBUG 1
#  else
#   define GMP_EXPRDEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define GMP_EXPRDEBUG 0
# endif /* ! defined YYDEBUG */
#endif  /* ! defined GMP_EXPRDEBUG */
#if GMP_EXPRDEBUG
extern int gmp_exprdebug;
#endif

/* Token kinds.  */
#ifndef GMP_EXPRTOKENTYPE
# define GMP_EXPRTOKENTYPE
  enum gmp_exprtokentype
  {
    GMP_EXPREMPTY = -2,
    GMP_EXPREOF = 0,               /* "end of file"  */
    GMP_EXPRerror = 256,           /* error  */
    GMP_EXPRUNDEF = 257,           /* "invalid token"  */
    NUMBER = 258,                  /* NUMBER  */
    OPEN_B = 259,                  /* OPEN_B  */
    CLOSE_B = 260,                 /* CLOSE_B  */
    ADD = 261,                     /* ADD  */
    SUB = 262,                     /* SUB  */
    MUL = 263,                     /* MUL  */
    DIV = 264,                     /* DIV  */
    EXP = 265                      /* EXP  */
  };
  typedef enum gmp_exprtokentype gmp_exprtoken_kind_t;
#endif

/* Value type.  */
#if ! defined GMP_EXPRSTYPE && ! defined GMP_EXPRSTYPE_IS_DECLARED
union GMP_EXPRSTYPE
{
#line 19 "parser.y"

char *str;
int lval;

#line 87 "bison.gmp_expr.tab.h"

};
typedef union GMP_EXPRSTYPE GMP_EXPRSTYPE;
# define GMP_EXPRSTYPE_IS_TRIVIAL 1
# define GMP_EXPRSTYPE_IS_DECLARED 1
#endif


extern GMP_EXPRSTYPE gmp_exprlval;


int gmp_exprparse (void);


#endif /* !YY_GMP_EXPR_BISON_GMP_EXPR_TAB_H_INCLUDED  */
