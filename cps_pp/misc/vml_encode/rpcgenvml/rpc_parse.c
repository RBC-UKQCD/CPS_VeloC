/* Hacked by Peter Boyle for VML 2004 *//*
 * Sun RPC is a product of Sun Microsystems, Inc. and is provided for
 * unrestricted use provided that this legend is included on all tape
 * media and as a part of the software program in whole or part.  Users
 * may copy or modify Sun RPC without charge, but are not authorized
 * to license or distribute it to anyone else except as part of a product or
 * program developed by the user or with the express written consent of
 * Sun Microsystems, Inc.
 *
 * SUN RPC IS PROVIDED AS IS WITH NO WARRANTIES OF ANY KIND INCLUDING THE
 * WARRANTIES OF DESIGN, MERCHANTIBILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE, OR ARISING FROM A COURSE OF DEALING, USAGE OR TRADE PRACTICE.
 *
 * Sun RPC is provided with no support and without any obligation on the
 * part of Sun Microsystems, Inc. to assist in its use, correction,
 * modification or enhancement.
 *
 * SUN MICROSYSTEMS, INC. SHALL HAVE NO LIABILITY WITH RESPECT TO THE
 * INFRINGEMENT OF COPYRIGHTS, TRADE SECRETS OR ANY PATENTS BY SUN RPC
 * OR ANY PART THEREOF.
 *
 * In no event will Sun Microsystems, Inc. be liable for any lost revenue
 * or profits or other special, indirect and consequential damages, even if
 * Sun has been advised of the possibility of such damages.
 *
 * Sun Microsystems, Inc.
 * 2550 Garcia Avenue
 * Mountain View, California  94043
 */

/*
 * From: @(#)rpc_parse.c 1.8 89/02/22 (C) 1987 SMI
 */
const char parse_rcsid[] =
  "$Id: rpc_parse.c,v 1.4.358.1 2012-11-15 18:17:08 ckelly Exp $";

/*
 * rpc_parse.c, Parser for the RPC protocol compiler
 * Copyright (C) 1987 Sun Microsystems, Inc.
 */
#include <stdio.h>
#include <string.h>
#include "rpc/types.h"
#include "rpc_scan.h"
#include "rpc_parse.h"
#include "rpc_util.h"
#include "proto.h"

#define ARGNAME "arg"

static void isdefined (definition * defp);
static void def_struct (definition * defp);
static void def_program (definition * defp);
static void def_enum (definition * defp);
static void def_const (definition * defp);
static void def_union (definition * defp);
static void check_type_name (const char *name, int new_type);
static void def_typedef (definition * defp);
static void get_declaration (declaration * dec, defkind dkind);
static void get_prog_declaration (declaration * dec, defkind dkind, int num);
static void get_type (const char **prefixp, const char **typep, defkind dkind);
static void unsigned_dec (const char **typep);
/*
 * PAB Class definition parsing
 */
static void def_class (definition * defp);

/* CK include pragma parsing */
static void def_include (definition * defp);

/*
 * return the next definition you see
 */
definition *
get_definition (void)
{
  definition *defp;
  token tok;

  defp = ALLOC (definition);
  get_token (&tok);
  switch (tok.kind)
    {
    case TOK_STRUCT:
      def_struct (defp);
      break;
    case TOK_UNION:
      def_union (defp);
      break;
    case TOK_TYPEDEF:
      def_typedef (defp);
      break;
    case TOK_ENUM:
      def_enum (defp);
      break;
    case TOK_PROGRAM:
      def_program (defp);
      break;
    case TOK_CONST:
      def_const (defp);
      break;
    case TOK_CLASS: /*PAB*/
      def_class (defp);
      break;
    case TOK_INCLUDEPRAGMA: /*CK*/
      def_include (defp);
      break;
    case TOK_EOF:
      return (NULL);
    default:
      error ("definition keyword expected");
    }
  scan (TOK_SEMICOLON, &tok);
  isdefined (defp);
  return (defp);
}

static void
isdefined (definition * defp)
{
  STOREVAL (&defined, defp);
}

/*CK*/
static void
def_include (definition * defp)
{
  token tok;
  defp->def_kind = DEF_INCLUDEPRAGMA;
  
  defp->def.in.is_relative = 0;
  peek(&tok);
  
  if(tok.kind==TOK_STRCONST) defp->def.in.is_relative = 1;

  //f_print (stderr,"def_include got relative %d\n",defp->def.in.is_relative);

  if(defp->def.in.is_relative) scan(TOK_STRCONST,&tok);
  else{ 
    scan (TOK_LANGLE, &tok);
    scan (TOK_IDENT, &tok);
  }

  defp->def_name = "";
  defp->def.in.file = tok.str;

  //f_print(stderr,"def_include got string '%s'\n",tok.str);
  
  if(!defp->def.in.is_relative) scan(TOK_RANGLE,&tok);
}
/*End CK*/

static void
def_struct (definition * defp)
{
  token tok;
  declaration dec;
  decl_list *decls;
  decl_list **tailp;

  defp->def_kind = DEF_STRUCT;

  scan (TOK_IDENT, &tok);
  defp->def_name = tok.str;
  scan (TOK_LBRACE, &tok);
  tailp = &defp->def.st.decls;
  do
    {
      get_declaration (&dec, DEF_STRUCT);
      decls = ALLOC (decl_list);
      decls->decl = dec;
      *tailp = decls;
      tailp = &decls->next;
      scan (TOK_SEMICOLON, &tok);
      peek (&tok);
    }
  while (tok.kind != TOK_RBRACE);
  get_token (&tok);
  *tailp = NULL;
}

static void
def_program (definition * defp)
{
  token tok;
  declaration dec;
  decl_list *decls;
  decl_list **tailp;
  version_list *vlist;
  version_list **vtailp;
  proc_list *plist;
  proc_list **ptailp;
  int num_args;
  bool_t isvoid = FALSE;	/* whether first argument is void */
  defp->def_kind = DEF_PROGRAM;
  scan (TOK_IDENT, &tok);
  defp->def_name = tok.str;
  scan (TOK_LBRACE, &tok);
  vtailp = &defp->def.pr.versions;
  tailp = &defp->def.st.decls;
  scan (TOK_VERSION, &tok);
  do
    {
      scan (TOK_IDENT, &tok);
      vlist = ALLOC (version_list);
      vlist->vers_name = tok.str;
      scan (TOK_LBRACE, &tok);
      ptailp = &vlist->procs;
      do
	{
	  /* get result type */
	  plist = ALLOC (proc_list);
	  get_type (&plist->res_prefix, &plist->res_type,
		    DEF_PROGRAM);
	  if (streq (plist->res_type, "opaque"))
	    {
	      error ("illegal result type");
	    }
	  scan (TOK_IDENT, &tok);
	  plist->proc_name = tok.str;
	  scan (TOK_LPAREN, &tok);
	  /* get args - first one */
	  num_args = 1;
	  isvoid = FALSE;
	  /* type of DEF_PROGRAM in the first
	   * get_prog_declaration and DEF_STURCT in the next
	   * allows void as argument if it is the only argument
	   */
	  get_prog_declaration (&dec, DEF_PROGRAM, num_args);
	  if (streq (dec.type, "void"))
	    isvoid = TRUE;
	  decls = ALLOC (decl_list);
	  plist->args.decls = decls;
	  decls->decl = dec;
	  tailp = &decls->next;
	  /* get args */
	  while (peekscan (TOK_COMMA, &tok))
	    {
	      num_args++;
	      get_prog_declaration (&dec, DEF_STRUCT,
				    num_args);
	      decls = ALLOC (decl_list);
	      decls->decl = dec;
	      *tailp = decls;
	      if (streq (dec.type, "void"))
		isvoid = TRUE;
	      tailp = &decls->next;
	    }
	  /* multiple arguments are only allowed in newstyle */
	  if (!newstyle && num_args > 1)
	    {
	      error ("only one argument is allowed");
	    }
	  if (isvoid && num_args > 1)
	    {
	      error ("illegal use of void in program definition");
	    }
	  *tailp = NULL;
	  scan (TOK_RPAREN, &tok);
	  scan (TOK_EQUAL, &tok);
	  scan_num (&tok);
	  scan (TOK_SEMICOLON, &tok);
	  plist->proc_num = tok.str;
	  plist->arg_num = num_args;
	  *ptailp = plist;
	  ptailp = &plist->next;
	  peek (&tok);
	}
      while (tok.kind != TOK_RBRACE);
      *ptailp = NULL;
      *vtailp = vlist;
      vtailp = &vlist->next;
      scan (TOK_RBRACE, &tok);
      scan (TOK_EQUAL, &tok);
      scan_num (&tok);
      vlist->vers_num = tok.str;
      /* make the argument structure name for each arg */
      for (plist = vlist->procs; plist != NULL;
	   plist = plist->next)
	{
	  plist->args.argname = make_argname (plist->proc_name,
					      vlist->vers_num);
	  /* free the memory ?? */
	}
      scan (TOK_SEMICOLON, &tok);
      scan2 (TOK_VERSION, TOK_RBRACE, &tok);
    }
  while (tok.kind == TOK_VERSION);
  scan (TOK_EQUAL, &tok);
  scan_num (&tok);
  defp->def.pr.prog_num = tok.str;
  *vtailp = NULL;
}


static void
def_enum (definition * defp)
{
  token tok;
  enumval_list *elist;
  enumval_list **tailp;

  defp->def_kind = DEF_ENUM;
  scan (TOK_IDENT, &tok);
  defp->def_name = tok.str;
  scan (TOK_LBRACE, &tok);
  tailp = &defp->def.en.vals;
  do
    {
      scan (TOK_IDENT, &tok);
      elist = ALLOC (enumval_list);
      elist->name = tok.str;
      elist->assignment = NULL;
      scan3 (TOK_COMMA, TOK_RBRACE, TOK_EQUAL, &tok);
      if (tok.kind == TOK_EQUAL)
	{
	  scan_num (&tok);
	  elist->assignment = tok.str;
	  scan2 (TOK_COMMA, TOK_RBRACE, &tok);
	}
      *tailp = elist;
      tailp = &elist->next;
    }
  while (tok.kind != TOK_RBRACE);
  *tailp = NULL;
}

static void
def_const (definition * defp)
{
  token tok;

  defp->def_kind = DEF_CONST;
  scan (TOK_IDENT, &tok);
  defp->def_name = tok.str;
  scan (TOK_EQUAL, &tok);
  scan2 (TOK_IDENT, TOK_STRCONST, &tok);
  defp->def.co = tok.str;
}

static void
def_union (definition *defp)
{
  token tok;
  declaration dec;
  case_list *cases;
/*  case_list *tcase; */
  case_list **tailp;
  int flag;

  defp->def_kind = DEF_UNION;
  scan (TOK_IDENT, &tok);
  defp->def_name = tok.str;

  /*CK add optional extra open brace to enclose other declarations within 'union' struct*/
  int other_decl = 0;
  if(peekscan(TOK_LBRACE,&tok)) other_decl = 1; //consumes the token if returns true
  /* end CK */

  scan (TOK_SWITCH, &tok);
  scan (TOK_LPAREN, &tok);
  get_declaration (&dec, DEF_UNION);
  defp->def.un.enum_decl = dec;
  tailp = &defp->def.un.cases;
  scan (TOK_RPAREN, &tok);
  scan (TOK_LBRACE, &tok);
  scan (TOK_CASE, &tok);
  while (tok.kind == TOK_CASE)
    {
      scan2 (TOK_IDENT, TOK_CHARCONST, &tok);
      cases = ALLOC (case_list);
      cases->case_name = tok.str;
      scan (TOK_COLON, &tok);
      /* now peek at next token */
      flag = 0;
      if (peekscan (TOK_CASE, &tok))
	{

	  do
	    {
	      scan2 (TOK_IDENT, TOK_CHARCONST, &tok);
	      cases->contflag = 1;	/* continued case statement */
	      *tailp = cases;
	      tailp = &cases->next;
	      cases = ALLOC (case_list);
	      cases->case_name = tok.str;
	      scan (TOK_COLON, &tok);

	    }
	  while (peekscan (TOK_CASE, &tok));
	}
      else if (flag)
	{

	  *tailp = cases;
	  tailp = &cases->next;
	  cases = ALLOC (case_list);
	};

      get_declaration (&dec, DEF_UNION);
      cases->case_decl = dec;
      cases->contflag = 0;	/* no continued case statement */
      *tailp = cases;
      tailp = &cases->next;
      scan (TOK_SEMICOLON, &tok);

      scan3 (TOK_CASE, TOK_DEFAULT, TOK_RBRACE, &tok);
    }
  *tailp = NULL;
  if (tok.kind == TOK_DEFAULT)
    {
      scan (TOK_COLON, &tok);
      get_declaration (&dec, DEF_UNION);
      defp->def.un.default_decl = ALLOC (declaration);
      *defp->def.un.default_decl = dec;
      scan (TOK_SEMICOLON, &tok);
      scan (TOK_RBRACE, &tok);
    }
  else
    {
      defp->def.un.default_decl = NULL;
    }

  /*CK extra declarations*/
  if(other_decl){
    decl_list **tailp = &defp->def.un.other_decls;
    declaration dec;
    decl_list *decls;
    do
      {
	get_declaration (&dec, DEF_UNION);
	decls = ALLOC (decl_list);
	decls->decl = dec;
	*tailp = decls;
	tailp = &decls->next;
	scan (TOK_SEMICOLON, &tok);
	peek (&tok);
      }
    while (tok.kind != TOK_RBRACE);
    get_token (&tok);
    *tailp = NULL;
  }
  /*End CK*/
}

static const char *reserved_words[] =
{
  "array",
  "bytes",
  "destroy",
  "free",
  "getpos",
  "inline",
  "pointer",
  "reference",
  "setpos",
  "sizeof",
  "union",
  "vector",
  NULL
};

static const char *reserved_types[] =
{
  "opaque",
  "string",
  NULL
};

/*
 * check that the given name is not one that would eventually result in
 * xdr routines that would conflict with internal XDR routines.
 */
static void
check_type_name (const char *name, int new_type)
{
  int i;
  char tmp[100];

  for (i = 0; reserved_words[i] != NULL; i++)
    {
      if (strcmp (name, reserved_words[i]) == 0)
	{
	  sprintf (tmp,
		"illegal (reserved) name :\'%s\' in type definition", name);
	  error (tmp);
	}
    }
  if (new_type)
    {
      for (i = 0; reserved_types[i] != NULL; i++)
	{
	  if (strcmp (name, reserved_types[i]) == 0)
	    {
	      sprintf (tmp,
		"illegal (reserved) name :\'%s\' in type definition", name);
	      error (tmp);
	    }
	}
    }
}



static void
def_typedef (definition * defp)
{
  declaration dec;

  defp->def_kind = DEF_TYPEDEF;
  get_declaration (&dec, DEF_TYPEDEF);
  defp->def_name = dec.name;
  check_type_name (dec.name, 1);
  defp->def.ty.old_prefix = dec.prefix;
  defp->def.ty.old_type = dec.type;
  defp->def.ty.rel = dec.rel;
  defp->def.ty.array_max = dec.array_max;
  if ( dec.array_max && strlen(dec.array_max) ) {
    fprintf(stderr,"Typedef %s is %s[%s]\n",defp->def_name,defp->def.ty.old_type,dec.array_max);
  } else { 
    fprintf(stderr,"Typedef %s is %s\n",defp->def_name,defp->def.ty.old_type);
  }
}

static void
get_declaration (declaration * dec, defkind dkind)
{
  token tok;

  get_type (&dec->prefix, &dec->type, dkind);

  dec->array_max = "";

  dec->rel = REL_ALIAS;
  if (streq (dec->type, "void"))
    {
      return;
    }


  if ( streq (dec->type,"memfun") ) { 

    dec->name = dec->prefix;
    dec->type = "";
    dec->prefix = "";
    return;
  }

  if ( streq (dec->type,"rpccommand") ) {
    /* Added by CK */
    dec->name = dec->prefix;
    dec->prefix = "";
    return;
  }

  if ( streq (dec->type,"//") ) { 
    dec->name = dec->prefix;
    dec->type = "";
    dec->prefix = "//";
    return;
  }
  check_type_name (dec->type, 0);

  scan2 (TOK_STAR, TOK_IDENT, &tok);
  if (tok.kind == TOK_STAR)
    {
      dec->rel = REL_POINTER;
      scan (TOK_IDENT, &tok);
    }
  dec->name = tok.str;
  if (peekscan (TOK_LBRACKET, &tok))
    {
      if (dec->rel == REL_POINTER)
	{
	  error ("no array-of-pointer declarations -- use typedef");
	}
      dec->rel = REL_VECTOR;
      scan_num (&tok);
      dec->array_max = tok.str;
      scan (TOK_RBRACKET, &tok);
    }
  else if (peekscan (TOK_LANGLE, &tok))
    {
      if (dec->rel == REL_POINTER)
	{
	  error ("no array-of-pointer declarations -- use typedef");
	}
      dec->rel = REL_ARRAY;
      if (peekscan (TOK_RANGLE, &tok))
	{
	  dec->array_max = "~0";	/* unspecified size, use max */
	}
      else
	{
	  scan_num (&tok);
	  dec->array_max = tok.str;
	  scan (TOK_RANGLE, &tok);
	}
    }
  if (streq (dec->type, "opaque"))
    {
      if (dec->rel != REL_ARRAY && dec->rel != REL_VECTOR)
	{
	  error ("array declaration expected");
	}
    }
  else if (streq (dec->type, "string"))
    {
      if (dec->rel != REL_ARRAY)
	{
	  error ("variable-length array declaration expected");
	}
    }
}

static void
get_prog_declaration (declaration * dec, defkind dkind, int num /* arg number */ )
{
  token tok;
  char name[10];		/* argument name */

  if (dkind == DEF_PROGRAM)
    {
      peek (&tok);
      if (tok.kind == TOK_RPAREN)
	{			/* no arguments */
	  dec->rel = REL_ALIAS;
	  dec->type = "void";
	  dec->prefix = NULL;
	  dec->name = NULL;
	  return;
	}
    }
  get_type (&dec->prefix, &dec->type, dkind);
  dec->rel = REL_ALIAS;
  if (peekscan (TOK_IDENT, &tok))	/* optional name of argument */
    strcpy (name, tok.str);
  else
    sprintf (name, "%s%d", ARGNAME, num);	/* default name of argument */

  dec->name = (char *) strdup (name);

  if (streq (dec->type, "void"))
    {
      return;
    }

  if (streq (dec->type, "opaque"))
    {
      error ("opaque -- illegal argument type");
    }
  if (peekscan (TOK_STAR, &tok))
    {
      if (streq (dec->type, "string"))
	{
	  error ("pointer to string not allowed in program arguments\n");
	}
      dec->rel = REL_POINTER;
      if (peekscan (TOK_IDENT, &tok))	/* optional name of argument */
	dec->name = strdup (tok.str);
    }
  if (peekscan (TOK_LANGLE, &tok))
    {
      if (!streq (dec->type, "string"))
	{
	  error ("arrays cannot be declared as arguments to procedures -- use typedef");
	}
      dec->rel = REL_ARRAY;
      if (peekscan (TOK_RANGLE, &tok))
	{
	  dec->array_max = "~0";	/* unspecified size, use max */
	}
      else
	{
	  scan_num (&tok);
	  dec->array_max = tok.str;
	  scan (TOK_RANGLE, &tok);
	}
    }
  if (streq (dec->type, "string"))
    {
      if (dec->rel != REL_ARRAY)
	{			/* .x specifies just string as
				 * type of argument
				 * - make it string<>
				 */
	  dec->rel = REL_ARRAY;
	  dec->array_max = "~0";	/* unspecified size, use max */
	}
    }
}

static void
get_type (const char **prefixp, const char **typep, defkind dkind)
{
  token tok;
  char *memfun_string;

  *prefixp = NULL;
  get_token (&tok);
  switch (tok.kind)
    {
    case TOK_COMMENT:
      *typep  = "// ";
      *prefixp= tok.str;
      break;
    case TOK_IDENT:
      *typep = tok.str;
      break;
    case TOK_STRUCT:
    case TOK_CLASS:
    case TOK_ENUM:
    case TOK_UNION:
      *prefixp = tok.str;
      scan (TOK_IDENT, &tok);
      *typep = tok.str;
      break;
    case TOK_UNSIGNED:
      unsigned_dec (typep);
      break;
    case TOK_SHORT:
      *typep = "short";
      (void) peekscan (TOK_INT, &tok);
      break;
    case TOK_LONG:
      *typep = "long";
      (void) peekscan (TOK_INT, &tok);
      break;
    case TOK_VOID:
      if (dkind != DEF_UNION && dkind != DEF_PROGRAM)
	{
	  error ("voids allowed only inside union and program definitions with one argument");
	}
      *typep = tok.str;
      break;
    case TOK_STRING:
    case TOK_OPAQUE:
    case TOK_CHAR:
    case TOK_INT:
    case TOK_FLOAT:
    case TOK_DOUBLE:
    case TOK_BOOL:
      *typep = tok.str;
      break;
    case TOK_INHERITANCE:
      if (dkind != DEF_CLASS ) error("Inheritance only inside a class");
      error("Inheritance not implemented");
      *typep = tok.str;
      break;
    case TOK_MEMFUN:
      /*   if (dkind != DEF_CLASS ) error("Member functions only inside a class"); */
      memfun_string = (char *)malloc(512);
      memfun_string[0]='\0';
      *typep = tok.str;

      do { 
	get_token (&tok);
	if ( tok.kind == TOK_LPAREN ) { 
	  strcat(memfun_string," ( ");
	} else if ( tok.kind == TOK_RPAREN ) { 
	  strcat(memfun_string," ) ");
	} else if ( tok.kind == TOK_SEMICOLON ) { 
	  strcat(memfun_string," ; ");
	} else if ( tok.kind == TOK_COMMA ) { 
	  strcat(memfun_string," , ");
	} else if ( tok.kind == TOK_STAR ) {
	  strcat(memfun_string," * "); /* Added by CK to allow pointers in memfun */
	} else if ( tok.kind == TOK_AMPERSAND ) {
	  strcat(memfun_string," & "); /* Added by CK to allow references in memfun */
	} else {
	  strcat(memfun_string," ");
	  strcat(memfun_string,tok.str);
	}
      } while ( tok.kind != TOK_RPAREN);

      *prefixp=memfun_string;
      break;
    case TOK_RPCCOMMAND:
      /* Added by CK to allow user to ask for specific hard-coded operations to be added to the source code */
      *typep = tok.str;
      get_token (&tok);
      *prefixp=tok.str;
      break;  
    default:
      error ("expected type specifier");
    }
}

static void
unsigned_dec (const char **typep)
{
  token tok;

  peek (&tok);
  switch (tok.kind)
    {
    case TOK_CHAR:
      get_token (&tok);
      *typep = "u_char";
      break;
    case TOK_SHORT:
      get_token (&tok);
      *typep = "u_short";
      (void) peekscan (TOK_INT, &tok);
      break;
    case TOK_LONG:
      get_token (&tok);
      *typep = "u_long";
      (void) peekscan (TOK_INT, &tok);
      break;
    case TOK_INT:
      get_token (&tok);
      *typep = "u_int";
      break;
    default:
      *typep = "u_int";
      break;
    }
}

/*
 * PAB ... parse a class.
 */
static void
def_class (definition * defp)
{
  token tok;
  declaration dec;
  decl_list *decls;
  decl_list **tailp;

  defp->def_kind = DEF_CLASS;

  scan (TOK_IDENT, &tok);
  defp->def_name = tok.str;
  scan (TOK_LBRACE, &tok);
  tailp = &defp->def.ct.decls;
  do
    {
      get_declaration (&dec, DEF_CLASS);
      decls = ALLOC (decl_list);
      decls->decl = dec;
      *tailp = decls;
      tailp = &decls->next;
      scan (TOK_SEMICOLON, &tok);
      peek (&tok);
    }
  while (tok.kind != TOK_RBRACE);
  get_token (&tok);
  *tailp = NULL;
}
