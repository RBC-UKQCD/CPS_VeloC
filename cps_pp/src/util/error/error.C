#include<config.h>
#include<string>
CPS_START_NAMESPACE
/*!\file 
  \brief   Definition of Error class methods.
*/

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <stdarg.h>
#include <execinfo.h>
#include <stdlib.h>  //exit
#include <util/error.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Constructor
/*!
   The exit codes and standard messages are initialised here.
*/
//------------------------------------------------------------------
Error ERR;

static inline void ErrorPrint(FILE *fp){
#if 0
   const int BT_BUF_SIZE=20;
   void *buffer[BT_BUF_SIZE];
   char **strings;
   int npts=backtrace(buffer,BT_BUF_SIZE);
   strings = backtrace_symbols(buffer,npts);
   for(int j=0;j<npts;j++)fprintf(fp,"Node %d: backtrace %d: %s\n",j,strings[j]);
#endif
}

static inline void Exit(int status){
#if 1
   Float  *tmp = (Float *)NULL;
   *tmp = 1.;
#endif
   exit(status);
}

Error::Error() {

    error_class_name = "Error";
    
// Error message format strings

    error_string[pointer] = "Error in %s::%s :\n\tpointer %s is not initialized.\n";
    error_string[file_r] = "Error in %s::%s :\n\tcan not open file %s to read.\n";
    error_string[file_w] = "Error in %s::%s :\n\tcan not open file %s to write.\n";
    error_string[file_a] = "Error in %s::%s :\n\tcan not open file %s to append.\n";
    error_string[not_implemented] = "Error in %s::%s :\n\tnot implemented.\n";
    error_string[hardware] = "Hardware error in %s::%s :\n\t";
    error_string[general] = "Error in %s::%s : ";
 

    // Error exit values
    
    exit_value[pointer] = -1;
    exit_value[file_r] = -2;
    exit_value[file_w] = -3;
    exit_value[file_a] = -4;
    exit_value[not_implemented] = -5;
    exit_value[hardware] = -6;
    exit_value[general] = -7;

// File to which error messages are written
    
    error_file_name = "phys.error";

}



//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Error::~Error() {}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> : </tt>\n
  <tt>pointer <i>pointer name</i> is not initialized.</tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -1.

  \param class_name The name of a class whose method is being entered.
  \param func_name The name of a function or method.
  \param ptr_name The variable name of the pointer.
*/
//------------------------------------------------------------------
void Error::Pointer(const char *class_name, const char *func_name, 
		    const char *ptr_name)
{
  FILE *fp;
  error_func_name = "Pointer";

  printf(error_string[pointer], class_name, func_name, ptr_name); 

  if( (fp = Fopen(ADD_ID,error_file_name, "w")) == NULL ) { 
    printf(error_string[file_w], error_class_name, error_func_name, error_file_name); 
    Exit(exit_value[file_w]); 
  } 
  Fprintf(fp, error_string[pointer], class_name, func_name, ptr_name); 
  ErrorPrint(fp);
  Fclose(fp); 
  
  Exit(exit_value[pointer]); 
  
}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> : </tt>\n
  <tt>can not open file <i>file name</i> to read.</tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -2.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param file_name The name of the file.
*/
//------------------------------------------------------------------
void Error::FileR(const char *class_name, const char *func_name,
		  const char *file_name)   // file name
{
  FILE *fp;
  error_func_name = "FileR";

  printf(error_string[file_r], class_name, func_name, file_name);

  if( (fp = Fopen(ADD_ID,error_file_name, "w")) == NULL ) {
    printf(error_string[file_w], error_class_name, error_func_name, error_file_name);
    Exit(exit_value[file_w]);
  }
  Fprintf(fp, error_string[file_r], class_name, func_name, file_name);
  ErrorPrint(fp);
  Fclose(fp);

  Exit(exit_value[file_r]);
}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> : </tt>\n
  <tt>can not open file <i>file name</i> to write.</tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -3.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param file_name The name of the file.
*/
//------------------------------------------------------------------
void Error::FileW(const char *class_name, const char *func_name,
		  const char *file_name)   // file name
{
  FILE *fp;
  error_func_name = "FileW";

  printf(error_string[file_w],class_name, func_name, file_name);

  if( (fp = Fopen(ADD_ID,error_file_name, "w")) == NULL ) {
    printf(error_string[file_w], error_class_name, error_func_name, error_file_name);
    Exit(exit_value[file_w]);
  }
  Fprintf(fp, error_string[file_w],class_name, func_name, file_name);
  ErrorPrint(fp);
  Fclose(fp);

  Exit(exit_value[file_w]);
}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> :</tt>\n
  <tt>can not open file <i>file name</i> to append.</tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -4.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param file_name The name of the file.
*/
//------------------------------------------------------------------
void Error::FileA(const char *class_name, const char *func_name, 
		  const char *file_name)   // file name
{
  FILE *fp;
  error_func_name = "FileA";

  printf(error_string[file_a], class_name, func_name, file_name);

  if( (fp = Fopen(ADD_ID,error_file_name, "w")) == NULL ) {
    printf(error_string[file_w], error_class_name, error_func_name, error_file_name);
    Exit(exit_value[file_w]);
  }
  Fprintf(fp, error_string[file_a], class_name, func_name, file_name);
  ErrorPrint(fp);
  Fclose(fp);

  Exit(exit_value[file_a]);
}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> :</tt>\n
  <tt>not implemented.</tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -5.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
*/
//------------------------------------------------------------------
void Error::NotImplemented(const char *class_name, const char *func_name) 
{
  FILE *fp;
  error_func_name = "NotImplemented";

  printf(error_string[not_implemented], class_name, func_name);

  if( (fp = Fopen(ADD_ID,error_file_name, "w")) == NULL ) {
    printf(error_string[file_w], error_class_name, error_func_name, error_file_name);
    Exit(exit_value[file_w]);
  }
  Fprintf(fp, error_string[not_implemented], class_name, func_name);
  ErrorPrint(fp);
  Fclose(fp);

  Exit(exit_value[not_implemented]);
}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> :</tt>\n
  <tt>not implemented.</tt>\n
  <tt><i>message</i></tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -5.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.

  \todo Why not use vprintf rather than vsprintf followed by printf?      
*/
//------------------------------------------------------------------
void Error::NotImplemented(const char *class_name, const char *func_name,
			   const char *format,  // format of message
			   ...)                 // argument list of message
{
  FILE *fp;
  error_func_name = "NotImplemented";
  va_list args;
  va_start(args, format);

  printf(error_string[not_implemented], class_name, func_name);
  vprintf(format, args);

  if( (fp = Fopen(ADD_ID,error_file_name, "w")) == NULL ) {
    printf(error_string[file_w], error_class_name, error_func_name, error_file_name);
    Exit(exit_value[file_w]);
  }
  Fprintf(fp, error_string[not_implemented], class_name, func_name);
  Vfprintf(fp, format, args);
  ErrorPrint(fp);
  Fclose(fp);

  Exit(exit_value[not_implemented]);
}
        


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Hardware error in <i>class name</i>::<i>function name</i> :</tt>\n
  <tt><i>message</i></tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -6.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.
*/
//------------------------------------------------------------------
void Error::Hardware(const char *class_name, const char *func_name,
		    const char *format,  // format of message
		    ...)                 // argument list of message
{
  FILE *fp;
  error_func_name = "Hardware";
  va_list args;
  va_start(args, format);

  printf(error_string[hardware], class_name, func_name);
  vprintf(format, args);

  if( (fp = Fopen(ADD_ID,error_file_name, "w")) == NULL ) {
    printf(error_string[file_w], error_class_name, error_func_name, error_file_name);
    Exit(exit_value[file_w]);
  }
  Fprintf(fp, error_string[hardware], class_name, func_name);
  Vfprintf(fp, format, args);
  ErrorPrint(fp);
  Fclose(fp);

  Exit(exit_value[hardware]);
}

        

//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> :</tt>\n
  <tt><i>message</i></tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -7.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.
*/
//------------------------------------------------------------------
void Error::General(const char *class_name, const char *func_name,
		    const char *format,  // format of message
		    ...)                 // argument list of message
{
  FILE *fp;
  error_func_name = "General";
  va_list args;
  va_start(args, format);

  printf(error_string[general], class_name, func_name);
  vprintf(format, args);

  if( (fp = Fopen(ADD_ID,error_file_name, "w")) == NULL ) {
    printf(error_string[file_w], error_class_name, error_func_name, error_file_name);
    Exit(exit_value[file_w]);
  }
  Fprintf(fp, error_string[general], class_name, func_name);
  Vfprintf(fp, format, args);
  ErrorPrint(fp);
  Fclose(fp);

  Exit(exit_value[general]);
}

void Error::General(std::string &cname, std::string &fname,
		    const char *format,  // format of message
		    ...)                 // argument list of message
{
  FILE *fp;
  error_func_name = "General";
  va_list args;
  va_start(args, format);

  if( (fp = Fopen(ADD_ID,error_file_name, "w")) == NULL ) {
    printf(error_string[file_w], error_class_name, error_func_name, error_file_name);
    Exit(exit_value[file_w]);
  }
  Fprintf(fp, error_string[general], cname.c_str(), fname.c_str());
  Vfprintf(fp, format, args);
  ErrorPrint(fp);
  Fclose(fp);

  Exit(exit_value[general]);
}

CPS_END_NAMESPACE
