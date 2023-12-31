/*
 * Please do not edit this file.
 * It was generated using PAB's VML system.
 */

#include <alg/kinetic_smear_arg.h>
CPS_START_NAMESPACE
	 bool KineticSmearArg::Encode(char *filename,char *instance){
		 VML vmls;
		 if ( !vmls.Create(filename,VML_ENCODE)) return false;
		 if ( !Vml(&vmls,instance) ) return false;
		 vmls.Destroy(); return true;
	 }

	 bool KineticSmearArg::Decode(char *filename,char *instance){
		 VML vmls;
		 if ( !vmls.Create(filename,VML_DECODE)) return false;
		 if ( !Vml(&vmls,instance)) return false;
		 vmls.Destroy(); return true;
	 }
	 bool KineticSmearArg::Vml(VML *vmls,char *instance){
		 if(!vml_KineticSmearArg(vmls,instance,this)) return false;
	 return true;
	}


bool_t
vml_KineticSmearArg (VML *vmls, char *name,KineticSmearArg *objp)
{
	 vml_class_begin(vmls,"KineticSmearArg",name);
	 if (!vml_int (vmls, "orthog", &objp->orthog))
		 return FALSE;
	 if (!vml_Float (vmls, "single_link", &objp->single_link))
		 return FALSE;
	 if (!vml_Float (vmls, "three_link", &objp->three_link))
		 return FALSE;
	 if (!vml_Float (vmls, "five_link", &objp->five_link))
		 return FALSE;
	 if (!vml_Float (vmls, "seven_link", &objp->seven_link))
		 return FALSE;
	 if (!vml_Float (vmls, "lepage", &objp->lepage))
		 return FALSE;
	 vml_class_end(vmls,"KineticSmearArg",name);
	return TRUE;
}
CPS_END_NAMESPACE
