/* 
   


 */


//extern "C" {
#include <R.h>
//#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>
//}

/**
Defines:
*/

#define LOCUSLEN         5
#define ALLELELEN        4
#define NONGENOTYPECOLS  6

///names of parameters in the simulation.  Used as names in lists

#define INTEGERPARAMS    "intparam"
#define SWITCHPARAMS     "switchparam"
#define FLOATPARAMS      "floatparam"
#define DEMOPARAMS       "demography"
#define LOCIPARAMS       "loci"
#define INDPARAMS        "individuals"

#define HABNAMES         "habitats"     
#define STAGENAME 	 "stages"
#define XDIMNAME 	 "xdim"
#define YDIMNAME 	 "ydim"       
#define	LNUMNAME	 "locusnum"     
#define	ENUMNAME	 "numepochs"    
#define	CGNAME		 "currentgen"   
#define	CENAME		 "currentepoch" 
#define	FINALAGE	 "totalgens"    
#define	DNUMNAME	 "numdemos"     
#define	MAXLANDNAME	 "maxlandsize"  
#define NEXTIDNAME       "nextid"

#define TYPENAME	 "type"    
#define	PLOIDYNAME	 "ploidy"  
#define	TRANSNAME	 "trans"   
#define	RATENAME	 "rate"    
#define	ALISTNAME 	 "alleles" 

#define AINDXNAME        "aindex"  
#define	ABIRTHNAME	 "birth"   
#define	PROPNAME	 "prop"    
#define	STATENAME	 "state"   

#define LOCALDEMNM       "localdem"
#define LOCALDEMKNM      "localdemK"
#define EPOCHDEMNM       "epochs"

#define LCLSMATNM        "LocalS"
#define LCLRMATNM        "LocalR"
#define LCLMMATNM        "LocalM"

#define RNDCHSNAME       "RndChooseProb"
#define	SGENAME	         "StartGen"     
#define	EXTINCTNAME	 "Extinct"      
#define	CARRYNAME	 "Carry"        
#define	LPNAME		 "Localprob"    
#define SNAME            "S"    
#define	RNAME		 "R"            
#define	MNAME      	 "M"            

#define RANDEPOCHN       "randepoch"
#define RANDDEMON 	 "randdemo" 
#define MULTPNAME        "multp"
///KKM 5.20.05..............................................................
#define DENSDEP     "densdepdemo"
///.........................................................................
#define SELFRATENAME     "selfing"

///input/output of landscapes to/from metasim lib 
extern "C" SEXP read_landscape(SEXP fn); //DEPRECATED/NO CALLS FROM R NOW
//extern "C" SEXP convert_metasim_to_R(Landscape_statistics &L);
extern "C" SEXP write_landscape(SEXP fn, SEXP Rland); //DEPRECATED/NO CALLS FROM R NOW
//extern "C" void convert_R_to_metasim(SEXP Rland, Landscape_statistics &L);
extern "C" SEXP getListElement(SEXP list, const char *str);

///simulations
///run metasim on the landscape a certain number of times
extern "C" SEXP iterate_landscape(SEXP numit, SEXP Rland, SEXP cmpress, SEXP bypop);
///perform survival step on the landscape
extern "C" SEXP survive_landscape(SEXP Rland);
///perform reproduce step on the landscape
extern "C" SEXP reproduce_landscape(SEXP Rland);
///perform carry step on the landscape
extern "C" SEXP carry_landscape(SEXP Rland);
///perform extinct step on the landscape
extern "C" SEXP extinct_landscape(SEXP Rland);
///advance the landscape through time
extern "C" SEXP advance_landscape(SEXP Rland);


extern "C" SEXP populate_Rland(SEXP Rland, SEXP Population_sizes);

extern "C" SEXP clean_landscape(SEXP Rland);
extern "C" SEXP compress_landscape(SEXP Rland);

extern "C" SEXP num_demo_cols();
///utility functions
///convert a landscape into a format that the weir fst calculations in R can use.
extern "C" SEXP l2w(SEXP Rland, SEXP numind); //DEPRECATED NO CALL FROM R

///return the maximum loci compiled into this version of rmetasim
extern "C" SEXP num_loci_poss();


///
///this is C code to calculate relatedness
///

extern "C" SEXP relateinternal(SEXP ind, SEXP acnp);
 


extern "C" SEXP test(SEXP mat1, SEXP mat2);



static const R_CallMethodDef CallEntries[] = {
    {"advance_landscape",   (DL_FUNC) &advance_landscape,   1},
    {"carry_landscape",     (DL_FUNC) &carry_landscape,     1},
    {"clean_landscape",     (DL_FUNC) &clean_landscape,     1},
    {"compress_landscape",  (DL_FUNC) &compress_landscape,  1},
    {"extinct_landscape",   (DL_FUNC) &extinct_landscape,   1},
    {"iterate_landscape",   (DL_FUNC) &iterate_landscape,   4},
    //    {"landlambda",          (DL_FUNC) &landlambda,          1},
    {"num_demo_cols",       (DL_FUNC) &num_demo_cols,       0},
    {"num_loci_poss",       (DL_FUNC) &num_loci_poss,       0},
    {"populate_Rland",      (DL_FUNC) &populate_Rland,      2},
    {"relateinternal",      (DL_FUNC) &relateinternal,      2},
    {"reproduce_landscape", (DL_FUNC) &reproduce_landscape, 1},
    {"survive_landscape",   (DL_FUNC) &survive_landscape,   1},
    {"test",                (DL_FUNC) &test,                2},
    {NULL, NULL, 0}
};
