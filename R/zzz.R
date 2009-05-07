#*** zzz.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
#* ********************************************************************************* */

.First.lib <- function(lib, pkg)
{
#   require(coda)
   library.dynam("mixAK", pkg, lib)

   invisible()
}

