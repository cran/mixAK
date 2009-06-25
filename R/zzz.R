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
   cat("\n")
   cat("### Mixture of methods including mixtures\n")
   cat("### Arnost Komarek\n\n")
   cat("### See citation(\"mixAK\") for the best way to cite\n")
   cat("### the package if you find it useful.\n\n")
   invisible()
}

