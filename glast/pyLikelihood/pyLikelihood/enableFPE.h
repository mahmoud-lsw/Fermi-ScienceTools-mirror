/**
 * @file enableFPE.h
 * @brief This header file defines a function that will allow one to
 * enable floating point exception handling on posix systems.  The
 * protecting #ifdef TRAP_FPE should be set a build time.
 *
 * @author J. Chiang
 *
 * $Header: /glast/ScienceTools/glast/pyLikelihood/pyLikelihood/enableFPE.h,v 1.1.1.2 2011/03/20 19:24:46 elwinter Exp $
 */

#ifndef _pyLikelihood_enableFPE_h
#define _pyLikelihood_enableFPE_h

namespace pyLikelihood {
   void enableFPE() {
#ifdef TRAP_FPE
      feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#else
      throw std::runtime_error("Floating point exception trapping "
                               "cannot be enabled for this build.");
#endif
   }

} // namespace pyLikelihood

#endif // _pyLikelihood_enableFPE_h
