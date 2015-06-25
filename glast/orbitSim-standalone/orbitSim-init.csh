#orbitSim-init.csh
#Setup script for standalone version of gtorbsim
#08/27/2008

if ( ${?ORBITSIMROOT} ) then
  # Prepend directory containing the executables to the path, weeding out redundancies.
  set tmp_path=`echo ":${PATH}:" | sed "s%:${ORBITSIMROOT}/bin:%:%g" | sed 's%::*$%%'`
  setenv PATH "${ORBITSIMROOT}/bin${tmp_path}"
  unset tmp_path

  # Append to, or set PFILES
  if ( ${?PFILES} ) then
      setenv PFILES "${PFILES}:${ORBITSIMROOT}/pfiles"
  else
      setenv PFILES "${HOME}/pfiles;${ORBITSIMROOT}/pfiles"
      if (! -e "${HOME}/pfiles" ) then
	  mkdir "${HOME}/pfiles"
      else if (( -e "${HOME}/pfiles" ) && ( ! -d "${HOME}/pfiles")) then
	echo "Error, ~/pfiles exists but is not a directory!"
      endifx
      endif
  endif
else
  echo "Set ORBITSIMROOT before sourcing orbitSim-init.csh"
endif
