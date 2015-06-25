#orbitSim-init.sh
#Setup script for standalone version of gtorbsim
#08/27/2008

if test "x${ORBITSIMROOT}" != x; then
  # Prepend directory containing the executables to the path, weeding out redundancies.
  tmp_path=`echo ":${PATH}:" | sed "s%:${ORBITSIMROOT}/bin:%:%g" | sed 's%::*$%%'`
  PATH="${ORBITSIMROOT}/bin${tmp_path}"
  unset tmp_path
  export PATH

  # Append to, or set PFILES
  if test "x${PFILES}" != x; then
    PFILES="${PFILES}:${ORBITSIMROOT}/pfiles"
  else
    PFILES="${HOME}/pfiles;${ORBITSIMROOT}/pfiles"
    # Make sure ~/pfiles/ exists and if not make it.
    if [ ! -e "$HOME/pfiles" ]; then
	mkdir "$HOME/pfiles"
    elif [ -e "$HOME/pfiles" ] && [ ! -d "$HOME/pfiles" ]; then
	echo "Error, ~/pfiles exists but is not a directory!"
    fi
  fi
  export PFILES
else
  echo "Set ORBITSIMROOT before sourcing orbitSim-init.sh" >&2
fi
