#  Filename: fermi-init.sh
# Description: Bourne-shell flavor initialization for all FERMI software.
#              Runs fermi-setup to generate a sh script tailored
#              specifically to this user and FERMI software
#              installation, then source that.
# Author/Date: James Peachey, HEASARC/GSFC/NASA, May 3, 1999
# Modified for HEADAS December 2001
# Adapted for FERMI September 2008
#
#if [ "x$HEADAS" = x ]; then
#  echo "fermi-init.sh: WARNING -- set HEADAS and source headas-init.sh before sourcing fermi-init.sh!"			
#  echo "Do you wish to proceed with sourcing?"
	#	select yn in "Yes" "No"; do
		#		case $yn in
			#			Yes ) 
				#				if [ "x$FERMI_DIR" = x ]; then 
    				#    				echo "fermi-init.sh: ERROR -- set FERMI_DIR before sourcing fermi-init.sh"
				#				elif [ -x "$FERMI_DIR/BUILD_DIR/fermi-setup" ]; then
    				#    				export FERMI_INST_DIR=${FERMI_DIR}
    				#    				fermi_init=`$FERMI_INST_DIR/BUILD_DIR/fermi-setup sh`
    				#    				if [ $? -eq 0 -a "x$fermi_init" != x ]; then
					#					if [ -f "$fermi_init" ]; then
	    				#	    				. $fermi_init
					#					fi
					#					rm -f $fermi_init
    				#    				fi
    				#    				unset fermi_init
				#				else
    				#    				echo "fermi-init.sh: ERROR -- cannot execute $FERMI_DIR/BUILD_DIR/fermi-setup"
				#				fi
				#				break;;
			#			No )
				#				break;;
		#		esac
#done
#else
  if [ "x$FERMI_DIR" = x ]; then 
     echo "fermi-init.sh: ERROR -- set FERMI_DIR before sourcing fermi-init.sh"
  elif [ -x "$FERMI_DIR/BUILD_DIR/fermi-setup" ]; then
      export FERMI_INST_DIR=${FERMI_DIR}
      fermi_init=`$FERMI_INST_DIR/BUILD_DIR/fermi-setup sh`
      if [ $? -eq 0 -a "x$fermi_init" != x ]; then
	  if [ -f "$fermi_init" ]; then
	      . $fermi_init
	  fi
	  rm -f $fermi_init
      fi
      unset fermi_init
  else
      echo "fermi-init.sh: ERROR -- cannot execute $FERMI_DIR/BUILD_DIR/fermi-setup"
fi
#fi