# $Id: GtpsfWidget.py,v 1.14 2008/11/05 21:16:26 elwinter Exp $

#******************************************************************************

# Import external modules.

# Standard modules
import os
import Tkinter
from tkFileDialog import Open

# Third-party modules
import Pmw
import pyfits
import pylab

# Project modules
from Parfile import Parfile

#******************************************************************************

class GtpsfWidget(Tkinter.Frame):

    # Define the valid response functions.
    responseFunctions = [
        'DC1',
        'DC1A',
        'DC1AB',
        'DC1AF',
        'DC1B',
        'DC1F',
        'DC2',
        'DC2BA',
        'DC2BB',
        'DC2FA',
        'DC2FB',
        'DC2_A',
        'G25',
        'G25B',
        'G25F',
        'HANDOFF',
        'HANDOFF::BACK',
        'HANDOFF::FRONT',
        'P5_v13_0_diff',
        'P5_v13_0_diff::BACK',
        'P5_v13_0_diff::FRONT',
        'P5_v13_0_source',
        'P5_v13_0_source::BACK',
        'P5_v13_0_source::FRONT',
        'P5_v13_0_trans',
        'P5_v13_0_trans::BACK',
        'P5_v13_0_trans::FRONT',
        'P6_V1_DIFFUSE',
        'P6_V1_DIFFUSE::BACK',
        'P6_V1_DIFFUSE::FRONT',
        'P6_V1_SOURCE',
        'P6_V1_SOURCE::BACK',
        'P6_V1_SOURCE::FRONT',
        'P6_V1_TRANSIENT',
        'P6_V1_TRANSIENT::BACK',
        'P6_V1_TRANSIENT::FRONT',
        'PASS4',
        'PASS4::BACK',
        'PASS4::FRONT',
        'PASS4_v2',
        'PASS4_v2::BACK',
        'PASS4_v2::FRONT',
        'PASS5_v0',
        'PASS5_v0::BACK',
        'PASS5_v0::FRONT',
        'PASS5_v0_DIFFUSE',
        'PASS5_v0_DIFFUSE::BACK',
        'PASS5_v0_DIFFUSE::FRONT',
        'PASS5_v0_TRANSIENT',
        'PASS5_v0_TRANSIENT::BACK',
        'PASS5_v0_TRANSIENT::FRONT',
        'TEST',
        'TESTB',
        'TESTF'
        ]

    #--------------------------------------------------------------------------

    def __init__(self, parent = None, *args, **kwargs):
        
        # Initialize the parent class.
        Tkinter.Frame.__init__(self, parent)

        # Build the widgets for this widget.
        self._makeWidgets()

    #--------------------------------------------------------------------------

    def _makeWidgets(self):

        # Create a Balloon object for help.
        self._balloon = Pmw.Balloon( )

        # Create a Parfile object for gtpsf.
        parfile = Parfile('gtpsf.par')

        # Initialize the row index for gridding the widgets.
        row = 0

        #----------------------------------------------------------------------

        # Create and grid an arbitrary string EntryField for the
        # livetime cube file. Initialize the contents of the
        # EntryField with the value of the 'expcube' parameter from
        # the parameter file. Then bind the appropriate help text for
        # the help balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    label_text = 'Livetime cube file')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('expcube'))
        self._balloon.bind(entryField,
                           'Enter the path to the livetime cube file ' +
                           'generated by gtltcube.')
        self._expcubeEntryField = entryField

        # Create and grid a button which summons an Open dialog to
        # select a new livetime cube file. Then bind the appropriate
        # help text for the help balloon.
        button = Tkinter.Button(self, text = 'Browse',
                                command = self._onBrowse)
        button.grid(row = row, column = 1, sticky = 'e')
        self._balloon.bind(button, 'Select a new livetime cube file ' +
                           'generated by gtltcube.')

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Find the index of the current response function in the list
        # of response functions.
        irfs = parfile.GetParameterValue('irfs')
	try:
            i_irfs = self.responseFunctions.index(irfs)
	except ValueError:
	    i_irfs = 0

        # Create and grid a OptionMenu for the instrument response
        # function.  Initialize the selection with the value of the
        # 'irfs' parameter from the parameter file. Then bind the
        # appropriate help text for the help balloon.
        optionMenu = Pmw.OptionMenu(self, labelpos = 'w',
                                    label_text = 'Response Function',
                                    items = self.responseFunctions,
                                    initialitem = i_irfs)
        optionMenu.grid(row = row, column = 0)
        self._balloon.bind(optionMenu, 'Select a response function.')
        self._irfsOptionMenu = optionMenu

        # Increment the row index for the next widget.
        row += 1
                                    
        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the source
        # right ascension. Initialize the contents of the EntryField
        # with the value of the 'ra' parameter from the parameter
        # file. Then bind the appropriate help text for the help
        # balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w', validate = 'real',
                                    label_text = 'Source RA (deg)')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('ra'))
        self._balloon.bind(entryField,
                           'Enter the right ascension (degrees) for ' +
                           'the center of the PSF.')
        self._raEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the source
        # declination. Initialize the contents of the EntryField with
        # the value of the 'dec' parameter from the parameter
        # file. Then bind the appropriate help text for the help
        # balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w', validate = 'real',
                                    label_text = 'Source DEC (deg)')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('dec'))
        self._balloon.bind(entryField,
                           'Enter the declination (degrees) for ' +
                           'the center of the PSF.')
        self._decEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the maximum PSF
        # angle. Initialize the contents of the EntryField with the
        # value of the 'thetamax' parameter from the parameter
        # file. Then bind the appropriate help text for the help
        # balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w', validate = 'real',
                                    label_text = 'Maximum angle (deg)')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('thetamax'))
        self._balloon.bind(entryField,
                           'Enter the maximum angle (degrees) between ' +
                           'the source and event to consider for the PSF.')
        self._thetamaxEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a integer EntryField for the number of angle
        # bins for the PSF. Initialize the contents of the EntryField
        # with the value of the 'ntheta' parameter from the parameter
        # file. Then bind the appropriate help text for the help
        # balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w', validate = 'integer',
                                    label_text = 'Number of angle bins')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('ntheta'))
        self._balloon.bind(entryField,
                           'Enter the number of angle bins for the PSF.')
        self._nthetaEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the minimum
        # photon energy to use in the PSF. Initialize the contents of
        # the EntryField with the value of the 'emin' parameter from
        # the parameter file. Then bind the appropriate help text for
        # the help balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w', validate = 'real',
                                    label_text = 'Minimum energy (MeV)')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('emin'))
        self._balloon.bind(entryField,
                           'Enter the minimum energy (MeV) of photons to ' +
                           'consider for the PSF.')
        self._eminEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the maximum
        # photon energy to use in the PSF. Initialize the contents of
        # the EntryField with the value of the 'emax' parameter from
        # the parameter file. Then bind the appropriate help text for
        # the help balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w', validate = 'real',
                                    label_text = 'Maximum energy (MeV)')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('emax'))
        self._balloon.bind(entryField,
                           'Enter the maximum energy (MeV) of photons to ' +
                           'consider for the PSF.')
        self._emaxEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a integer EntryField for the number of
        # logarithmic energy bins for the PSF. Initialize the contents
        # of the EntryField with the value of the 'nenergies'
        # parameter from the parameter file. Then bind the appropriate
        # help text for the help balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w', validate = 'integer',
                                    label_text = 'Number of log energy bins')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('nenergies'))
        self._balloon.bind(entryField,
                           'Enter the number of logarithmic energy bins for ' +
                           'the PSF.')
        self._nenergiesEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the power law
        # reference energy (MeV). Initialize the contents of the
        # EntryField with the value 200 MeV. Then bind the appropriate
        # help text for the help balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w', validate = 'real',
                                    label_text = 'Power law E0 (MeV)')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(200.0)
        self._balloon.bind(entryField,
                           'Enter the power law reference energy (MeV)')
        self._E0EntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the power law
        # spectral index. Initialize the contents of the EntryField
        # with the value 2.0. Then bind the appropriate help text for
        # the help balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w', validate = 'real',
                                    label_text = 'Power law spectral index')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(2.0)
        self._balloon.bind(entryField,
                           'Enter the power law spectral index')
        self._spectralIndexEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a frame to hold the control buttons.
        frame = Tkinter.Frame(self)
        frame.grid(row = row, column = 0, columnspan = 2)

        # Add a button to run the gtpsf tool. Then bind the
        # appropriate help text for the help balloon.
        button = Tkinter.Button(frame, text = 'Run gtpsf',
                                command = self._onRun)
        button.grid(row = 0, column = 0)
        self._balloon.bind(button, 'Run gtpsf to compute the predicted PSF.')

        # Create and grid a button to plot the output from the gtpsf
        # tool.
        button = Tkinter.Button(frame, text = 'Plot PSF',
                                command = self._onPlot)
        button.grid(row = 0, column = 1)
        self._balloon.bind(button, 'Plot the PSF from gtpsf.')

    #--------------------------------------------------------------------------

    # Use the Open dialog to select a new gtltcube-generated livetime
    # cube file. Set the EntryField for the livetime cube file to this
    # path.

    def _onBrowse(self):

        # Present the file/open dialog.
        path = Open().show()

        # If no file was selected, return.
        if path == '':
            return

        # Save the current path in the livetime cube file EntryField.
        self._expcubeEntryField.setvalue(path)

    #--------------------------------------------------------------------------

    def _onRun(self):

        # Assemble the command string to run gtpsf.
        cmd = 'gtpsf'
        cmd += ' expcube="' + self._expcubeEntryField.getvalue() + '"'
        cmd += ' outfile="gtpsf_out.fits"'
        cmd += ' irfs="' + self._irfsOptionMenu.getvalue() + '"'
        cmd += ' ra=' + self._raEntryField.getvalue()
        cmd += ' dec=' + self._decEntryField.getvalue()
        cmd += ' emin=' + self._eminEntryField.getvalue()
        cmd += ' emax=' + self._emaxEntryField.getvalue()
        cmd += ' nenergies=' + self._nenergiesEntryField.getvalue()
        cmd += ' thetamax=' + self._thetamaxEntryField.getvalue()
        cmd += ' ntheta=' + self._nthetaEntryField.getvalue()
        print cmd

        # Run the gtpsf tool.
        os.system(cmd)

        # Open the result file from gtpsf.
        gtpsf_results = pyfits.open('gtpsf_out.fits')

        #----------------------------------------------------------------------

        # PSF

        # Create an array to hold the aggregate PSF.
        ntheta = int(self._nthetaEntryField.getvalue())
        PSF = []
        for i in range(0, ntheta):
            PSF.append(0.0)

        # Fetch the PSF table.
        psf_data = gtpsf_results['PSF'].data

        # Fetch the energy bin values, and append the maximum energy.
        energies = []
        for record in psf_data:
            energies.append(record[0])
        emax = float(self._emaxEntryField.getvalue())
        energies.append(emax)

        # Retrieve the power law parameters.
        E0 = float(self._E0EntryField.getvalue())
        alpha = float(self._spectralIndexEntryField.getvalue())

        # Compute the weighting factors for each energy bin.
        W = []
        k1 = E0**alpha / (1 - alpha)
        for i in range(0, len(energies) - 1):
            E1 = energies[i]
            E2 = energies[i + 1]
            w = k1 * (E2**(1 - alpha) - E1**(1 - alpha))
            W.append(w)

        # Add up the PSF for each energy using the computed weights.
        for record in psf_data:
            for i in range(0, ntheta):
                PSF[i] += record[2][i] * W[i]

        # Normalize the PSF to unity at the center.
        PSF0 = PSF[0]
        for i in range(0, ntheta):
            PSF[i] /= PSF0

        #----------------------------------------------------------------------

        # Create an array to hold the PSF angles.
        angles = []

        # Fetch the angle table.
        theta_data = gtpsf_results['THETA'].data

        # Copy the bin angles.
        for record in theta_data:
            angles.append(record[0]);

        # Save the results.
        self._x = angles
        self._y = PSF

    #--------------------------------------------------------------------------

    def _onPlot(self):

        # Make a simple x-y plot of the PSF.
        pylab.ioff()
        pylab.figure(2)
        pylab.plot(self._x, self._y)
        pylab.xlabel('Angle (degrees)')
        pylab.ylabel('Normalized PSF')
        pylab.title('Predicted PSF from gtpsf')
        pylab.grid(True)
        pylab.show()

    #--------------------------------------------------------------------------

    # Return a tuple of arrays containing the results of the most
    # recent run of gtpsf.

    def get_psf(self):
        return (self._x, self._y)

#******************************************************************************

# Self-test code.

if __name__ == '__main__':

    # Create the root window.
    root = Tkinter.Tk()
    Pmw.initialise(root)
    root.title('GtpsfWidget')

    # Create a gtpsf widget.
    gtpsfWidget = GtpsfWidget(root)
    gtpsfWidget.grid(sticky = 'nsew')

    # Enter the main program loop.
    root.mainloop()