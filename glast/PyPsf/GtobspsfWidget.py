# $Id: GtobspsfWidget.py,v 1.16 2009/03/24 13:10:10 elwinter Exp $

#******************************************************************************

# Import external modules.

# Standard modules
import math
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

class GtobspsfWidget(Tkinter.Frame):

    #--------------------------------------------------------------------------

    def __init__(self, parent = None, *args, **kwargs):
        
        # Initialize the parent class.
        Tkinter.Frame.__init__(self, parent)

        # Build the widgets for this widget.
        self._makeWidgets()

    #--------------------------------------------------------------------------

    def _makeWidgets(self):

        # Create a Balloon object for help.
        self._balloon = Pmw.Balloon(self)

        # Create a Parfile object for gtobspsf.
        parfile = Parfile('gtobspsf.par')

        # Initialize the row index for gridding the widgets.
        row = 0

        #----------------------------------------------------------------------

        # Create and grid an arbitrary string EntryField for the path
        # to the event file. Initialize the contents of the EntryField
        # with the value of the 'evfile' parameter from the parameter
        # file. Then bind the appropriate help text for the help
        # balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    label_text = 'Event file')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('evfile'))
        self._balloon.bind(entryField,
                           'Enter the path to the FT1-format event file.')
        self._evfileEntryField = entryField

        # Create and grid a button which summons an Open dialog to
        # select a new event file. Then bind the appropriate help text
        # for the help balloon.
        button = Tkinter.Button(self, text = 'Browse',
                                command = self._onBrowse)
        button.grid(row = row, column = 1)
        self._balloon.bind(button, 'Select a new FT1 event file.')

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the source
        # right ascension. Initialize the contents of the EntryField
        # with the value of the 'src_ra' parameter from the parameter
        # file. Then bind the appropriate help text for the help
        # balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    validate = { 'validator' : 'real',
                                                 'min' : 0.0, 'max' : 360.0 },
                                    label_text = 'Source RA (deg)')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('src_ra'))
        self._balloon.bind(entryField,
                           'Enter the right ascension (degrees) for ' +
                           'the center of the PSF.')
        self._src_raEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the source
        # declination. Initialize the contents of the EntryField with
        # the value of the 'src_dec' parameter from the parameter
        # file. Then bind the appropriate help text for the help
        # balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    validate = { 'validator' : 'real',
                                                 'min' : -90.0, 'max' : 90.0 },
                                    label_text = 'Source DEC (deg)')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('src_dec'))
        self._balloon.bind(entryField,
                           'Enter the declination (degrees) for ' +
                           'the center of the PSF.')
        self._src_decEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the maximum PSF
        # angle. Initialize the contents of the EntryField with the
        # value of the 'maxangle' parameter from the parameter
        # file. Then bind the appropriate help text for the help
        # balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    validate = { 'validator' : 'real',
                                                 'min' : 0.0, 'max' : 90.0 },
                                    label_text = 'Maximum angle (deg)')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('maxangle'))
        self._balloon.bind(entryField,
                           'Enter the maximum angle (degrees) between ' +
                           'the source and event to consider for the PSF.')
        self._maxangleEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a integer EntryField for the number of angle
        # bins for the PSF. Initialize the contents of the EntryField
        # with the value of the 'nangles' parameter from the parameter
        # file. Then bind the appropriate help text for the help
        # balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    validate = { 'validator' : 'integer',
                                                 'min' : 1 },
                                    label_text = 'Number of angle bins')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('nangles'))
        self._balloon.bind(entryField,
                           'Enter the number of angle bins for the PSF.')
        self._nanglesEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the minimum
        # photon energy to use in the PSF. Initialize the contents of
        # the EntryField with the value of the 'emin' parameter from
        # the parameter file. Then bind the appropriate help text for
        # the help balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    validate = { 'validator' : 'real',
                                                 'min' : 0.0 },
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
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    validate = { 'validator' : 'real',
                                                 'min' : 0.0 },
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
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    validate = { 'validator' : 'integer',
                                                 'min' : 1 },
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

        # Create and grid a real-number EntryField for the constant
        # detector effective area. Initialize the contents of the
        # EntryField with the value of 'effarea' parameter from the
        # parameter file. Then bind the appropriate help text for the
        # help balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    validate = { 'validator' : 'real',
                                                 'min' : 1.0 },
                                    label_text =
                                    'Constant effective area (cm^2)')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('effarea'))
        self._balloon.bind(entryField,
                           'Enter a constant detector effective area (cm^2).')
        self._effareaEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the constant
        # background flux. Initialize the contents of the EntryField
        # with the value of the 'bgflux' parameter from the
        # parameter file. Then bind the appropriate help text for the
        # help balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    validate = { 'validator' : 'real',
                                                 'min' : 0.0 },
                                    label_text =
                                    'Constant background flux (ph/cm^2/s)')
        entryField.grid(row = row, column = 0, sticky = 'e')
        entryField.setvalue(parfile.GetParameterValue('bgflux'))
        self._balloon.bind(entryField,
                           'Enter a constant background flux (ph/cm^2/s)')
        self._bgfluxEntryField = entryField

        # Increment the row index for the next widget.
        row += 1

        #----------------------------------------------------------------------

        # Create and grid a real-number EntryField for the power law
        # reference energy (MeV). Initialize the contents of the
        # EntryField with the value 200 MeV. Then bind the appropriate
        # help text for the help balloon.
        entryField = Pmw.EntryField(self, labelpos = 'w',
                                    validate = { 'validator' : 'real',
                                                 'min' : 1.0 },
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

        # Create and grid a button to run the gtobspsf tool. Then bind
        # the appropriate help text for the help balloon.
        button = Tkinter.Button(frame, text = 'Run gtobspsf',
                                command = self._onRun)
        button.grid(row = 0, column = 0)
        self._balloon.bind(button,
                           'Run gtobspsf to compute the observed PSF.')

        # Create and grid a button to plot the output from the
        # gtobspsf tool.
        button = Tkinter.Button(frame, text = 'Plot PSF',
                                command = self._onPlot)
        button.grid(row = 0, column = 1)
        self._balloon.bind(button, 'Plot the PSF from gtobspsf.')

    #--------------------------------------------------------------------------

    # Use the Open dialog to select a new FT1 event file. Set the
    # EntryField for the event file to this path.

    def _onBrowse(self):

        # Present the file/open dialog.
        path = Open().show()

        # If no file was selected, return.
        if path == '':
            return

        # Save the current path in the event file EntryField.
        self._evfileEntryField.setvalue(path)

    #--------------------------------------------------------------------------

    # Run the gtobspsf tool using the current values of the parameters.

    def _onRun(self):

        # Assemble the command string to run gtobspsf.
        cmd = 'gtobspsf'
        cmd += ' evfile=' + self._evfileEntryField.getvalue()
        cmd += ' effarea=' + self._effareaEntryField.getvalue()
        cmd += ' bgflux=' + self._bgfluxEntryField.getvalue()
        cmd += ' outfile=gtobspsf_out.fits'
        cmd += ' src_ra=' + self._src_raEntryField.getvalue()
        cmd += ' src_dec=' + self._src_decEntryField.getvalue()
        cmd += ' maxangle=' + self._maxangleEntryField.getvalue()
        cmd += ' nangles=' + self._nanglesEntryField.getvalue()
        cmd += ' emin=' + self._eminEntryField.getvalue()
        cmd += ' emax=' + self._emaxEntryField.getvalue()
        cmd += ' nenergies=' + self._nenergiesEntryField.getvalue()
        print cmd

        # Run the gtobspsf tool.
        os.system(cmd)
    
        # Open the result file from gtobspsf.
        gtobspsf_results = pyfits.open('gtobspsf_out.fits')

        #----------------------------------------------------------------------

        # PSF and error

        # Create an array to hold the aggregate PSF.
        nangles = int(self._nanglesEntryField.getvalue())
        PSF = []
        for i in range(0, nangles):
            PSF.append(0.0)

        # Create an array to hold the aggregate PSF error.
        PSF_err = []
        for i in range(0, nangles):
            PSF_err.append(0.0)

        # Fetch the PSF table.
        psf_data = gtobspsf_results['Psf'].data

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

        # Add up the PSF for each energy, and the associated squared
        # errors, using the computed weights.
        for record in psf_data:
            for i in range(0, nangles):
                PSF[i] += record[2][i] * W[i]
                PSF_err[i] += record[3][i]**2 * W[i]

        # Normalize the PSF to unity at the center, and normalize the
        # errors by the same amount.
        PSF0 = PSF[0]
        for i in range(0, nangles):
            PSF[i] /= PSF0
            PSF_err[i] = math.sqrt(PSF_err[i]) / PSF0

        #----------------------------------------------------------------------

        # Create an array to hold the PSF angles.
        angles = []

        # Fetch the angle table.
        theta_data = gtobspsf_results['THETA'].data

        # Copy the bin angles.
        for record in theta_data:
            angles.append(record[0]);

        #----------------------------------------------------------------------

        # Save the results.
        self._x = angles
        self._y = PSF
        self._dy = PSF_err

    #--------------------------------------------------------------------------

    def _onPlot(self):

        # Make a simple x-y plot of the PSF, with errorbars.
        pylab.ioff()
        pylab.figure(1)
        pylab.errorbar(self._x, self._y, self._dy)
        pylab.xlabel('Angle (degrees)')
        pylab.ylabel('Normalized PSF')
        pylab.title('Observed PSF from gtobspsf')
        pylab.grid(True)
        pylab.show()

    #--------------------------------------------------------------------------

    # Return a tuple of arrays containing the results of the most
    # recent run of gtobspsf.

    def get_psf(self):
        return (self._x, self._y, self._dy)

#******************************************************************************

# Self-test code.

if __name__ == '__main__':

    # Create the root window.
    root = Tkinter.Tk()
    Pmw.initialise(root)
    root.title('GtobspsfWidget')

    # Create a gtobspsf widget.
    gtobspsfWidget = GtobspsfWidget(root)
    gtobspsfWidget.grid(sticky = 'nsew')

    # Enter the main program loop.
    root.mainloop()
