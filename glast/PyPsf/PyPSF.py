# $Id: PyPSF.py,v 1.18 2009/02/03 21:40:03 elwinter Exp $

# Python program to compute, display, and compare GLAST LAT Point
# Spread Functions (PSFs).

#******************************************************************************

# Import external modules.

# Standard modules
import os
import Tkinter

# Third-party modules
import Pmw
import pylab

# Project modules
import HelpWindow
import GtobspsfWidget
import GtpsfWidget
from PSFMath import PChiSq

#******************************************************************************

class PyPSFApp(Tkinter.Frame):

    #--------------------------------------------------------------------------

    def __init__(self, parent = None, *args, **kwargs):

        # Initialize the Frame base class.
        Tkinter.Frame.__init__(self, parent, *args, **kwargs)

        # Build the widgets for this widget.
        self._makeWidgets()

    #--------------------------------------------------------------------------

    def _makeWidgets(self):

        # Create the main menubar.
        self._createMenubar()

        # Create and grid the widget to run gtobspsf.
        self._gtobspsfWidget = GtobspsfWidget.GtobspsfWidget()
        self._gtobspsfWidget.grid(row = 0, column = 0)

        # Create and grid the widget to run gtpsf.
        self._gtpsfWidget = GtpsfWidget.GtpsfWidget()
        self._gtpsfWidget.grid(row = 1, column = 0)

        # Create and grid the comparison button.
        compareButton = Tkinter.Button(self, text = 'Compare',
                                       command = self._onCompare)
        compareButton.grid(row = 2, column = 0)

        # Create and grid the chi-square probability label.
        self._probLabel = Tkinter.Label(self, text = 'P(chi_sq) = 0.0')
        self._probLabel.grid(row=2, column = 1)

    #--------------------------------------------------------------------------

    def _createMenubar(self):
        menuBar = Tkinter.Menu(self.master)
        self.master.config(menu = menuBar)
        self._createFileMenu(menuBar)
        self._createEditMenu(menuBar)
        self._createHelpMenu(menuBar)

    #--------------------------------------------------------------------------

    def _createFileMenu(self, menuBar):
        fileMenu = Tkinter.Menu(menuBar)
        fileMenu.add_command(label = 'New', command = self._onFileNew)
        fileMenu.add_command(label = 'Open...', command = self._onFileOpen)
        fileMenu.add_command(label = 'Close', command = self._onFileClose)
        fileMenu.add_command(label = 'Save', command = self._onFileSave)
        fileMenu.add_command(label = 'Save As...',
                             command = self._onFileSaveAs)
        fileMenu.add_separator()
        fileMenu.add_command(label = 'Exit', command = self._onFileExit)
        menuBar.add_cascade(label = 'File', menu = fileMenu)

        # Disable most File menu commands for now.
        fileMenu.entryconfigure(1, state = 'disabled')
        fileMenu.entryconfigure(2, state = 'disabled')
        fileMenu.entryconfigure(3, state = 'disabled')
        fileMenu.entryconfigure(4, state = 'disabled')
        fileMenu.entryconfigure(5, state = 'disabled')

    def _onFileNew(self):
        pass

    def _onFileOpen(self):
        pass

    def _onFileClose(self):
        pass

    def _onFileSave(self):
        pass

    def _onFileSaveAs(self):
        pass

    def _onFileExit(self):

        # Exit from the application.
        self.master.destroy()

    #--------------------------------------------------------------------------

    def _createEditMenu(self, menuBar):
        editMenu = Tkinter.Menu(menuBar)
        editMenu.add_command(label = 'Cut', command = self._onEditCut)
        editMenu.add_command(label = 'Copy', command = self._onEditCopy)
        editMenu.add_command(label = 'Paste', command = self._onEditPaste)
        editMenu.add_command(label = 'Undo', command = self._onEditUndo)
        menuBar.add_cascade(label = 'Edit', menu = editMenu)

        # Disable all Edit menu commands for now.
        editMenu.entryconfigure(1, state = 'disabled')
        editMenu.entryconfigure(2, state = 'disabled')
        editMenu.entryconfigure(3, state = 'disabled')
        editMenu.entryconfigure(4, state = 'disabled')

    def _onEditCut(self):
        pass

    def _onEditCopy(self):
        pass

    def _onEditPaste(self):
        pass

    def _onEditUndo(self):
        pass
    
    #--------------------------------------------------------------------------

    def _createHelpMenu(self, menuBar):
        helpMenu = Tkinter.Menu(menuBar)
        helpMenu.add_command(label = 'Help', command = self._onHelpHelp)
        helpMenu.add_command(label = 'About', command = self._onHelpAbout)
        menuBar.add_cascade(label = 'Help', menu = helpMenu)

    def _onHelpHelp(self):
        path = os.environ['FERMI_INST_DIR'] + '/help/PyPSF.txt'
        helpWindow = HelpWindow.HelpWindow(path = path)

    def _onHelpAbout(self):

        # Create the 'about the program' dialog.
        Pmw.aboutversion('0.1')
        Pmw.aboutcopyright('Copyright NASA/GSFC 2008')
        Pmw.aboutcontact('For more information, contact ' + \
                         'Eric Winter (Eric.L.Winter@nasa.gov).')
        aboutDialog = Pmw.AboutDialog(self, applicationname = 'PyPSF')
    
    #--------------------------------------------------------------------------

    def _onCompare(self):

        # Fetch the PSF from gtpsf.
        (x_pred, y_pred) = self._gtpsfWidget.get_psf()

        # Fetch the PSF and error from gtobspsf.
        (x_obs, y_obs, dy_obs) = self._gtobspsfWidget.get_psf()

        # Compare the observed PSF to the model PSF by computing the
        # chi-square value.
        N = len(x_pred)
        chi_sq = 0.0
        for i in range(0, N):
            if y_pred[i] <= 0: continue
            chi_sq += (y_obs[i] - y_pred[i])**2 / y_pred[i]

        # The number of degrees of freedom (DOF) is one less than the
        # number of bins, since both PSFs are normalized to unit area.
        n_dof = N - 1

        # Compute the probability that the observed chi_sq should be
        # less than this value. Update the probability label.
        p = PChiSq(chi_sq, n_dof)
        self._probLabel.configure(text = 'P(chi_sq) = ' + str(p))

        # Make a simple x-y plot of the predicted and observed PSF,
        # using error basrs for the observed PSF.
        pylab.ioff()
        pylab.figure(3)
        pylab.clf()
        pylab.subplot(211)
        pylab.plot(x_pred, y_pred, label = 'Predicted')
        pylab.errorbar(x_obs, y_obs, dy_obs, label = 'Observed')
        pylab.ylabel('Normalized PSF')
        pylab.title('Comparison of predicted and observed PSF')
        pylab.grid(True)
        pylab.legend()

        # Compute and plot the residuals.
        residual = []
        for i in range(0, len(x_pred)):
            residual.append(y_obs[i] - y_pred[i])
        pylab.subplot(212)
        pylab.plot(x_pred, residual)
        pylab.xlabel('Angle (degrees)')
        pylab.ylabel('Residual = observed - predicted')
        pylab.title('PSF residuals')
        pylab.grid(True)
        pylab.show()
        
#******************************************************************************

# If run as a script, build and run the application.
if __name__ == '__main__':

    # Create the root window.
    root = Tkinter.Tk()
    Pmw.initialise(root)
    root.title('PyPSF')

    # Create and pack the application object, passing it the path to
    # the model file.
    pyPSF = PyPSFApp(root)
    pyPSF.grid(sticky = 'nsew')

    # Enter the main program loop.
    root.mainloop()
