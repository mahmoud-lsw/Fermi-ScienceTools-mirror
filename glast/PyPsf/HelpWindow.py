# $Id: HelpWindow.py,v 1.6 2009/02/03 21:40:03 elwinter Exp $

#******************************************************************************

# Import external modules.

# Standard modules
import os
import Tkinter

# Third-party modules
import Pmw

# Project modules

#******************************************************************************

class HelpWindow(Tkinter.Toplevel):

    #--------------------------------------------------------------------------

    def __init__(self, parent = None, path = None, *args, **kwargs):
        """Initialize this HelpWindow.

        Parameters:

        self: This object.

        parent: (Tkinter.Frame) Parent object for this widget.

        Return value:

        None.

        Description:
        """

        # Initialize the parent class.
        Tkinter.Toplevel.__init__(self, parent)

        # Create widgets.
        self._makeWidgets(path)

    #--------------------------------------------------------------------------

    def _makeWidgets(self, path):
        """Build child widgets.

        Parameters:

        self: This object.

        Return value:

        None.

        Description:

        Create all of the GUI widgets required by this object.
        """

        # Create the Text widget and read the help text.
        self._scrolledText = Pmw.ScrolledText(self)
        self._scrolledText.pack(expand = 'yes', fill = 'both')
        self._scrolledText.importfile(path)

        # Create the close button.
        self._closeButton = Tkinter.Button(self, text = 'Close',
                                           command = self._onClose)
        self._closeButton.pack()

    #--------------------------------------------------------------------------

    def _onClose(self):
        self.destroy()

#******************************************************************************

if __name__ == '__main__':

    # Create the root window.
    root = Tkinter.Tk()
    Pmw.initialise(root)
    root.title('HelpWindow test')

    # Create the help window.
    helpWindow = HelpWindow()

    # Enter the main program loop.
    root.mainloop()
