SIMPLE   =                    T / File conforms to NOST standard
BITPIX   =                    8 / Bits per pixel
NAXIS    =                    0 / No data is associated with this header
EXTEND   =                    T / Extensions may be present
DATE     =                      / Date file was made
FILENAME =                      / Name of this file
TELESCOP =                GLAST / Name of telescope generating data
INSTRUME =                  LAT / Name of instrument generating data
DATE-OBS                        / Start Date and Time of the observation (UTC)
DATE-END                        / End Date and Time of the observation (UTC)
OBSERVER = 'Michelson'          / PI name
CREATOR  =                      / Software and version creating file
HISTORY                     LatEnergyBinDef.tpl,v 1.2 2005/04/05 21:06:39 peachey Exp
END


##################################################################################
XTENSION =             BINTABLE / Binary table extension
BITPIX   =                    8 / Bits per pixel
NAXIS    =                    2 / Required value
NAXIS1   =                    0 / Number of bytes per row
NAXIS2   =                    0 / Number of rows
PCOUNT   =                    0 / Normally 0 (no varying arrays)
GCOUNT   =                    1 / Required value
TFIELDS  =                    0 / Number of columns in table
EXTNAME  =           ENERGYBINS / Extension name
TELESCOP =                GLAST / Telescope or mission name
INSTRUME =                  LAT / Instrument name
FILTER   =                 NONE / Instrument filter in use (if any)
CHANTYPE =                   PI / Energy channel type
HDUCLASS =                 OGIP / Format confirms to OGIP standard
HDUCLAS1 =             RESPONSE / Extension contains response data
HDUCLAS2 =              EBOUNDS / Extension contains response matrix
HDUVERS  =                1.2.0 / Version number of the format
DATE-OBS                        / Start Date and Time of the observation (UTC)
DATE-END                        / End Date and Time of the observation (UTC)
CREATOR  =                      / Software and version creating file
HISTORY                   LatEnergyBinDef.tpl,v 1.2 2005/04/05 21:06:39 peachey Exp


TTYPE#  CHANNEL        / Raw channel number
TFORM#  I              / Data format of this field

TTYPE#  E_MIN          / Lower boundary of the detector channel
TFORM#  1E             / Data format of this field
TUNIT#  keV            / Unit of this field

TTYPE#  E_MAX          / Upper boundary of the detector channel
TFORM#  1E             / Data format of this field
TUNIT#  keV            / Unit of this field
