/**
 * @file   catalog.h
 * @brief  Declaration for the Catalog class.
 * Four symbolic constants are defined (used for C array size).
 *
 * @author A. Sauvageon
 *
 * $Header $
 */

#ifndef catalogAccess_cat_h
#define catalogAccess_cat_h

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/Header.h"
#include "catalogAccess/quantity.h"
// can compile without first three
//#include <cctype>      //for toupper, tolower
#include <algorithm>   //for transform
//#include <cstdio>      //for sprintf
#include <fstream>    //for ifstream, ofstream
//#include <dirent.h>   //for DIR type 
#include <iomanip>    //for setprecision, _Ios_Fmtflags, ...
#include <stdexcept>  //for std::runtime_error

#define MAX_CAT 11              // number of known catalogs
#define MAX_GEN  6              // number of generic quantities
#define MAX_URL  9              // number of known VizieR web address
#define MAX_LINE 1025           // (maximum number of char)+1 with getline
namespace catalogAccess {

/**
 * @class   Catalog
 *
 * @brief  Provide methods to define a catalog and access its data.
 * Only inline methods (default constructor, the destructor and 6 private
 * checking methods) are implemented here.
 *
 * @author A. Sauvageon
 *
 * $Header $
 */

class Catalog {

public:

  Catalog();                    // Default constructor
  ~Catalog();                   // Destructor needed to free memory
  Catalog(const Catalog & );    // Copy constructor needed


  // Methods giving general information
  //-----------------------------------

  static void getCatList(std::vector<std::string> *names,
                         const bool isCode=true);
      // return a list of all supported catalog names
  static void getWebList(std::vector<std::string> *names,
                         const bool isCode=true);
      // return a list of all supported web site names

  int getMaxNumRows(long *nrows, const std::string &fileName,
                    const std::string ext="1");
  int getMaxNumRowsWeb(long *nrows, const std::string catName,
                       const std::string urlCode="cds"); 
      // to know the maximal number of rows in the file or in the CDS catalog 


  // Methods for importing, saving, loading
  //---------------------------------------

  long getRAMsize(const long numRows=44000, const bool writeLog=false);
      // return the number of RAM bytes needed for numRows catalog rows
      // if writeLog is true (as called by import below) write info in Log

  int importDescription(const std::string &fileName, const std::string ext="1");
  int importDescriptionWeb(const std::string catName,
                           const std::string urlCode="cds",
                           const std::string &fileName="");
      // read only the catalog description (i.e. fill the
      // m_quantities vector and m_code, etc.)
      // the method returns the number of quantities,
      // -1 if importDescription already done, -3 if catName is unknown
      // other negative number for loading error


  int import(const std::string &fileName, const long maxRows=0,
             const std::string ext="1");
  int importWeb(const std::string catName, const std::string urlCode="cds",
                const long maxRow=44000, const std::string &fileName="");
      // import method for loading an entire catalog without selection
      // * recognizes catalog name and identifies the files to load
      // * translates the generic quantities
      // * loads non-generic quantities
      // the method returns the number of loaded rows,
      // -1 if successful import already done, -3 if catName is unknown,
      // other negative number for loading error

  int importSelected(std::string &filter);
      // if a catalog description was already loaded, this method does
      // the same as import(). However, it applies selection criteria
      // such that quantities which are not passing the criteria are not 
      // loaded;
      // the method returns the number of loaded rows,
      // -1 if successful import already done, -2 if importDescription not done
      // other negative number for loading error      

  int saveText(const std::string &fileName, bool clobber=false);
      // save the catalog information presently in memory to a text file
      // the method returns 1 if successful, negative number otherwise

  int saveSelectedText(const std::string &fileName, bool clobber=false);
      // like saveText(), however, storing only the selected rows
      // the method returns 1 if successful, negative number otherwise

  int saveFits(const std::string &fileName, const std::string &extName,
               bool clobber=false, bool append=false);
      // save the catalog information presently in memory to a FITS file
      // the method returns 1 if successful, negative number otherwise

  int saveSelectedFits(const std::string &fileName, const std::string &extName,
                       bool clobber=false, bool append=false);
      // like saveFits(), however, storing only the selected rows
      // the method returns 1 if successful, negative number otherwise


  // Methods for accessing data
  //---------------------------

  // accessing the Catalog definition

  void getCatalogTitles(std::vector<std::string> *titles);
      // get the 6 definition strings if vector size is enough
      // otherwise get only m_code (size 0 or 1) or m_URL, etc.
/*  void getQuantityIter(std::vector<Quantity>::const_iterator *iter);*/
      // get an iterator on the quantities
  int getQuantityDescription(std::vector<Quantity> *myQuantities);
      // get a copy of the quantity vector
  int getQuantityNames(std::vector<std::string> *names);  // get only the names
  int getQuantityUnits(std::vector<std::string> *units);  // get only the units
  int getQuantityUCDs (std::vector<std::string> *ucds);   // get only the UCDs
  int getQuantityTypes(std::vector<Quantity::QuantityType> *types);
      // get only the types

  int getStatErrorName(const std::string name, std::string *statErrName);
      // get the member "statErrName" of given quantity "name",
      // identifying the statistical error quantity
  int getSysErrorName(const std::string name, std::string *sysErrName);
      // get the member "sysErrName" of given quantity "name",
      // identifying the systematic error quantity
/*  int getVecQNames(const std::string name,
                   std::vector<std::string> *vecNames);*/
      // get the list of the names of the quantities which
      // form the elements of the given VECTOR "name"


  // accessing all the Catalog contents in memory
  //       (ignoring in-memory selection criteria)


  void deleteContent();
      // erase m_strings, m_numericals but keep catalog definition

  void getNumRows(long *nrows);  // get the number of rows in the catalog
  int getSValue(const std::string name, const long row, std::string *stringVal);
      // get the value of the given string quantity "name"
      // in the given catalog row
  int getNValue(const std::string name, const long row, double *realVal);
      // get the value of the given numerical quantity "name"
      // in the given catalog row
/*  int getVecValues(const std::string name, const  int row,
                   std::vector<double> *vecVal);*/
      // get the values of the given vector quantity "name"
      // in the given catalog row
  int getStatError(const std::string name, const long row, double *realValStat);
      // get the value of the statistical error associated with quantity "name"
      // in the given catalog row;
      // value is negative if unavailable
  int getSysError(const std::string name, const long row, double *realValSys);
      // get the value of the sytematic error associated with quantity "name"
      // in the given catalog row;
      // value is negative if unavailable
/*  int getVecStatErrors(const std::string name, const long row,
                       std::vector<double> *vecValStat);*/
      // get the statistical errors associated with vector quantity "name"
      // in the given catalog row
/*  int getVecSysErrors(const std::string name, const long row,
                      std::vector<double> *vecValSys);*/
      // get the systematic errors associated with vector quantity "name"
      // in the given catalog row

  // possible values or range of a given quantity "name"

  int getSValues(const std::string name, std::vector<std::string> *values);
      // for numerical quantity: empty
      // for string quantity: the list of different values assumed
  int minVal(const std::string name, double *realVal);
      // for string quantity: undefined; default == NO_SEL_CUT
      // for numerical quantity: the minimum value in the catalog
  int maxVal(const std::string name, double *realVal);
      // for string quantity: undefined; default == NO_SEL_CUT
      // for numerical quantity: the maximum value in the catalog


  // accessing the Catalog contents with in-memory selection criteria applied
  // the row index relates to the selected rows, i.e. is continuous!

  void getNumSelRows(long *nrows);
      // get the number of selected rows in the catalog
  int getSelSValue(const std::string name, const long srow,
                   std::string *stringVal);
  int getSelNValue(const std::string name, const long srow, double *realVal);
/*  int getSelVecValues(const std::string name, const long srow,
                      std::vector<double> *vecVal);*/
  int getSelStatError(const std::string name, const long srow,
                      double *realValStat);
  int getSelSysError(const std::string name, const long srow,
                     double *realValSys);
/*  int getSelVecStatErrors(const std::string name, const long srow,
                          std::vector<double> *vecValStat);
  int getSelVecSysErrors(const std::string name, const long srow,
                       std::vector<double> *vecValSys);*/

  // possible values or range of a given quantity "name" in the selected rows

  int getSelSValues(const std::string name, std::vector<std::string> *values);
  int minSelVal(const std::string name, double *realVal);
  int maxSelVal(const std::string name, double *realVal);


  // quick access to generic quantities with or without in-memory selection;

  int ra_deg(const long row, double *realVal, const bool inSelection=false);
      // access to generic quantity for RA  (degrees)
  int dec_deg(const long row, double *realVal, const bool inSelection=false);
      // access to generic quantity for DEC (degrees)
  int posError_deg(const long row, double *realVal,
                   const bool inSelection=false);
      // access to generic quantity for position uncertainty (degrees)

  // 6 methods returning the name of generic quantities (to use any methods)

  std::string getNameRA();
  std::string getNameDEC();
  std::string getNamePosErr();
  std::string getNameL();
  std::string getNameB();
  std::string getNameObjName();


  // Methods for selecting data
  //---------------------------

  // All setting or unsetting of cuts immediately takes effect,
  // i.e. the value of the vector m_rowIsSelected is recalculated.
  // Except for the 2 following methods with no impact on m_rowIsSelected:

  int selectQuantity(const std::string name, const bool toBeLoaded);
      // select or unselect a given quantity before loading,
      // importSelected() will actually load selected quantities
  int selectAllQuantities(const bool toBeLoaded);
      // idem on all quantities at once

  // comment: all 'unset' methods have only an effect if the cuts were applied
  // after loading the data AND the eraseNonSelected method was not called

  int unsetCuts(const std::string name);
      // unset all selection criteria relating to quantity "name"
  int unsetCuts();
      // unset all cuts on all quantities except the selection ellipse;
      // this also deletes the selection string

  int setCriteriaORed(const bool bitOR=true);
      // if true, criteria between quantities are ORed instead of ANDed

  int setSelEllipse(const double centRA_deg, const double centDEC_deg,
                    const double majAxis_deg, const double minAxis_deg,
                    const double rot_deg=0.);
      // set and apply an elliptical selection region
      // (box cuts of constant size CANNOT be achieved)
  int unsetSelEllipse();        // remove the effects of the ellipse selection

  int setLowerCut(const std::string name, double cutVal);
      // set and apply a cut on quantity "name" (all values >= cutVal pass)
      // double is not const because it can be locally modified to NO_SEL_CUT
  int setUpperCut(const std::string name, double cutVal);
      // set and apply a cut on quantity "name" (all values <= cutVal pass)
      // double is not const because it can be locally modified to NO_SEL_CUT
/*  int setLowerVecCuts(const std::string name,
                      const std::vector<double> &cutValues);*/
      // set and apply a cut on quantities in VECTOR type quantity "name"
      // such that all values >= cutValues[i] pass
/*  int setUpperVecCuts(const std::string name,
                      const std::vector<double> &cutValues);*/
      // set and apply a cut on quantities in VECTOR type quantity "name"
      // such that all values <= cutValues[i] pass
  int excludeS(const std::string name, const std::vector<std::string> &slist,
               const bool exact=false);
      // exclude all rows which have string value in the given list
  int useOnlyS(const std::string name, const std::vector<std::string> &slist,
               const bool exact=false);
      // only include rows which have string value in the given list
  int excludeN(const std::string name, const std::vector<double> &listVal);
      // include rows which value is NOT around one numerical in list
      // AND value is inside the cut interval
  int useOnlyN(const std::string name, const std::vector<double> &listVal);
      // only include rows which value is around one numerical in list
      // AND is inside the cut interval
  int includeN(const std::string name, const std::vector<double> &listVal);
      // only include rows which value is around one numerical in list
      // OR is inside the cut interval
  int setRejectNaN(const std::string name, const bool rejectNaN);
      // set quantity member m_rejectNaN and apply
      // (if different from previous and quantity is selected)  
  int setMatchPercent(const std::string name, double percent);
  int setMatchEpsilon(const std::string name, const unsigned long step);
      // set quantity member m_precision and apply
      // (if quantity has a non empty selection list)   
  int checkValWithinCut(const std::string name, const double value,
                        bool *cutLow, bool *cutUp, int *checkResult);
      // check if value (that could be used for list selection) is incompatible
      // with the cut of given quantity
      // Is given value within interval cut (boolean true if cut exist) ?
      // 0 means within interval (or no cut, both boolean false)
      // negative if below LowerCut, positive if above UpperCut;
      // -2 or 2 if value and precision sphere outside cut (useOnlyN can't
      //         select anything, value useless for excludeS)
      // -1 or 1 if value outside cut but precision sphere cross cut

  int eraseNonSelected();
      // erase all non-selected rows from memory
  int eraseSelected();
      // erase all selected rows from memory


  // general cut described by string to be parsed

/*  int setCutString(const std::string stringVal);*/
      // set the selection string m_selection
      // Syntax:
      // a) quantities are described by giving their name
      //    in square brackets, e.g. [f6cm]
      // b) otherwise FORTRAN syntax is used, e.g.
      //   "[f6cm].geq.[f12cm].and.([hr1]-[hr2]).lt.2.3"
      //
      // returns 1 if expression is correct, negative number otherwise

/*  void getCutString(std::string *stringVal);*/
      // get a copy of m_selection
    
/*  int selStringTrue(const long row, bool *isSelected);*/
      // return true if the condition described by m_selection is true
      // for the given row, otherwise false


  // Methods for sorting
  //--------------------

  // will use sorting methods present in the standard template library

//  void sortAscend(const std::string quantityName);
      // reassign the row numbers by sorting by the given
      // quantity in ascending order

//  void sortDecend(const std::string quantityName);
      // reassign the row numbers by sorting by the given
      // quantity in decending order

/**********************************************************************/
private:

  std::string m_code;
  std::string m_URL;
  std::string m_catName;
  std::string m_catRef;
  std::string m_tableName;
  std::string m_tableRef;
  std::string m_filename; // input file (with extension for fits)
  long        m_filePos;  // position where data start for importSelected()
  double m_posErrSys;
  double m_posErrFactor;

  std::vector<Quantity> m_quantities;      // the definition of the catalog
  std::vector<bool>     m_loadQuantity;    // which quantities to load ?

  std::vector<std::vector<std::string> > m_strings;
      // stores all string contents of the catalog;
      // e.g. m_strings[3] gives you the vector containing
      // the values of quantity 3 from the catalog;
      // m_strings[3][25] gives you the value of quantity 3 for catalog entry 25

  // comment: the number string quantities in the catalog == m_strings.size()

  std::vector<std::vector<double> > m_numericals;
      // stores all numerical contents of the catalog

  // comment: the number of numerical quantities == m_numericals.size()

  long m_numOriRows;
      // maximal number of rows in CDS catalog
  long m_numRows;
      // number of catalog rows loaded into memory
      // ( == m_strings[0].size() if m_strings.size() != 0 which should always
      //   be true due to the existence of generic quantities )
      // ( also == m_numericals[0].size if m_numericals.size != 0 which should
      //   always be true due to the existence of generic quantities )

  std::vector<std::vector<unsigned long> > m_rowIsSelected; 
      // false by default, i.e. all bits to 0;
      // each vector column has m_numRows elements; for a given row:
      // first bit to 1 if all selection criteria are met
      // each following bits for one quantity criteria

  std::string m_selection;      // to contain a general cut which is parsed
                                // by the method setCutString()
  bool m_criteriaORed;
      // if true: OR instead of AND between bits of m_rowIsSelected

  bool m_selRegion;
      // true if an elliptical region is to be selected;
      // default == false
  double m_selEllipseCentRA_deg;
      // the RA of the center of the ellipse (degrees)
  double m_selEllipseCentDEC_deg;
      // the DEC of the center of the ellipse (degrees)
  double m_selEllipseMinAxis_deg;
      // the size of the minor axis (degrees)
  double m_selEllipseMajAxis_deg;
      // the size of the major axis (degrees)
  double m_selEllipseRot_deg;
      // the rotation angle, i.e. the angle between major axis and the
      // celestial equator (degrees); default == 0

  // following four data members needed for efficient selection
  long m_numSelRows;            // for quick test: 0 = nothing selected
  int m_indexErr;               // index for position error in Quantity vector
  int m_indexRA;                // index for RA  in Quantity vector
  int m_indexDEC;               // index for DEC in Quantity vector
  std::vector<double> m_selEllipse;
      // the sinus and cosinus of the 2 spherical angles + ellipse size

  void deleteDescription();
      // erase the changes made by "importDescription"
      // must call first "deleteContent" if "import" or "importSelected" done
  void deleteQuantities();
      // erase elements in m_quantities according to m_loadQuantity,
      // deleteContent() MUST be done before

  bool existCriteria(std::vector<bool> *quantSel);
      // return true if at least one criteria or selection region exist

  bool checkRegion(const long row, const int nRA, const int nDEC);
      // check if given row is inside the elliptical region,
      // nRA and nDEC are the position inside m_numericals.
  bool checkNUM(const double r, const int index, const bool miss,
                const bool reject, const double precis);
      // check if value pass criteria for given quantity index (cut AND list)
  bool checkNUMor(const double r, const int index,
                  const bool reject, const double precis);
      // check if value pass criteria for given quantity index (cut OR list)
  bool rowSelect(const long row, const std::vector<bool> &quantSel);
      // compute the global row selection from bits in m_rowIsSelected
  void unsetCuts(const int index);
      // unset cut on quantity found by its existing index

  int doSelS(const std::string name, const int index, const int code,
             const std::vector<std::string> &list, const bool exact);
      // select rows depending on the given string list
  int doSelN(const std::string name, const int index, const int code,
             const std::vector<double> &listVal);
      // select rows depending on the given numerical list

  void setPosErrFactor(const int index);
      // called by setGeneric() to set m_posErrFactor
  void setGeneric(const int whichCat);
      // called just after "import..." to set m_isGeneric, m_index*
  void create_tables(const int nbQuantNum, const long maxRows);
      // creates a new column in m_strings, m_numericals
  void add_rows(const long maxRows);
      // creates a new row in m_strings, m_numericals
  void translate_cell(std::string mot, const int index);
      // loads one quantity at last row (m_numRows);

  int analyze_fits(const tip::Table *myDOL, const bool getDescr,
                   const std::string origin, long *maxRows);
  int analyze_head(unsigned long *tot, int *what, bool *testCR, std::fstream*);
  int analyze_body(unsigned long *tot, int *what, const bool testCR,
                   const bool getDescr, std::fstream*, long *maxRows);
      // 3 methods read file for import or importDescription (getDescr=true)
      // returns IS_OK for completion, otherwise strictly negative number

  int load(const std::string &fileName, const std::string ext,
           const bool getDescr, long *maxRows);
      // common code between import and importDescription
  int loadWeb(const std::string catName, const std::string urlCode,
              const std::string &fileName, const long maxRow);
      // common code between importWeb and importDescriptionWeb
  int loadSelectFits(const std::string &fileName, const std::string ext,
                     long *maxRows, std::string &filter);
      // code for fits file importSelected()
  int loadSelected(unsigned long *tot, std::fstream*, long *maxRows);
      // code for ASCII file importSelected()

  // create catalog header from memory to a text file
  int createText(const std::string &fileName, bool clobber,
                 const std::string origin);
  // create catalog header from memory to a FITS file
  int createFits
     (const std::string &fileName, const std::string &extName, bool clobber,
      bool append, const std::string origin, tip::Table **ptrTable);

  // inline private methods

  int checkImport(const std::string origin, const bool isDone);
      // return -1 or -2 if problem, 0 or positive number otherwise
      // if isDone is true, cannot return 0=IS_VOID
      // if isDone is true and import is done, return m_quantities.size
      // if isDone is false and import NOT done, return 0
  int checkCatName(const std::string origin, const std::string catName);
      // return -3 if catName do not exist, its index otherwise
  int checkSize_row(const std::string origin, const long row);
      // return strictly positive number (m_numRows) if row exist
  int checkQuant_name(const std::string origin, const std::string name);
      // return negative number if problem, quantity index otherwise
  int checkSel_row(const std::string origin, const long srow);
      // return strictly positive number (m_numSelRows) if selected row exist
  unsigned long bitPosition(const int index, int *k);
      // return long int to test m_quantities[index] bit in m_rowIsSelected[k]

  // constant members
  static const char *s_CatalogURL[MAX_URL];
  static const char *s_CatalogList[2*MAX_CAT];
  static const char *s_CatalogGeneric[MAX_CAT][MAX_GEN];
  static const std::string s_genericL;
  static const std::string s_genericB;

}; // end top class definition


/**********************************************************************/
/*  DEFINING inline FUNCTION MEMBERS                                  */
/**********************************************************************/

// Default constructor
inline Catalog::Catalog() {

  m_code="";
  m_URL ="";
  m_catName   ="";
  m_catRef    ="";
  m_tableName ="";
  m_tableRef  ="";
  m_filename  ="";
  m_filePos   =0;
  m_posErrSys = -1.0;
  m_posErrFactor=1.0;  // "deg" by default

  m_numRows   =IMPORT_NEED;
  m_numOriRows=0;
  m_selection ="";
  m_criteriaORed=false;
  m_selRegion =false;
  // following four data members needed for efficient selection
  m_numSelRows=0;
  m_indexErr= -1;
  m_indexRA = -1;
  m_indexDEC= -1;
  try { m_selEllipse.assign(7, 0.0); }
  catch (const std::exception &err) {
    std::string errText;
    errText=std::string("EXCEPTION on creating m_selEllipse: ")+err.what();
    printErr("Catalog constructor", errText);
    throw;
  }
  verbosity=3;
  #ifdef DEBUG_CAT
  std::cout << "!! DEBUG Catalog constructor (sizes=" <<  m_quantities.size()
            <<","<< m_strings.size() <<","<< m_numericals.size() <<","
            << m_rowIsSelected.size() <<", m_selEllipse "<< m_selEllipse.size()
            << ")" << std::endl;
  #endif
}

/**********************************************************************/
// Destructor needed to free memory
inline Catalog::~Catalog() {

  #ifdef DEBUG_CAT
  std::cout << "!! DEBUG Catalog destructor on: "
            << m_tableName << std::endl;
  #endif
  deleteContent();
  deleteDescription();
}

/**********************************************************************/
// check if import or importDescription already done
inline int Catalog::checkImport(const std::string origin, const bool isDone) {

  int quantSize=m_quantities.size();
  if (quantSize > 0) {
    if (isDone) return quantSize;
    if (m_numRows <= 0) {
      printLog(2, "deleting previous Catalog description");
      deleteDescription();
      return IS_VOID;
    }
    printWarn(origin, "call 'deleteContent' before importing again");
    return IMPORT_BIS;
  }
  else {
    if (! isDone) return IS_VOID;
    printWarn(origin, "must first use one 'import' method");
    return IMPORT_NEED;
  }
}

/**********************************************************************/
// check if catName is valid 
inline int Catalog::checkCatName(const std::string origin,
                                 const std::string catName) {

  for (int i=0; i<MAX_CAT; i++) {
    if (Catalog::s_CatalogList[2*i] == catName) return i;
  }
  std::string errText;
  errText="given Catalog name ("+catName+") do not exist";
  printWarn(origin, errText);
  return BAD_CATNAME;
}

/**********************************************************************/
// check if catalog was succesfully loaded and row exist
inline int Catalog::checkSize_row(const std::string origin, const long row) {

  if (m_numRows <= 0) {
    printWarn(origin, "catalog is empty");
    return IS_VOID;
  }
  if ((row < 0) || (row >= m_numRows)) {
    std::ostringstream sortie;
    sortie << "row must be within [ 0, " << m_numRows-1 << "]";
    printWarn(origin, sortie.str());
    return BAD_ROW;
  }
  return m_numRows;
}

/**********************************************************************/
// check if quantity name is valid
inline int Catalog::checkQuant_name(const std::string origin,
                                    const std::string name) {

  int quantSize=m_quantities.size();
  for (int i=0; i<quantSize; i++) {
    if (m_quantities[i].m_name == name) return i;
  }
  std::string errText;
  errText="given Quantity name ("+name+") do not exist";
  printWarn(origin, errText);
  return BAD_QUANT_NAME;
}

/**********************************************************************/
// check if catalog was succesfully loaded and selected row exist
inline int Catalog::checkSel_row(const std::string origin, const long srow) {

  if (m_numRows <= 0) {
    printWarn(origin, "catalog is empty");
    return IS_VOID;
  }
  if (m_numSelRows == 0) {
    printWarn(origin, "no row is selected");
    return IS_VOID;
  }
  if ((srow < 0) || (srow >= m_numSelRows)) {
    std::ostringstream sortie;
    sortie << "row must be within [ 0, " << m_numSelRows-1 << "]";
    printWarn(origin, sortie.str());
    return BAD_ROW;
  }
  return m_numSelRows;
}

/**********************************************************************/
// return long int to test m_quantities[index] bit in m_rowIsSelected[k]
inline unsigned long Catalog::bitPosition(const int index, int *k) {

  unsigned long test=1ul;
  int pos=index+1; // position of given criteria in isSelected vector
  *k=(pos+1)/(sizeof(long)*8);
  // index of m_rowIsSelected with given criteria
  int i=pos+1-(*k)*sizeof(long)*8;
  // bit pos.of given criteria in m_rowIsSelected[k]
  int j=0;
  while (j < i) {test*=2ul; j++;}

  return test;
}

} // namespace catalogAccess
#endif // catalogAccess_cat_h
