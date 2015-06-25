/**
 * @file   catalog_io.cxx
 * @brief  Read/write (or import/export) routines for Catalog class.
 * Methods exist to access both files: ASCII or fits (writing only BINARY table,
 * reading any table format).
 * Methods to load catalog through web are placeholders.
 *
 * @author A. Sauvageon
 *
 * $Header $
 *
 */

#include "catalogAccess/catalog.h"

namespace catalogAccess {
  using namespace tip;
/**********************************************************************/
/*  DEFINING CLASS CONSTANTS                                          */
/**********************************************************************/

// BEWARE: the two constant below are UCD1 standard
//         need to be changed for UCD1+

// the order is important and must be compatible with s_CatalogGeneric
static const char *UCD_List[MAX_GEN]={
       "ID_MAIN", "POS_EQ_RA_MAIN", "POS_EQ_DEC_MAIN",
        "",       "POS_GAL_LON",    "POS_GAL_LAT"};
// "ERROR", "POS_GAL_LON", "POS_GAL_LAT"
// ERROR do not only refer to the position error
// ERROR is for any "Uncertainty in Measurements"

static const char *UCD_Added[MAX_GEN]={
       "", "_RAJ2000", "_DEJ2000", "", "_Glon", "_Glat"};

// the 2 constants above are only defined in this file,
// need to set catalog static members for global availability

const std::string Catalog::s_genericL = UCD_List[4];
const std::string Catalog::s_genericB = UCD_List[5];
const std::string KeyCDS[2]={"CDS-CAT","CDS-NAME"};
const std::string KeyCom[2]={"Catalogue designation in CDS nomenclature",
                             "Table used in a VizieR Query"};
const std::string Key_UCD="TBUCD";
const char SepNull = '!';

/**********************************************************************/
/*  METHODS for IMPORTING, SAVING, LOADING                            */
/**********************************************************************/

// set m_posErrFactor, suppose that Quantity index exist (private method)
void Catalog::setPosErrFactor(const int index) {

  std::string text=m_quantities[index].m_unit;
  if ( text.empty() ) {
    std::ostringstream sortie;
    sortie << "Generic position error (" << m_quantities[index].m_name
           << ") has no unit, by default multiply by: 1/" << m_posErrFactor;
    printWarn("private setGeneric", sortie.str());
  }
  else {
    // use explicit cast to be sure compiler choose the right tolower()
    std::transform(text.begin(), text.end(), text.begin(),
                   (int(*)(int)) tolower);
    if (text == "deg") m_posErrFactor=1.0;
    else if (text == "arcsec") m_posErrFactor=3600.0;
    else if (text == "arcmin") m_posErrFactor=60.0;
    else if (text == "rad") m_posErrFactor=1.0/Angle_Conv;
    else {
      std::ostringstream sortie;
      sortie << "Generic position error (" << m_quantities[index].m_name
             << ") has unknown unit, by default multiply by: 1/" << m_posErrFactor;
      printWarn("private setGeneric", sortie.str());
    }
  }
}
// search for generic quantities, private: suppose import done (private method)
void Catalog::setGeneric(const int whichCat) {

  const std::string epoch="J2000";
  int max=m_quantities.size();
  int j;
  std::string text, name;
  m_indexErr= -1;
  m_indexRA = -1;
  m_indexDEC= -1;
  std::vector<Quantity>::iterator itQ=m_quantities.begin();
  for (int i=0; itQ != m_quantities.end(); ++itQ, ++i) {

    if ((whichCat >= 0) && (whichCat < MAX_CAT)) {
      // if it is a known catalog, the generic are defined
      text=itQ->m_name;
      for (j=0; j<MAX_GEN; j++) {
        name=Catalog::s_CatalogGeneric[whichCat][j];
        if (name == "+") name=UCD_Added[j];
        if (text == name) {
          itQ->m_isGeneric=true;
          switch (j) {
            case 1: m_indexRA=i;
            break;
            case 2: m_indexDEC=i;
            break;
            case 3: m_indexErr=i;
                    setPosErrFactor(i);
            break;
            default:
              if (  itQ->m_ucd.empty() && j ) itQ->m_ucd=UCD_List[j];
            break;
          }
        }
      }// loop on generic
      text=itQ->m_ucd;
    }

    else {
      // otherwise, take first matching UCD
      text=itQ->m_ucd;
      if (text != "") for (j=0; j<MAX_GEN; j++) { 
        if (text == UCD_List[j]) {
          itQ->m_isGeneric=true;
          if ((j == 1) || (j == 2)) {
            // check that epoch is J2000
            name=itQ->m_name;
            name.erase(0, name.length()-epoch.length());
            // format contains at least 1 char
            if ((name == epoch) && (itQ->m_type == Quantity::NUM)) {
              if (j == 1) m_indexRA=i;
              else        m_indexDEC=i;
            }
            else itQ->m_isGeneric=false;
          }
          break;
        }
      }// loop on generic
    }

    // flag error quantities, useless for already found generic
    if ((text == "ERROR") && (!itQ->m_isGeneric)) {
      // error column name is e_"associated column" except for position error
      text=itQ->m_name;
      if (text.find("e_") == 0) {
        text.erase(0, 2);
        #ifdef DEBUG_CAT
        std::cout << "ERROR column (" << text << ")" << std::endl;
        #endif
        for (j=0; j<max; j++) if (m_quantities.at(j).m_name == text) {
          itQ->m_statError=text;
          break;
        }
      }
      else if ((text == "PosErr") || (text == "ErrorRad")) {
        itQ->m_isGeneric=true;
        m_indexErr=i;
        setPosErrFactor(i);
      }
    }

  }// loop on quantities

  /* since RA, DEC columns can be erased while importSelected() */
  /* must update the existence of elliptical region */
  if ((m_indexRA < 0) || (m_indexDEC < 0)) m_selRegion=false;
}

/**********************************************************************/
// return the number of RAM bytes needed for numRows catalog rows
long Catalog::getRAMsize(const long numRows, const bool writeLog) {

  std::string mot, text;
  int quantSize=m_quantities.size();
  if (quantSize == 0) return IMPORT_NEED;

  int nD=0, nS=0, nchar=0;
  int i, j;
  std::vector<Quantity>::iterator itQ=m_quantities.begin();
  for (; itQ != m_quantities.end(); ++itQ) {
    if (itQ->m_type == Quantity::NUM) nD++;
    else if  (itQ->m_type == Quantity::STRING) {
      nS++;
      text=itQ->m_format;
      i=text.length();
      if (i > 0) {
        if (m_URL == "CDS") text.erase(0, 1);
        else { /* from binary, first char can be letter because 1 is optional */
          if (i > 1) text.erase(i-1); else text="1";
        }
        j=std::atoi(text.c_str());
        if (j > 0) nchar+=j;
      }
    }
    else if (itQ->m_type == Quantity::VECTOR) {
    }
    else { //Quantity::LOGICALS)
      nS++;
      text=itQ->m_format;
      j=std::atoi(text.c_str());
      if (j > 0) nchar+=j;
    }
  }
  char buffer[50];
  long sizeD=nD*sizeof(double)*numRows;
  long sizeS=nchar*sizeof(char)*numRows;
  i=1+(quantSize+1)/(sizeof(long)*8);
  long sizeB=i*sizeof(long)*numRows;
  if (writeLog && (m_numOriRows > 0)) {
    sprintf(buffer, "%6ld", m_numOriRows);
    mot=buffer; /* convert C string to C++ string */
    text="Original whole catalog number of rows = "+mot;
    printLog(1, text);
  }
  if (numRows <= 0) return 0;
  if (writeLog) {
    sprintf(buffer, "%6ld", numRows);
    mot=buffer; /* convert C string to C++ string */
    text="Needed RAM space (MB) for "+mot;
    sprintf(buffer, "%5.1f", (sizeD+sizeS+sizeB)/(1024.*1024.));
    mot=buffer;
    text=text+" data rows = "+mot+"\n";
    sprintf(buffer, "%5.0f kB for numericals (%3d double per row)",
                   sizeD/1024., nD);
    mot=buffer;
    text=text+mot+"\n";
    sprintf(buffer, "%5.0f kB for %2d strings (%3d char per row)",
                   sizeS/1024., nS, nchar);
    mot=buffer;
    text=text+mot+"\n";
    sprintf(buffer, "%5.0f kB for select bits (%2d long per row)",
                   sizeB/1024., i);
    mot=buffer;
    text=text+mot;
    printLog(1, text);
  }
  return (sizeD+sizeS+sizeB);
}


/**********************************************************************/
// creates a new column in m_strings, m_numericals (private method)
void Catalog::create_tables(const int nbQuantNum, const long maxRows) {

  int vecSize;
  std::string errText;
  if (nbQuantNum > 0) {
    try {
      m_numericals.resize(nbQuantNum);
    }
    catch (const std::exception &err) {
      errText=std::string("EXCEPTION on m_numericals: ")+err.what();
      printErr("private create_tables", errText);
      throw;
    }
  }
  vecSize=m_quantities.size()-nbQuantNum;
// printf("sizes = %d , %d\n", nbQuantNum, vecSize);
  if (vecSize > 0) {
    try {
      m_strings.resize(vecSize);
    }
    catch (const std::exception &err) {
      errText=std::string("EXCEPTION on m_strings: ")+err.what();
      printErr("private create_tables", errText);
      throw;
    }
  }
  if (maxRows > 0) add_rows(maxRows);

}
/**********************************************************************/
// creates a new row in m_strings, m_numericals (private method)
void Catalog::add_rows(const long maxRows) {

  int i, vecSize;
  std::string errText;
// printf("!! %ld / %ld \n", m_numRows, maxRows);
  vecSize=m_strings.size();
  if (vecSize > 0 ) {
    try {
      for (i=0; i<vecSize; i++) {
        m_strings[i].resize(maxRows);
//        for (j=0; j<maxRows; j++) m_strings[i].at(j).resize(20);
      }
    }
    catch (const std::exception &err) {
      errText=std::string("EXCEPTION on m_strings: ")+err.what();
      printErr("private add_rows", errText);
      throw;
    }
  }
  vecSize=m_numericals.size();
  if (vecSize > 0 ) {
    try {
      for (i=0; i<vecSize; i++) m_numericals[i].resize(maxRows);
    }
    catch (const std::exception &err) {
      errText=std::string("EXCEPTION on m_numericals: ")+err.what();
      printErr("private add_rows", errText);
      throw;
    }
  } 

}


/**********************************************************************/
// read the catalog from fits file (only description if getDescr is true)
int Catalog::analyze_fits(const Table *myDOL, const bool getDescr,
                          const std::string origin, long *maxRows) {
  std::string text, mot;
  int   i, j, max,
        maxLogSize=0,
        maxVecSize=0,
        nbQuantNum=0;
/*  char  name[9];  8 char maximum for header key */
  bool  binary=true;
  short test;
  unsigned int pos;
  int (*pfunc)(int)=toupper; // function used by transform
  const Header &header=myDOL->getHeader();

  try {
    text=myDOL->getName();
    header.getKeyword("XTENSION", mot);
  /* text.rfind('_');
    since '/' are replaced by '_' in EXTNAME
    and  is ambiguous with '_' in catalog name:
    must search for keywords CDS-CAT  CDS-NAME
  */
  }
  catch (const TipException &x) {
    text=": fits EXTENSION, cannot get name";
    printErr(origin, text);
    return BAD_FITS;
  }
  std::ostringstream sortie;
  sortie << "fits extension " << mot << " name = " << text;
  printLog(1, sortie.str());
  sortie.str(""); // Will empty the string.

  m_catName="";   m_catRef="";
  m_tableName=""; m_tableRef="";
  m_URL="";
  if (mot == "TABLE") {
    binary=false;
    m_URL="CDS"; /* to know that format string need to be changed for fits */
  }
  try {
    test=0;
//    Header::Iterator itor_found;
    for (Header::ConstIterator itor=header.begin();itor!=header.end();++itor) {
      mot=itor->getName();
      // long comment from wanted key
      if ( mot.empty() ) { switch (test) {
        case 1: text=itor->getComment();
/*                m_catRef=itor_found->getComment()+text;*/
                m_catRef=text;
                break;
        case 2: text=itor->getComment();
/*                m_tableRef=itor_found->getComment()+text;*/
                m_tableRef=text;
                break;
        default: break;
      } }
      else {
        if (mot == KeyCDS[0]) {
          test=1; //itor_found=itor;
          m_catName=itor->getValue();
        }
        else if (mot == KeyCDS[1]) {
          test=2; //itor_found=itor;
          m_tableName=itor->getValue();
        }
        else if (test == 2) break;
      }// end of non-empty key
    }
  }
  catch (const TipException &x) {
    text=": fits EXTENSION, cannot read header";
    printErr(origin, text);
    return BAD_FITS;
  }
  if (test < 2) {
    if (binary) {
      printWarn(origin, "fits EXTENSION, cannot get CDS header keys");
    }
    else {
      printErr(origin, ": fits EXTENSION, cannot get CDS header keys");
      return BAD_FITS;
    }
  }
//  const Table::FieldCont &allCol=myDOL->getValidFields();
  try {
    max=myDOL->getValidFields().size();
    m_numOriRows=myDOL->getNumRecords();
  }
  catch (const TipException &x) {
    printErr(origin, ": fits EXTENSION, cannot get number of rows or columns");
    return BAD_FITS;
  }
  if (max < 1) {
    printErr(origin, ": fits EXTENSION, need at least 1 column");
    return BAD_FILETYPE;
  }
  text="";
  std::vector<double> colNull(max, 0.0);
  std::vector<int>    colSize(max, 0);
  try {
    const IColumn *myCol = 0;
    for (i=0; i < max; i++) {
      text="";
      myCol=myDOL->getColumn(i);
      Quantity readQ;
      readQ.m_name=myCol->getId();
      readQ.m_unit=myCol->getUnits();
//      sprintf(name, "TFORM%d", i+1);
//      mot=name; /* convert C string to C++ string */
      myCol->getColumnKeyword("TFORM").get(text);
//      header.getKeyword(mot, text);
      std::transform(text.begin(), text.end(), text.begin(), pfunc);
      j=text.length();
      if (j == 0) {
        text="unauthorized FITS empty format";
        sortie << ": FITS format is empty, column#" << i+1;
        throw std::runtime_error(text);
      }
      readQ.m_format=text;
      text="";
      if (binary) {
        if ( !myCol->isScalar() ) {
          if (readQ.m_format[j-1] == 'L') {
            readQ.m_type=Quantity::LOGICALS;
            j=atoi(readQ.m_format.c_str());
            colSize[i]=j;
            if (j > maxLogSize) maxLogSize=j;
          }
          else {
            readQ.m_type=Quantity::VECTOR;
            j=atoi(readQ.m_format.c_str());
            colSize[i]=j;
            if (j > maxVecSize) maxVecSize=j;
            text="skipping FITS TABLE vector";
            sortie << "VECTOR not supported, unusable column#" << i+1;
            printWarn(origin, sortie.str() );
            sortie.str("");
          }
        }
        else if ((readQ.m_format[j-1] == 'A') || (readQ.m_format[j-1] == 'L'))
          readQ.m_type=Quantity::STRING;
        else readQ.m_type=Quantity::NUM;
/*        sprintf(name, "TNULL%d", i+1); mot=name;
        try { header.getKeyword(mot, readQ.m_null); } */
        try { myCol->getColumnKeyword("TNULL").get(readQ.m_null); }
        catch (const TipException &x) {}
        colNull.at(i)=atof(readQ.m_null.c_str());
        text="";
/*          sprintf(name, "TSCAL%d", i+1); mot=name;
          try { header.getKeyword(mot, text); } */
        try { myCol->getColumnKeyword("TSCAL").get(text); }
        catch (const TipException &x) {}
        readQ.m_null+=SepNull+text;
        if ( !text.empty() && (readQ.m_null[0]!=SepNull) )
          colNull.at(i)*=atof(text.c_str());
        text="";
/*          sprintf(name, "TZERO%d", i+1); mot=name;
          try { header.getKeyword(mot, text); } */
        try { myCol->getColumnKeyword("TZERO").get(text); }
        catch (const TipException &x) {}
        readQ.m_null+=SepNull+text;
        if ( !text.empty() && (readQ.m_null[0]!=SepNull) )
          colNull.at(i)+=atof(text.c_str());
        text="";
        try {
          myCol->getColumnKeyword(Key_UCD).get(text);
          readQ.m_comment=myCol->getColumnKeyword(Key_UCD).getComment();
        } 
        catch (const TipException &x) {
          sortie << "missing keyword " << Key_UCD << " for column#" << i+1;
          printWarn(origin, sortie.str() );
          sortie.str(""); // Will empty the string.
        }
        /* could also read TDISP keyword, to have better ASCII output */
      }
      else {
        readQ.m_null=SepNull;
        readQ.m_null+=SepNull;
        /* according to standards: Aw or Iw or Fw.d or Ew.d or Dw.d */
        if (readQ.m_format[0] == 'A') readQ.m_type=Quantity::STRING;
        else readQ.m_type=Quantity::NUM;
/*        sprintf(name, "TBCOL%d", i+1); mot=name;
        text=header.getKeyComment(mot); */
        text=myCol->getColumnKeyword("TBCOL").getComment();
        j=4;
        if (text.substr(0,j) == "UCD=") text.erase(0,j);
        pos=text.find('.');
        if (pos != std::string::npos) text.erase(pos);
        /* ONLY  readQ.m_comment  IS NOT SET */
      }
      std::transform(text.begin(), text.end(), text.begin(), pfunc);
      readQ.m_ucd=text;
      if (readQ.m_type == Quantity::NUM) {
        readQ.m_index=nbQuantNum;
        nbQuantNum++;
      }
      else readQ.m_index=m_quantities.size()-nbQuantNum;
      // for the moment, VECTOR considered as STRING
      try { m_quantities.push_back(readQ); }
      catch (const std::exception &prob) {
        text=std::string("EXCEPTION filling m_quantities: ")+prob.what();
        sortie << text;
        throw;
      }
    }/* loop on columns */
  }
  catch (...) {
    m_quantities.clear();
    j=BAD_FILETYPE;
    m_numRows=j;
    if ( text.empty() ) {
      sortie << ": fits EXTENSION, error reading TABLE column#";
      sortie << i+1 << " description";
      printErr(origin, sortie.str() );
    }
    else {
      printErr(origin, sortie.str() );
      // exception sent only for vector allocation
      if (text == sortie.str() ) throw;
    }
    return j;
  }
// printf("numRows=%ld\n", m_numOriRows );
  m_numRows=0;
  if ( !*maxRows ) *maxRows=m_numOriRows;
  if ((getDescr) || (m_numOriRows == 0)) return IS_OK;

  create_tables(nbQuantNum, *maxRows);
  try {
    double rowVal;
/*when reading a Logical, string gives "T" or "F" while char gives '0' or '1'*/
    std::vector<char>  logic(maxLogSize, ' ');
    std::vector<double> vect(maxVecSize);
//printf("maxVecSize = %d\n", maxVecSize);
    std::vector<Quantity>::iterator itQ;
    // Loop over all records (rows) and extract values
    for (Table::ConstIterator itor=myDOL->begin(); itor != myDOL->end();
         ++itor) {
      i=0;
      for (itQ=m_quantities.begin(); itQ != m_quantities.end(); ++itQ, ++i) {
        j=itQ->m_index;
        // double variable to hold the value of all the numeric fields
        if (itQ->m_type == Quantity::NUM) {
          rowVal=(*itor)[itQ->m_name].get();
          if ( itQ->m_null[0]==SepNull ) m_numericals[j].at(m_numRows)=rowVal;
          else {
            if (rowVal == colNull.at(i)) m_numericals[j].at(m_numRows)=MissNAN;
            else m_numericals[j].at(m_numRows)=rowVal;
          }
        }
        else if (itQ->m_type == Quantity::STRING)
          (*itor)[itQ->m_name].get(m_strings[j].at(m_numRows));
        else if (itQ->m_type == Quantity::VECTOR) {
        }
        else { //Quantity::LOGICALS)
          (*itor)[itQ->m_name].get(logic); // function resizes the vector
          if (colSize.at(i) != (int)logic.size()) {
            text="read vector size differs from expected";
            printErr(origin, text);
            throw std::runtime_error(text);
          }
          text="";
          for (nbQuantNum=0; nbQuantNum < colSize.at(i); nbQuantNum++) {
//            printf("%d", logic[nbQuantNum]);
            switch (logic[nbQuantNum]) {
              case 0: text+="F"; break; case 1: text+="T"; break;
              default: text+=" ";
            }
          }
          m_strings[j].at(m_numRows)=text;
//printf("logicals[%d] row%03ld = %s.\n", colSize.at(i),m_numRows, text.c_str() );
        }
      }
      if (++m_numRows == *maxRows) break;
    }/* loop on rows*/
  }
//catch (const TipException &x) {
  catch (...) {
    sortie << ": fits EXTENSION, cannot read after row#" << m_numRows+1;
    printErr(origin, sortie.str() );
    return BAD_ROW;
  }
  return IS_OK;
}
/**********************************************************************/
/* PRIVATE METHODS analyze_head, analyze_body in file "catalog_ioText.cxx" */
/**********************************************************************/
// common code between import and importDescription (private method)
int Catalog::load(const std::string &fileName, const std::string ext,
                  const bool getDescr, long *maxRows) {

  int i, err=IS_OK;
  std::string origin, text;
  if (getDescr) origin="importDescription"; else origin="import";

  m_numRows=0;
  i=fileName.length();
  if (i == 0) {
    text=": FILENAME is EMPTY";
    printErr(origin, text);
    err=BAD_FILENAME;
    m_numRows=err;
    return err;
  }
  std::fstream myFile (fileName.c_str(), std::ios::in);
  // file can be opened ?
  if ( !myFile.is_open() ) {
    text=": FILENAME \""+fileName+"\" cannot be opened";
    printErr(origin, text);
    err=BAD_FILENAME;
    m_numRows=err;
    return err;
  }
/*  DIR *testDir=opendir(fileName.c_str());
  if (testDir != NULL) {
    closedir(testDir);
    text=": FILENAME \""+fileName+"\" is a directory";
    printErr(origin, text);
    err=BAD_FILETYPE;
    m_numRows=err;
    return err;
  }*/
  myFile.close();
  std::ostringstream sortie;
  bool  myTest;
  const Extension *myEXT = 0;
  // cannot use readTable because FITS IMAGE returns also error 1
  try {
    myEXT=IFileSvc::instance().readExtension(fileName, ext);
  }
  catch (const TipException &x) {
    err=x.code();
    if (err == 1) {
      // This non-cfitsio error number means the file does not exist
      // or is not a table, or is not either Root nor Fits format. 
      err=BAD_FITS;
/*      printWarn(origin, "FILENAME is NOT fits");*/
    }
    else {
      // Other errors come from cfitsio, but apply to a file
      // which is in FITS format but has some sort of format error.
      sortie << ": FILENAME is FITS, but cfitsio returned error=" << err;
      printErr(origin, sortie.str());
      err=BAD_FITS;
      m_numRows=err;
      delete myEXT;
      return err;
    }
  }
  if (err == IS_OK) myTest=myEXT->isTable();
  delete myEXT;
  if (err == IS_OK) {
 
    if (!myTest) {
      text=": FILENAME is FITS, but NOT a TABLE";
      printErr(origin, text);
      err=BAD_FITS;
      m_numRows=err;
      return err;
    }
    const Table *myDOL = 0;
    try {
      myDOL=IFileSvc::instance().readTable(fileName, ext);
    }
    catch (const TipException &x) {
      sortie << ": FITS is TABLE, but cfitsio returned error=" << x.code();
      printErr(origin, sortie.str());
      err=BAD_FITS;
      m_numRows=err;
      delete myDOL;
      return err;
    }
    const char lf = 0x0A;
    m_filename=fileName+lf+ext;
    err=analyze_fits(myDOL, getDescr, origin, maxRows);
    delete myDOL;

  }
  else { // (err == BAD_FITS)

    err=IS_OK;
    unsigned long tot=0ul; // number of lines read
    bool testCR=false; // true if lines ends with CR (on windows)
    int what=0;        // 0 for unkwown, >0 for csv, <0 to tsv
                       // fabs()=1 for standard, =2 for meta QUERY

    myFile.open(fileName.c_str(), std::ios::in);
    err=analyze_head(&tot, &what, &testCR, &myFile);
    if (!tot) {
      text=": FILENAME \""+fileName+"\" is fits without extension[] specified";
      printErr(origin, text);
      text="input file";
      if (err == IS_OK) err=BAD_FILENAME;
      // in case, for strange reason, file.good() fails at first time
    }
    else if (err < IS_OK) {
      sortie << ": FILENAME \"" << fileName << "\" is empty or "
             << "has unknown structure (stopped step " << -1*err << ")";
      printErr(origin, sortie.str());
      sortie.str(""); // Will empty the string.
      text="input file";
      if (err == 0) err=BAD_FILENAME; /* nothing is read */
      else err=BAD_FILETYPE; /* can search if catalog name exist */
    }
    else {

      //what is 1 or 2, then negate if TSV type
      if (what == 1) text="input text file";
      else text="input META file";

      // get columns description and units, data
      err=analyze_body(&tot, &what, testCR, getDescr, &myFile, maxRows);
      if (err == BAD_FILELINE) {
        // stopped before reading needed stuff (with getDescr false)
        text=": FILENAME \""+fileName+"\" couldn't be read (line too long)";
        printErr(origin, text);
      }
      else if (err < IS_OK) {
        sortie << ": FILENAME \"" << fileName
             << "\" has wrong type (stopped step " << 5-1*err << ")";
        printErr(origin, sortie.str());
        sortie.str(""); // Will empty the string.
        err=BAD_FILETYPE;
      }
      else {
        if (what > 0) sortie << text << " is CSV type (; separator)";
        else sortie << text << " is TSV type (Tab=0x09 separator)";
        printLog(1, sortie.str());
        sortie.str(""); // Will empty the string.
      }

    }
    myFile.close();
    sortie << text << ": " << tot << " lines read";
    printLog(0, sortie.str());
    sortie.str(""); // Will empty the string.
    m_filename=fileName;
    m_URL="CDS"; /* to know that format string need to be changed for fits */

  }// file is read
  if (err < 0) {
    deleteContent(); // to erase data already loaded in memory
    m_numRows=err;
    // exit only if nothing is read
    if ((err == BAD_FILENAME) || (err == BAD_FITS)) return err;
  }
  // m_quantities vector maybe filled, MUST fill m_selEllipse, m_loadQuantity
  try {
    m_selEllipse.assign(7, 0.0);
    m_loadQuantity.assign(m_quantities.size(), true);
  }
  catch (const std::exception &prob) {
    text=std::string("EXCEPTION on creating m_selEllipse, m_loadQuantity: ")
        +prob.what();
    printErr(origin, text);
    throw;
  }
  // check if tableName is known to set generic
  for (i=0; i<MAX_CAT; i++) {
    if (Catalog::s_CatalogList[2*i+1] == m_tableName) break;
  }
  if (i == MAX_CAT) {
    text="Unknown table name, all generic quantities may be not found";
    printWarn(origin, text);
  }
  else m_code=Catalog::s_CatalogList[2*i];
  setGeneric(i);
  return err;
}
/**********************************************************************/
// read only the catalog description from file
int Catalog::importDescription(const std::string &fileName,
                               const std::string ext) {
  int err;
  err=checkImport("importDescription", false);
  // must be called before load() which set m_numRows to 0
  if (err < IS_VOID) return err;

  long maxR=0l;
  err=load(fileName, ext, true, &maxR);
  if (err < IS_OK) return err;

  return m_quantities.size();
}
/**********************************************************************/
// read from file an entire catalog without selection
int Catalog::import(const std::string &fileName, const long maxRows,
                    const std::string ext) {
  int err;
  long maxR=maxRows;
  err=checkImport("import", false);
  if (err < IS_VOID) return err;

  if (maxR <= 0) {
    printWarn("import", "trying to get whole catalog file");
    maxR=0;
  }

  err=load(fileName, ext, false, &maxR);
  if (err < IS_OK) return err;
  getRAMsize(m_numRows, true);
  try {
    if (m_numRows < maxR) {
      int i;
      err=m_strings.size();
      // erase unused memory  printf("ERASING\n");
      for (i=0; i<err; i++) m_strings[i].resize(m_numRows);
      err=m_numericals.size();
      for (i=0; i<err; i++) m_numericals[i].resize(m_numRows);
    }
    err=m_quantities.size()+2;
    // number of required bits including global and region
    err=1+(err-1)/(sizeof(long)*8);
    m_rowIsSelected.resize(err);
    #ifdef DEBUG_CAT
    std::cout << "Number of unsigned long required for m_rowIsSelected = "
              << err << std::endl;
    #endif
    if (m_numRows)
      for (int j=0; j<err; j++) m_rowIsSelected[j].assign(m_numRows, 0);
    // above lines can be commented to test the try catch mechanism
  }
  catch (const std::exception &prob) {
    std::string text;
    text=std::string("EXCEPTION filling m_rowIsSelected: ")+prob.what();
    printErr("import", text);
    throw;
  }
  return m_numRows;
}


/**********************************************************************/
int Catalog::loadSelectFits(const std::string &fileName, const std::string ext,
                            long *maxRows, std::string &filter)
{
  const std::string  origin="importSelected";
  std::ostringstream sortie;
  unsigned int       i, max;
  const Table *myDOL = 0;

  /* create the selection string */
  char value[20];
  bool test, probCase=false;
  std::vector<bool> isSelected;
  if ( existCriteria(&isSelected) && (m_URL.empty()) ) {
    if (m_selRegion) {
      /* !! cannot concatenate Cstring and char !! */
      filter+='('+m_quantities[m_indexDEC].m_name+" >= ";
      sprintf(value, "%9.5f", m_selEllipseCentDEC_deg-m_selEllipseMajAxis_deg);
      filter+=value;
      filter+=" && "+m_quantities[m_indexDEC].m_name+" <= ";
      sprintf(value, "%9.5f", m_selEllipseCentDEC_deg+m_selEllipseMajAxis_deg);
      filter+=value;
      /* AND has higher priority than OR, needn't parenthesis (isnull...)*/
      if ( !m_quantities[m_indexRA].m_rejectNaN )
        filter+=" || isnull("+m_quantities[m_indexRA].m_name+')';
      if ( !m_quantities[m_indexDEC].m_rejectNaN )
        filter+=" || isnull("+m_quantities[m_indexDEC].m_name+')';
      filter+=')';
    }
    max=isSelected.size();
    unsigned int j, nbV;
    Quantity readQ;
    double   rVal;
    // first boolean is for elliptical region (m_selRegion)
    for (i=1; i<max; i++) if ( isSelected[i] ) {
      if ( !filter.empty() ) {
        if (m_criteriaORed) filter+=" || "; else filter+=" && ";
      }
      filter+='(';
      readQ=m_quantities[i-1];
      if (readQ.m_type == Quantity::STRING) {
        nbV=readQ.m_listValS.size(); /* nbV is >0 because criteria exist */
        test=readQ.m_excludeList;
        filter+='('+readQ.m_name;
        if (test) filter+=" != '"; else filter+=" == '";
        filter+=readQ.m_listValS[0]+"')";
        for (j=1; j<nbV; j++) {
          if (test) filter+=" && "; else filter+=" || ";
          filter+='('+readQ.m_name;
          if (test) filter+=" != '"; else filter+=" == '";
          filter+=readQ.m_listValS[j]+"')";
        }
        if ((!m_criteriaORed) && (!readQ.m_cutORed)) probCase=true;
      }
      else if (readQ.m_type == Quantity::NUM) {
        if (readQ.m_rejectNaN == false) {
          filter+="isnull("+readQ.m_name+") || ";  
          /* can put OR because criteria exist */
          /* AND has higher priority than OR, needn't parenthesis (isnull...)*/
        }
        if (readQ.m_lowerCut < NO_SEL_CUT) {
          filter+='('+readQ.m_name+" >= ";          
          sprintf(value, "%.9E", readQ.m_lowerCut);
          filter+=value; filter+=')';
          test=true;
        }
        else test=false;
        if (readQ.m_upperCut < NO_SEL_CUT) {
          if (test) filter+=" && ";
          filter+='('+readQ.m_name+" <= ";          
          sprintf(value, "%.9E", readQ.m_upperCut);
          filter+=value; filter+=')';
          test=true;
        }
        nbV=readQ.m_listValN.size();
        if (nbV) { 
          if (test) {
            if (readQ.m_cutORed) filter+=" || ("; else filter+=" && (";
          }
          if (readQ.m_excludeList) filter+="!(";
          filter+="near("+readQ.m_name+", ";
          rVal=readQ.m_listValN[0];
          sprintf(value, "%.9E", rVal);
          filter+=value; filter+=", ";
          if (fabs(rVal) > NearZero) rVal*=readQ.m_precision;
          else rVal=readQ.m_precision;
          sprintf(value, "%.9E", rVal);
          filter+=value; filter+=')';
          for (j=1; j<nbV; j++) {
            filter+=" || near("+readQ.m_name+", ";
            rVal=readQ.m_listValN[j];
            sprintf(value, "%.9E", rVal);
            filter+=value; filter+=", ";
            if (fabs(rVal) > NearZero) rVal*=readQ.m_precision;
            else rVal=readQ.m_precision;
            sprintf(value, "%.9E", rVal);
            filter+=value; filter+=')';
          }
          if (readQ.m_excludeList) filter+=")";
          if (test) filter+=")";
        }
      }
      filter+=')';
    }
    if (probCase) printWarn(origin,
      "filter (for binary fits) is NOT case sensitive for string");
  }
  try {
    myDOL=IFileSvc::instance().readTable(fileName, ext, filter);
  }
  catch (const TipException &x) {
    sortie << ": FITS is TABLE, but cfitsio returned error=" << x.code();
    printErr(origin, sortie.str());
    delete myDOL;
    return BAD_FITS;
  }
  try {
    max=myDOL->getValidFields().size();
    m_numOriRows=myDOL->getNumRecords();
  }
  catch (const TipException &x) {
    printErr(origin, ": fits EXTENSION, cannot get number of rows or columns");
    delete myDOL;
    return BAD_FITS;
  }
  if (max < 1) {
    printErr(origin, ": fits EXTENSION, need at least 1 column");
    delete myDOL;
    return BAD_FILETYPE;
  }
  /* beware some initial quantities can be skipped */
  if (m_loadQuantity.size() != max) {
    sortie << ": fits EXTENSION, number of fields (" << max;
    sortie << ") differ from previous import call (";
    sortie << m_loadQuantity.size() << ")";
    printErr(origin, sortie.str());
    delete myDOL;
    return BAD_FILETYPE;
  }
//std::cout << "Initial number of COL = " << max << std::endl;
  std::string text, unit, rescale, form;
  int  err=0;
/*  const Header &header=myDOL->getHeader();
    char name[9];  8 char maximum for header key */
  std::vector<double> colNull(max, 0.0);
  std::vector<Quantity>::iterator itQ=m_quantities.begin();
  try {
    const IColumn *myCol = 0;
    for (i=0; i < max; i++) {
      text="start";
      rescale="";
      myCol=myDOL->getColumn(i);
      text=myCol->getId();
      unit=myCol->getUnits();
/*      sprintf(name, "TNULL%d", i+1); mot=name;
      try { header.getKeyword(mot, rescale); } */
      try { myCol->getColumnKeyword("TNULL").get(rescale); }   
      catch (const TipException &x) {}
      colNull.at(i)=atof(rescale.c_str());
      form="";
/*          sprintf(name, "TSCAL%d", i+1); mot=name; */
      try { myCol->getColumnKeyword("TSCAL").get(form); }
      catch (const TipException &x) {}
      rescale+=SepNull+form;
      if ( !form.empty() && (rescale[0]!=SepNull) ) colNull.at(i)*=atof(form.c_str());
      form="";
/*          sprintf(name, "TZERO%d", i+1); mot=name; */
      try { myCol->getColumnKeyword("TZERO").get(form); }
      catch (const TipException &x) {}
      rescale+=SepNull+form;
      if ( !form.empty() && (rescale[0]!=SepNull) ) colNull.at(i)+=atof(form.c_str());
      myCol->getColumnKeyword("TFORM").get(form);
      if (m_loadQuantity[i]) {
//std::cout << itQ->m_index << " ";
        if ( (text!=itQ->m_name) || (unit!=itQ->m_unit)
            || (rescale!=itQ->m_null) || (form!=itQ->m_format) ) {
          text="";
          throw std::runtime_error("fileDiff");
        }
        if (itQ->m_type == Quantity::NUM) err++;
        itQ++;
      }
    }/* loop on columns */
  }
  catch (...) {
    if ( text.empty() ) {
      sortie << ": fits EXTENSION, column#";
      sortie << i+1 << " differ from previous import call";
      printErr(origin, sortie.str() );
      err=BAD_FILETYPE;
    }
    else {
      sortie << ": fits EXTENSION, error reading TABLE column#";
      sortie << i+1 << " description";
      printErr(origin, sortie.str() );
      err=BAD_FITS;
    }
    delete myDOL;
    return err;
  }
//std::cout << std::endl;
//std::cout << m_numOriRows << " OriRows" << std::endl;
  m_numRows=0;
  if ( !*maxRows ) *maxRows=m_numOriRows;
  if ( !m_numOriRows ) return IS_OK;

  create_tables(err, *maxRows);
  try {
    double rowVal;
    std::vector<char>  logic;
    std::vector<double> vect;
    for (Table::ConstIterator itor=myDOL->begin(); itor != myDOL->end();
         ++itor, ++m_numRows) {
      err=0;
      for (itQ=m_quantities.begin(); itQ != m_quantities.end(); ++itQ, ++err) {
        i=itQ->m_index;
        if (itQ->m_type == Quantity::NUM) {
          rowVal=(*itor)[itQ->m_name].get();
          if ( itQ->m_null[0]==SepNull ) m_numericals[i].at(m_numRows)=rowVal;
          else {
            if (rowVal == colNull[err]) m_numericals[i].at(m_numRows)=MissNAN;
            else m_numericals[i].at(m_numRows)=rowVal;
          }
        }
        else if (itQ->m_type == Quantity::STRING)
          (*itor)[itQ->m_name].get(m_strings[i].at(m_numRows) );
        else if (itQ->m_type == Quantity::VECTOR) {
        }
        else { //Quantity::LOGICALS)
          (*itor)[itQ->m_name].get(logic); // function resizes the vector
          text="";
          for (max=0; max < logic.size(); max++) {
            switch (logic[max]) {
              case 0: text+="F"; break; case 1: text+="T"; break;
              default: text+=" ";
            }
          }
          m_strings[i].at(m_numRows)=text;
//printf("logicals[%ud] row%03ld = %s.\n", max,m_numRows, text.c_str() );
        }
      }
    }/* loop on rows*/
  }
  catch (const TipException &x) {
    sortie << ": fits EXTENSION, cannot read after row#" << m_numRows+1;
    printErr(origin, sortie.str() );
    return BAD_ROW;
  }
  delete myDOL;
  return IS_OK;
}
/**********************************************************************/
/* PRIVATE METHOD loadSelected is in file "catalog_ioText.cxx"        */
/**********************************************************************/
// if a catalog description was already loaded, this method does
// the same as import() adding selection criteria for loading
int Catalog::importSelected(std::string &filter) {

  const std::string origin="importSelected";
  int quantSize=checkImport(origin, true);
  if (quantSize < IS_VOID) return quantSize;

  if (m_numRows > 0) {
    printWarn(origin, "call 'deleteContent' before 'importSelected'");
    return IMPORT_BIS;
  }
  int i, err=0,
      maxSize=m_loadQuantity.size();
  for (i=0; i<maxSize; i++) if (m_loadQuantity[i]) err++;
  if (err < 2) {
    printErr(origin, "at least 2 quantities must be selected for import");
    return BAD_SEL_QUANT;
  }
  filter="";
  std::ostringstream sortie; 
  sortie << err << " quantities (over " << maxSize << ") selected for import";
  printLog(1, sortie.str());
  sortie.str("");
  if (err != quantSize) {
    // in fact, due to select*Quantities() restrictions,
    // quantSize can only be greater than err
    printLog(2, "Erasing unwanted quantities for importSelected()");
    deleteQuantities();
  }
  err=IS_VOID;
  std::string text;
  long maxRows=0l;
  unsigned int pos=m_filename.find(0x0A);
  if (pos == std::string::npos) {

    if ( m_filename.empty() ) {
      // data must be querried via web


    }
    else {
      if (!m_filePos) {
        printErr(origin, "file import was not succesful");
        return IMPORT_NEED;
      }
      std::fstream myFile (m_filename.c_str(), std::ios::in);
      // file can be opened ?
      if ( !myFile.is_open() ) {
        text=": FILENAME \""+m_filename+"\" cannot be opened";
        printErr(origin, text);
        err=BAD_FILENAME;
        m_numRows=err;
        return err;
      }
      // go back to initial position
      myFile.clear(); // needed to reset flags before seekg
      myFile.seekg(m_filePos);
      unsigned long tot=0ul; // number of data lines read
      err=loadSelected(&tot, &myFile, &maxRows);
      if (err == IS_OK) {
        sortie << tot << " data lines read for importSelected()";
        printLog(0, sortie.str());
        sortie.str(""); // Will empty the string.
      }
      myFile.close();
    }

  }
  else { //fits file

    text=m_filename;
    std::string fileName=text.substr(0, pos);
    text.erase(0, pos+1);
    const Extension *myEXT = 0;
    // cannot use readTable because FITS IMAGE returns also error 1
    try {
      myEXT=IFileSvc::instance().readExtension(fileName, text);
    }
    catch (const TipException &x) {
      err=x.code();
      if (err == 1) {
        // This non-cfitsio error number means the file does not exist
        // or is not a table, or is not either Root nor Fits format. 
        sortie << ": FILENAME \"" << fileName;
        sortie << "\" cannot be opened or is NOT fits";
        err=BAD_FILENAME;
      }
      else {
        // Other errors come from cfitsio, but apply to a file
        // which is in FITS format but has some sort of format error.
        sortie << ": FILENAME is FITS, but cfitsio returned error=" << err;
        err=BAD_FITS;
      }
    }
    if (err >= IS_VOID) {
      if ( !myEXT->isTable() ) {
        sortie << ": FILENAME is FITS, but NOT a TABLE";
        err=BAD_FITS;
      }
    }
    delete myEXT;
    if (err < IS_VOID) {
      printErr(origin, sortie.str());
      m_numRows=err;
      return err;
    }
    err=loadSelectFits(fileName, text, &maxRows, filter);

  }
  if (err < IS_VOID) {
    deleteContent(); // to erase data already loaded in memory
    m_numRows=err;
    return err;
  }
  try {
    if (m_numRows < maxRows) {
      err=m_strings.size();
      // erase unused memory  printf("ERASING\n")
      for (i=0; i<err; i++) m_strings[i].resize(m_numRows);
      err=m_numericals.size();
      for (i=0; i<err; i++) m_numericals[i].resize(m_numRows);
    }
    err=m_quantities.size()+2;
    // number of required bits including global and region
    maxSize=sizeof(long)*8;
    err=1+(err-1)/maxSize;
    m_rowIsSelected.resize(err);
    #ifdef DEBUG_CAT
    std::cout << "Number of unsigned long required for m_rowIsSelected = "
              << err << std::endl;
    #endif
    if (m_numRows) {
      std::vector<bool> isSelected;
      if ( existCriteria(&isSelected) ) {

        unsigned long quantBit;
        int j, k;
        quantSize=isSelected.size();
        if ((!m_URL.empty()) || (m_criteriaORed)) {
          /* if region selected, must suppose that criteria is fulfilled */
          if (isSelected[0]) quantBit=2ul; else quantBit=0ul;
          m_rowIsSelected[0].assign(m_numRows, quantBit);
          for (j=1; j<err; j++) m_rowIsSelected[j].assign(m_numRows, 0);
          /* for the moment, no selection applied */
          std::vector<std::string> list;
          std::string mot;
          int nbV, (*pfunc)(int)=tolower; // function used by transform
          Quantity readQ;
          double   myVal, precis;
          bool     check, reject;
          for (i=1; i<quantSize; i++) if (isSelected[i]) {
            readQ=m_quantities[i-1];
            quantBit=bitPosition(i-1, &k);
//std::cout << m_criteriaORed <<"byte #"<< k << " unsigned long="<<quantBit;
            if (readQ.m_type == Quantity::STRING) {
              nbV=readQ.m_listValS.size(); /* >0 because criteria exist */
              list.assign(nbV, "");
              for (j=0; j<nbV; j++) {
                mot=readQ.m_listValS.at(j);
                if (readQ.m_cutORed)
                  std::transform(mot.begin(), mot.end(), mot.begin(), pfunc);
                list.at(j)=mot;
              }
              for (maxRows=0; maxRows<m_numRows; maxRows++) {
                mot=m_strings[readQ.m_index].at(maxRows);
                if (readQ.m_cutORed)
                  std::transform(mot.begin(), mot.end(), mot.begin(), pfunc);
                check=readQ.m_excludeList;
                for (j=0; j<nbV; j++) {
                  if (mot == list[j]) {check=!readQ.m_excludeList; break;}
                }
                // setting required bit
                if (check)
                  m_rowIsSelected[k].at(maxRows)|= quantBit;
                else
                  m_rowIsSelected[k].at(maxRows)&= (Max_Test-quantBit);
              }
            }
            else { // numerical quantity, VECTOR cannot exist in string data
              j=i-1;
              nbV=readQ.m_index;
              reject=readQ.m_rejectNaN;
              precis=readQ.m_precision;
              for (maxRows=0; maxRows<m_numRows; maxRows++) {
                myVal=m_numericals[nbV].at(maxRows);
                if (!readQ.m_cutORed) // usual case
                  check=checkNUM(myVal,j, readQ.m_excludeList,reject,precis);
                else
                  check=checkNUMor(myVal,j, reject,precis);
                // setting required bit
                if (check)
                  m_rowIsSelected[k].at(maxRows)|= quantBit;
                else
                  m_rowIsSelected[k].at(maxRows)&= (Max_Test-quantBit);
              }// loop on rows
            }
          }/*loop on selected quantities */
          m_numSelRows=0;
          for (maxRows=0; maxRows<m_numRows; maxRows++) {
            if (rowSelect(maxRows, isSelected) == true) m_numSelRows++;
          }
        }
        else {// from binary fits, only if criteria ANDed
          m_numSelRows=m_numRows;
          /* easy, as criteria are ANDed: bit to 1 */
          std::vector<unsigned long> selBits(err, 0);
          selBits[0]=1ul;
          j=0;
          quantBit=2ul;
          for (i=0; i<quantSize; ) {
            if (isSelected[i]) selBits[j]|=quantBit;
            if ((++i % maxSize) == 0) {
              quantBit=1ul; // first bit
              j++;          // of next vector index.
            }
            else quantBit*=2ul;
          }
//std::cout << m_criteriaORed <<"first unsigned long selection="<< selBits[0];
          for (j=0; j<err; j++)
            m_rowIsSelected[j].assign(m_numRows, selBits[j]);
        }
        if (m_selRegion) {
          sortie << origin << ", selecting region from " << m_numRows;
          sortie << " loaded rows (" << m_numSelRows << " already selected)";
          printLog(1, sortie.str()); sortie.str(""); // Will empty the string.
          int nRA =m_quantities[m_indexRA].m_index,
              nDEC=m_quantities[m_indexDEC].m_index;
          bool check;
          for (maxRows=0; maxRows<m_numRows; maxRows++) {
            check=checkRegion(maxRows, nRA, nDEC);
            if (!check) {// setting 2nd bit to 0
              m_rowIsSelected[0].at(maxRows)&=(Max_Test-2ul);
              if (!m_criteriaORed) {// setting 1st bit to 0
                m_rowIsSelected[0].at(maxRows)&=(Max_Test-1ul);
                m_numSelRows--;
              }
              else if (rowSelect(maxRows, isSelected) == false) m_numSelRows--;
            }
          }/* loop on rows */
        }
        eraseNonSelected();

      }
      else  // lines can be commented to test the try catch mechanism
        for (int j=0; j<err; j++) m_rowIsSelected[j].assign(m_numRows, 0);
    }/* row exist ==> MUST create selection array */
  }
  catch (const std::exception &prob) {
    std::string text;
    text=std::string("EXCEPTION filling m_rowIsSelected: ")+prob.what();
    printErr("import", text);
    throw;
  }
  return m_numRows;
}


/**********************************************************************/
// create catalog header from memory to a FITS file
int Catalog::createFits(const std::string &fileName, const std::string &extName,
                        bool clobber, bool append, const std::string origin,
                        tip::Table **ptrTable) {

/*  char name[9];  8 char maximum for header key */
  std::string text, newName;
  int  j, err;
  bool create=true;

  err=checkImport(origin, true);
  if (err < IS_VOID) return err;

  err=fileName.length();
  if (err == 0) {
    text=": FILENAME is EMPTY";
    printErr(origin, text);
    return BAD_FILENAME;
  }
  if (extName.length() > 68) {
    text=": EXTENSION name is TOO long (limited to 68 characters)";
    printErr(origin, text);
    return BAD_FILENAME;
  }
  // overwrite existing file ?
  if (!clobber) {
    std::fstream file (fileName.c_str(), std::ios::in);
    if ( file.is_open() ) {
      file.close();
      create=false;
      if (!append) {
        text=": FILENAME \""+fileName+"\" exist (clobber=no, append=no)";
        printErr(origin, text);
        return BAD_FILENAME;
      }
      else {
        text="appending extension to existing FILENAME (method "+origin+")";
        printLog(1, text);
      }
    }
  }
  if ( extName.empty() ) {
    /* max length is 68, as keyword value starts/ends with ' in 11/80 */
    newName=m_tableName.substr(0,68);
    /* replace forbidden char. */
    err=newName.length();
    for (j=0; j < err; j++) {
      if (newName[j] == '/') newName[j]='_';
    }
    text="EXTENSION name set to '"+newName+"'";
    printWarn(origin, text);
  }
  else newName=extName;

  try {
    if (create) IFileSvc::instance().createFile(fileName);
    IFileSvc::instance().appendTable(fileName, newName);
    (*ptrTable)=IFileSvc::instance().editTable(fileName, newName);
  }
  catch (const TipException &x) {
    text=": cannot append and edit fits EXTENSION";
    printErr(origin, text);
    return BAD_FITS;
  }
  Header &header=(*ptrTable)->getHeader();
  try {
/*    header.setKeyword(KeyCDS[0], m_catName);
    header.setKeyComment(KeyCDS[0], KeyCom[0]);
UNUSBALE, otherwise append of long comment done after EXTNAME
*/
    if ( !m_catName.empty() ) {
      header.append( KeyRecord(KeyCDS[0],m_catName,KeyCom[0]) );
      if ( !m_catRef.empty() )  {
        text="         "+m_catRef; /* at least 8 leading blanks */
        header.append(text);
      }
    }
    if ( !m_tableName.empty() ) {
      header.append( KeyRecord(KeyCDS[1],m_tableName,KeyCom[1]) );
      if ( !m_tableRef.empty() ) {
        text="         "+m_tableRef; /* at least 8 leading blanks */
        header.append(text);
      }
    }
  }
  catch (const TipException &x) {
    text="fits EXTENSION, cannot add CDS header keys";
    printWarn(origin, text);
  }
  try {
    unsigned int pos;
    IColumn *myCol = 0;

    err=m_quantities.size();
    for (j=0; j < err; j++) {
      const Quantity &readQ=m_quantities[j];
      text=readQ.m_format; // from binary: same format
      if ( !m_URL.empty() ) {
        /* from ASCII: Iw or Fw.d or Ew.d or Dw.d (Aw for string) */
        if (readQ.m_type == Quantity::NUM) {
//          if (text[0]=='I') text="1J"; else text="1D";
          text="1D"; // to keep NaN values, improve with TNULL for integers
        }
        else if (readQ.m_type == Quantity::STRING) {
          text.erase(0,1);
          text=text+"A";
        }
      }
      (*ptrTable)->appendField(readQ.m_name, text);
      myCol=(*ptrTable)->getColumn(j);
      text=readQ.m_null;
      if ( !text.empty() ) {
        if ( text[0]!=SepNull ) {
          myCol->getColumnKeyword("TNULL").set( atoi(text.c_str()) );
          myCol->getColumnKeyword("TNULL").setComment("Undefined value of field");
        }
        pos=text.find(SepNull); /* always exist */
        text.erase(0, pos+1);
        pos=text.find(SepNull); /* always exist */
        if (pos > 0) {
          myCol->getColumnKeyword("TSCAL").set( atoi(text.c_str()) );
          myCol->getColumnKeyword("TSCAL").setComment("to rescale data");
        }
        text.erase(0, pos+1);
        if ( !text.empty() ) {
          myCol->getColumnKeyword("TZERO").set( atoi(text.c_str()) );
          myCol->getColumnKeyword("TZERO").setComment("offset for unsigned integers");
        }
      }
      text=readQ.m_unit;
      if ( !text.empty() ) {
        myCol->getColumnKeyword("TUNIT").set(readQ.m_unit);
        myCol->getColumnKeyword("TUNIT").setComment("physical unit of field");
      }
      myCol->getColumnKeyword(Key_UCD).set(readQ.m_ucd);
      myCol->getColumnKeyword(Key_UCD).setComment(readQ.m_comment);
    }/* loop on quantities */
  }
  catch (const TipException &x) {
    text=": fits EXTENSION, cannot add column";
    printErr(origin, text);
    return BAD_FITS;
  }
/*  try {
    for (Header::Iterator itor=header.begin(); itor!=header.end(); ++itor) {
printf("%s\n", itor->getName().c_str() );
      if ( (!m_catRef.empty()) && (itor->getName()==KeyCDS[1]) ) {
        text="         "+m_catRef;
        header.insert(itor, text);
      }
    }
  }
  catch (const TipException &x) {
    text="fits EXTENSION, cannot add CDS header long comment";
    printWarn(origin, text);
  }
*/
  return IS_OK;
}
/**********************************************************************/
// save whole catalog from memory to a FITS file
int Catalog::saveFits(const std::string &fileName, const std::string &extName,
                      bool clobber, bool append) {

  const std::string origin="saveFits";
  Table *myDOL=0;
  int err;
  err=createFits(fileName, extName, clobber, append, origin, &myDOL);
  if (err < IS_VOID) return err;

  std::ostringstream sortie;
  long tot=0l;
  int  i, j;

  if (m_numRows > 0) {
  try {
    myDOL->setNumRecords(m_numRows);
    Quantity readQ;
    double rowVal;
    err=m_quantities.size();
    std::vector<long> colNull(err, 0);
    std::vector<bool> no_Null(err, true);
    for (j=0; j < err; j++) {
      readQ=m_quantities[j];
      if ( !readQ.m_null.empty() ) {
        colNull[j]=atoi( readQ.m_null.c_str() );
        if (readQ.m_null[0]!=SepNull) no_Null[j]=false;
      }
    }
    // Loop over all records (rows) and set values
    for (Table::Iterator itor=myDOL->begin(); itor != myDOL->end(); ++itor) {
      // double variable to hold the value of all the numeric fields
      for (j=0; j < err; j++) {
        readQ=m_quantities[j];
        i=readQ.m_index;
        if (readQ.m_type == Quantity::NUM) {
          rowVal= m_numericals[i].at(tot);
          if ( no_Null[j] ) (*itor)[readQ.m_name].set(rowVal);
          else {
#ifdef WIN32
            if ( _isnan(rowVal) ) (*itor)[readQ.m_name].set(colNull[j]);
#else
            if ( std::isnan(rowVal) ) (*itor)[readQ.m_name].set(colNull[j]);
#endif
            else (*itor)[readQ.m_name].set(rowVal);
          }
        }
        else if (readQ.m_type == Quantity::STRING)
          (*itor)[readQ.m_name].set( m_strings[i].at(tot) );
      }
      tot++;
    }/* loop on rows */
  }
  catch (const TipException &x) {
    sortie << ": fits EXTENSION, cannot write at row#" << tot;
    printErr(origin, sortie.str() );
    return BAD_ROW;
  }/* end of try */
  }/* at least 1 row */

  delete myDOL; myDOL=0;
  sortie << "output fits is closed ( " << tot << " rows written)";
  printLog(0, sortie.str());
  return IS_OK;
}
/**********************************************************************/
// save selected rows of catalog from memory to a FITS file
int Catalog::saveSelectedFits(const std::string &fileName,
                              const std::string &extName,
                              bool clobber, bool append) {

  const std::string origin="saveSelectedFits";
  Table *myDOL=0;
  int err;
  err=createFits(fileName, extName, clobber, append, origin, &myDOL);
  if (err < IS_VOID) return err;

  std::ostringstream sortie;
  long tot=0l;
  int  i, j;

  if (m_numSelRows > 0) {
  try {
    myDOL->setNumRecords(m_numSelRows);
    Quantity readQ;
    double rowVal;
    err=m_quantities.size();
    std::vector<long> colNull(err, 0);
    std::vector<bool> no_Null(err, true);
    for (j=0; j < err; j++) {
      readQ=m_quantities[j];
      if ( !readQ.m_null.empty() ) {
        colNull[j]=atoi( readQ.m_null.c_str() );
        if (readQ.m_null[0]!=SepNull) no_Null[j]=false;
      }
    }
    Table::Iterator itor=myDOL->begin();
    // Loop over all selected records (rows) and set values
    // first bit indicates global selection
    for (long k=0; k<m_numRows; k++) if (m_rowIsSelected[0].at(k) & 1) {
      for (j=0; j < err; j++) {
        readQ=m_quantities[j];
        i=readQ.m_index;
        if (readQ.m_type == Quantity::NUM) {
          rowVal= m_numericals[i].at(tot);
          if ( no_Null[j] ) (*itor)[readQ.m_name].set(rowVal);
          else {
#ifdef WIN32
            if ( _isnan(rowVal) ) (*itor)[readQ.m_name].set(colNull[j]);
#else
            if ( std::isnan(rowVal) ) (*itor)[readQ.m_name].set(colNull[j]);
#endif
            else (*itor)[readQ.m_name].set(rowVal);
            // to be improved for integers
          }
        }
        else if (readQ.m_type == Quantity::STRING)
          (*itor)[readQ.m_name].set( m_strings[i].at(tot) );
      }
      tot++;
      itor++;
      if (tot == m_numSelRows) break;
    }
  }
  catch (const TipException &x) {
    sortie << ": fits EXTENSION, cannot write at row#" << tot;
    printErr(origin, sortie.str() );
    return BAD_ROW;
  }/* end of try */
  }/* at least 1 selected row */

  delete myDOL; myDOL=0;
  sortie << "output fits is closed ( " << tot << " rows written)";
  printLog(0, sortie.str());
  return IS_OK;
}

} // namespace catalogAccess
