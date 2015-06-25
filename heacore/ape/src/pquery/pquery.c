/** \file pquery.c
    \brief
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_binary.h"
#include "ape/ape_error.h"

int main(int argc, char ** argv) {
  return eOK != ape_binary_pquery(argc - 1, argv + 1) ? 1 : 0;
}

/*
 * $Log: pquery.c,v $
 * Revision 1.1  2012/04/16 18:53:25  irby
 * New APE-based pquery binary (to replace the old XPI version in HEASoft).
 *
 */
