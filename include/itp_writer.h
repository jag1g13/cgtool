#ifndef ITP_WRITER_H_
#define ITP_WRITER_H_

#include <string>
#include <vector>

#include <stdio.h>

#include "cg_map.h"
#include "frame.h"

/** \brief Contains all functions necessary to print mapping to a GROMACS ITP file
* Handles the [moleculetype], [atoms], [bonds], [angles] and [dihedrals] sections
*/
class ITPWriter{
protected:
    /** Pointer to the ITP file being written to */
    FILE *itp_;
    /** The name of the ITP file */
    std::string name_;
    /** Residue name */
    std::string resName_;
    /** Which section of the file are we writing at the moment? */
    std::string section_;
    /** Header printed at the top of every ITP file */
    const std::string header_ =
            ";\n"
            "; Topology prepared automatically using CGTOOL\n"
            "; James Graham <J.A.Graham@soton.ac.uk> 2015\n"
            "; University of Southampton\n"
            "; https://bitbucket.org/jag1g13/cgtool\n"
            ";\n";

public:
    /** Create an ITP file and prepare to write */
    ITPWriter(const std::string &resName, std::string itpname="");

    /** Close output file in destructor */
    ~ITPWriter();

    /** Create a new section in the ITP file */
    void newSection(const std::string &section_name);

    /** \brief Print atoms section of itp.
    * If isMartini is true, don't print a masses column
    * and only put charges on 'charged/Q' beads. */
    void printAtoms(const CGMap &map, const bool isMartini=true);

    /** Print bond params to itp */
    void printBonds(const BondSet &bond_set);
};

#endif