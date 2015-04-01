#ifndef ITP_WRITER_H_
#define ITP_WRITER_H_

#include <string>
#include <vector>

#include <stdio.h>

#include "cg_map.h"
#include "frame.h"
#include "file_io.h"

/** \brief Contains all functions necessary to print mapping to a
 * GROMACS/LAMMPS force field file.  Handles the [moleculetype],
 * [atoms], [bonds], [angles] and [dihedrals] sections
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
    /** Which output format should be written? */
    FileFormat format_;
    /** Header printed at the top of every ITP file */
    const std::vector<std::string> header_ =
        {
            "Topology prepared automatically using CGTOOL\n",
            "James Graham <J.A.Graham@soton.ac.uk> 2015\n",
            "University of Southampton\n",
            "https://bitbucket.org/jag1g13/cgtool\n",
        };

    /** Create a new section in the ITP file */
    void newSection(const std::string &section_name);

public:
    /** Create an ITP file and prepare to write */
    ITPWriter(const std::string &resName, FileFormat format=FileFormat::GROMACS, std::string itpname="");

    /** Close output file in destructor */
    ~ITPWriter();

    /** \brief Print atoms section of itp.
    * If isMartini is true, don't print a masses column
    * and only put charges on 'charged/Q' beads.
    */
    void printAtoms(const CGMap &map, const bool isMartini=true);

    /** Print bond params to itp */
    void printBonds(const BondSet &bond_set, const bool round=false);
};

#endif