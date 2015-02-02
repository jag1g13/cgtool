#ifndef ITP_WRITER_H_
#define ITP_WRITER_H_

#include <string>
#include <vector>

#include <stdio.h>

#include "cg_map.h"

/** \brief Contains all functions necessary to print mapping to a GROMACS ITP file
* Handles the [moleculetype], [atoms], [bonds], [angles] and [dihedrals] sections
*/
class ITPWriter{
protected:
    /** Pointer to the ITP file being written to */
    FILE *itp_;
    /** The name of the ITP file */
    std::string name_;
    /** Which section of the file are we writing at the moment? */
    std::string section_;
    /** Header printed at the top of every ITP file */
    const std::string header_ =
            ";\n"
            "; Topology prepared automatically using CGTOOL\n"
            "; James Graham <J.A.Graham@soton.ac.uk> 2015\n"
            "; University of Southampton\n"
            ";\n; This isn't quite a valid GROMACS topology yet\n"
            ";\n";

public:
    /** Create an ITP file and prepare to write */
    ITPWriter(std::string name);

    /** Create a new section in the ITP file */
    void newSection(std::string section_name);

    /** Print atoms to itp */
//    void printAtoms(const std::vector<BeadMap> &mapping);
    void printAtoms(const CGMap &map);

    /** Print bond params to itp */
    void printBonds(const BondSet &bond_set);
};

#endif