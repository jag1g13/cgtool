//
// Created by james on 31/03/15.
//

#ifndef CGTOOL_FILE_IO_H
#define CGTOOL_FILE_IO_H

#include <map>

enum class FileFormat{GROMACS, LAMMPS};
enum class FieldFormat{MARTINI, ELBA, OTHER};
enum class PotentialType{HARMONIC, COS, COSSQUARED};

const std::map<std::string, FileFormat> getFileFormat =
        {{"GROMACS", FileFormat::GROMACS},
         {"LAMMPS",  FileFormat::LAMMPS }};

const std::map<std::string, FieldFormat> getFieldFormat =
        {{"MARTINI", FieldFormat::MARTINI},
         {"ELBA",    FieldFormat::ELBA   },
         {"OTHER",   FieldFormat::OTHER  }};

const std::map<std::string, PotentialType> getPotential =
        {{"HARMONIC",   PotentialType::HARMONIC},
         {"COS",        PotentialType::COS},
         {"COSSQUARED", PotentialType::COSSQUARED}};

#endif //CGTOOL_FILE_IO_H
