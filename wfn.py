#!/usr/bin/python

import numpy

def load_wfn(filename):

    def helper_num(f):
        '''Reads number of orbitals, primitives and atoms'''
        line = f.readline()
        assert line.startswith('GAUSSIAN')
        return [int(i) for i in line.split() if i.isdigit()]

    def helper_coordinates(f):
        '''Reads the coordiantes of the atoms'''
        numbers = numpy.empty(num_atoms, dtype=numpy.int32)
        coordinates = numpy.empty((num_atoms, 3))
        symbols = []
        for atom in range(num_atoms):
            line = f.readline()
            line = line.split()
            symbols.append(line[0])
            numbers[atom] = int(float(line[9]))
            coordinates[atom, :] = [line[4], line[5], line[6]]
        symbols = numpy.asarray(symbols, dtype=numpy.character)
        return symbols, numbers, coordinates

    def helper_section(f, start, skip):
        '''Reads CENTRE ASSIGNMENTS, TYPE ASSIGNMENTS, and EXPONENTS sections'''
        section = []
        while len(section) < num_primitives:
            line = f.readline()
            assert line.startswith(start)
            line = line.split()
            section.extend(line[skip:])
        assert len(section) == num_primitives
        return section

    def helper_mo(f):
        '''Reads all MO information'''
        line = f.readline()
        assert line.startswith('MO')
        line = line.split()
        count = line[1]
        occ, energy = line[-5], line[-1]
        coeffs = helper_section(f, ' ', 0)
        coeffs = [i.replace('D', 'E') for i in coeffs]
        return count, occ, energy, coeffs

    with open(filename) as f:
        title = f.readline().strip()
        num_mo, num_primitives, num_atoms = helper_num(f)
        symbols, numbers, coordinates = helper_coordinates(f)
        centers = helper_section(f, 'CENTRE ASSIGNMENTS', 2)
        centers = numpy.array([int(i) for i in centers], dtype=numpy.int32)
        type_assignment = helper_section(f, 'TYPE ASSIGNMENTS', 2)
        type_assignment = numpy.array([int(i) for i in type_assignment],
                                      dtype=numpy.int32)
        exponents = helper_section(f, 'EXPONENTS', 1)
        exponents = numpy.array([float(i.replace('D', 'E')) for i in exponents])
        mo_count = numpy.empty(num_mo, dtype=numpy.int32)
        mo_occ = numpy.empty(num_mo)
        mo_energy = numpy.empty(num_mo)
        coefficients = numpy.empty([num_primitives, num_mo])
        for mo in range(num_mo):
            mo_count[mo], mo_occ[mo], mo_energy[
                mo], coefficients[:, mo] = helper_mo(f)
        line = f.readline()

    return title, symbols, numbers, coordinates, centers, type_assignment, exponents, \
        mo_count, mo_occ, mo_energy, coefficients

