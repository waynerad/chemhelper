package main

import (
	"fmt"
)

type element struct {
	name               string
	symbol             string
	atomicNumber       int
	relativeAtomicMass float64
	radioactive        bool
}

type isotope struct {
	atomicNumber int
	massNumber int
	isotopicMass float64
	abundance float64
}

type molele struct {
	elementAtomicNumber int
	numberOfAtoms int
}

type molecule struct {
	components []molele
}

const avogadrosNumber = 6.0221420e23

func constructElementRelativeAtomicMassesTable() []element {
	relativeAtomicMassesTable := []element{
		{"Actinium", "Ac", 89, 227.0277, true},
		{"Aluminum", "Al", 13, 26.981538, false},
		{"Americium", "Am", 95, 243.0614, true},
		{"Antimony", "Sb", 51, 121.760, false},
		{"Argon", "Ar", 18, 39.948, false},
		{"Arsenic", "As", 33, 74.92160, false},
		{"Astatine", "At", 85, 209.9871, true},
		{"Barium", "Ba", 56, 137.327, false},
		{"Berkelium", "Bk", 97, 247.0703, true},
		{"Beryllium", "Be", 4, 9.012182, false},
		{"Bismuth", "Bi", 83, 208.98038, false},
		{"Bohrium", "Bh", 107, 264.12, false},
		{"Boron", "B", 5, 10.811, false},
		{"Bromine", "Br", 35, 79.904, false},
		{"Cadmium", "Cd", 48, 112.411, false},
		{"Calcium", "Ca", 20, 40.078, false},
		{"Californium", "Cf", 98, 251.0796, true},
		{"Carbon", "C", 6, 12.0107, false},
		{"Cerium", "Ce", 58, 140.116, false},
		{"Cesium", "Cs", 55, 132.90545, false},
		{"Chlorine", "Cl", 17, 35.4527, false},
		{"Chromium", "Cr", 24, 51.9961, false},
		{"Cobalt", "Co", 27, 58.99320, false},
		{"Copper", "Cu", 29, 63.456, false},
		{"Curium", "Cm", 96, 247.0703, true},
		{"Dubnium", "Db", 105, 262.1144, false},
		{"Dysprosium", "Dy", 66, 162.50, false},
		{"Einsteinium", "Es", 99, 252.0830, true},
		{"Erbium", "Er", 68, 167.26, false},
		{"Europium", "Eu", 63, 151.964, false},
		{"Fermium", "Fm", 100, 257.0951, true},
		{"Fluorine", "F", 9, 18.9984032, false},
		{"Francium", "Fr", 87, 223.0197, true},
		{"Gadolinium", "Gd", 64, 157.25, false},
		{"Gallium", "Ga", 31, 69.723, false},
		{"Germanium", "Ge", 32, 72.64, false},
		{"Gold", "Au", 79, 196.96655, false},
		{"Hafnium", "Hf", 72, 178.49, false},
		{"Hassium", "Hs", 108, 277.0, true},
		{"Helium", "He", 2, 4.002602, false},
		{"Holmium", "Ho", 67, 164.93032, false},
		{"Hydrogen", "H", 1, 1.00794, false},
		{"Indium", "In", 49, 114.818, false},
		{"Iodine", "I", 53, 126.90447, false},
		{"Iridium", "Ir", 77, 192.217, false},
		{"Iron", "Fe", 26, 55.845, false},
		{"Krypton", "Kr", 36, 83.798, false},
		{"Lanthanum", "La", 57, 138.9055, false},
		{"Lawrencium", "Lr", 103, 262.1097, true},
		{"Lead", "Pb", 82, 207.2, false},
		{"Lithium", "Li", 3, 6.941, false},
		{"Lutetium", "Lu", 71, 174.967, false},
		{"Magnesium", "Mg", 12, 24.3050, false},
		{"Manganese", "Mn", 25, 54.938049, false},
		{"Meitnerium", "Mt", 109, 268.1388, true},
		{"Mendelevium", "Md", 101, 258.0984, true},
		{"Mercury", "Hg", 80, 200.59, false},
		{"Molybdenum", "Mo", 42, 95.94, false},
		{"Neudymium", "Nd", 60, 144.24, false},
		{"Neon", "Ne", 10, 20.1797, false},
		{"Neptunium", "Np", 93, 237.0482, true},
		{"Nickel", "Ni", 28, 58.6934, false},
		{"Niobium", "Nb", 41, 92.90618, false},
		{"Nitrogen", "N", 7, 14.00674, false},
		{"Nobelium", "No", 102, 259.1011, true},
		{"Osmium", "Os", 76, 190.23, false},
		{"Oxygen", "O", 8, 15.9994, false},
		{"Palladium", "Pd", 46, 106.42, false},
		{"Phosphorus", "P", 15, 30.973761, false},
		{"Platnum", "Pt", 78, 195.078, false},
		{"Plutonium", "Pu", 94, 244.0642, true},
		{"Polonium", "Po", 84, 208.9824, true},
		{"Potassium", "K", 19, 39.0983, false},
		{"Praseodymium", "Pr", 59, 140.90765, false},
		{"Promethium", "Pm", 61, 144.9127, true},
		{"Protactinium", "Pa", 91, 231.0359, false},
		{"Radium", "Ra", 88, 226.0254, true},
		{"Radon", "Rn", 86, 222.0176, true},
		{"Rhenium", "Re", 75, 186.207, false},
		{"Rhodium", "Rh", 45, 102.90550, false},
		{"Rubidium", "Rb", 37, 85.4678, false},
		{"Ruthenium", "Ru", 44, 101.07, false},
		{"Rutherfordium", "Rf", 104, 261.1089, true},
		{"Samarium", "Sm", 62, 150.36, false},
		{"Scandium", "Sc", 21, 44.955910, false},
		{"Seaborgium", "Sg", 106, 266.1219, true},
		{"Selenium", "Se", 34, 78.96, false},
		{"Silicon", "Si", 14, 28.0855, false},
		{"Silver", "Ag", 47, 107.8682, false},
		{"Sodium", "Na", 11, 22.989770, false},
		{"Strontium", "Sr", 38, 87.82, false},
		{"Sulfur", "S", 16, 32.065, false},
		{"Tantalum", "Ta", 73, 180.9479, false},
		{"Technetium", "Tc", 43, 97.9072, true},
		{"Tellurium", "Te", 52, 127.60, false},
		{"Terbium", "Tb", 65, 158.92534, false},
		{"Thallium", "Tl", 81, 204.3833, false},
		{"Thorium", "Th", 90, 232.0381, false},
		{"Thulium", "Tm", 69, 168.93421, false},
		{"Tin", "Sn", 50, 118.710, false},
		{"Titanium", "Ti", 22, 47.867, false},
		{"Tungsten", "W", 74, 183.84, false},
		{"Copernicium", "Cn", 112, 285.0, true},
		{"Darmstadtium", "Ds", 110, 281.0, true},
		{"Flerovium", "Fl", 114, 289.0, true},
		{"Roentgenium", "Rg", 111, 272, true},
		{"Uranium", "U", 92, 238.0289, false},
		{"Vanadium", "V", 23, 50.9415, false},
		{"Xenon", "Xe", 54, 131.29, false},
		{"Ytterbium", "Yb", 70, 173.04, false},
		{"Yttrium", "Y", 39, 88.90585, false},
		{"Zinc", "Zn", 30, 65.409, false},
		{"Zirconium", "Zr", 40, 91.224, false},
		{"Nihonium", "Nh", 113, 286.0, true}}

	return relativeAtomicMassesTable
}

func constructIsotopeRelativeMassesTable() []isotope {
	relativeAtomicMassesTable := []isotope{
		{6, 12, 12.0, 0.98892},
		{6, 13, 13.003354, .01108} }
	return relativeAtomicMassesTable
}

func testCheckForDuplicates(relativeAtomicMassesTable []element) {
	dupName := make(map[string]bool)
	dupSymbol := make(map[string]bool)
	dupAtomicNumber := make(map[int]bool)
	dupRelativeAtomicMass := make(map[float64]bool)
	for _, elemtInfo := range relativeAtomicMassesTable {
		_, already := dupName[elemtInfo.name]
		if already {
			fmt.Println("Duplicate name:", elemtInfo.name)
		}
		_, already = dupSymbol[elemtInfo.symbol]
		if already {
			fmt.Println("Duplicate symbol:", elemtInfo.symbol)
		}
		_, already = dupAtomicNumber[elemtInfo.atomicNumber]
		if already {
			fmt.Println("Duplicate atomic number:", elemtInfo.atomicNumber)
		}
		_, already = dupRelativeAtomicMass[elemtInfo.relativeAtomicMass]
		if already {
			fmt.Println("Duplicate relative atomic mass:", elemtInfo.relativeAtomicMass)
		}
	}
}

func makeRelativeAtomicMassesTableIndexes(relativeAtomicMassesTable []element) ([]int, map[string]int) {
	atomicNumber := make([]int, 115)
	symbol := make(map[string]int)
	for idx, elemtInfo := range relativeAtomicMassesTable {
		atomicNumber[elemtInfo.atomicNumber] = idx
		symbol[elemtInfo.symbol] = idx
	}
	return atomicNumber, symbol
}

func makeIsotopicMassesTableIndexes(isos []isotope) map[int]map[int]int {
	isolookup := make(map[int]map[int]int)
	for idx, isoInfo := range isos {
		_, exists := isolookup[isoInfo.atomicNumber]
		if !exists {
		isolookup[isoInfo.atomicNumber] = make(map[int]int)
		}
		isolookup[isoInfo.atomicNumber][isoInfo.massNumber] = idx
	}
	return isolookup
}

func testMain() {
	elem := constructElementRelativeAtomicMassesTable()
	testCheckForDuplicates(elem)
	atomicNumber, symbol := makeRelativeAtomicMassesTableIndexes(elem)
	for num := 0; num < 115; num++ {
		fmt.Println(num, elem[atomicNumber[num]])
	}
	fmt.Println("H", elem[symbol["H"]])
	fmt.Println("O", elem[symbol["O"]])
	fmt.Println("C", elem[symbol["C"]])
	fmt.Println("Si", elem[symbol["Si"]])
}

func main() {
	elem := constructElementRelativeAtomicMassesTable()
	atomicNumber, symbol := makeRelativeAtomicMassesTableIndexes(elem)
	isos := constructIsotopeRelativeMassesTable()
	// fmt.Println("iso",isos)

	isolookup := makeIsotopicMassesTableIndexes(isos) 

	fmt.Println("isolookup",isolookup)

	fmt.Println(elem[atomicNumber[1]])

	fmt.Println("Relative atomic mass of carbon (calculated from isotopes)", (isos[isolookup[elem[symbol["C"]].atomicNumber][12]].isotopicMass * isos[isolookup[elem[symbol["C"]].atomicNumber][12]].abundance) + (isos[isolookup[elem[symbol["C"]].atomicNumber][13]].isotopicMass * isos[isolookup[elem[symbol["C"]].atomicNumber][13]].abundance) )

	fmt.Println("Relative atomic mass of carbon (from table)", elem[symbol["C"]].relativeAtomicMass)

	fmt.Println("Relative molecular mass of water:", (2*elem[symbol["H"]].relativeAtomicMass)+elem[symbol["O"]].relativeAtomicMass)

	fmt.Println("Mass of one uranium atom:", elem[symbol["U"]].relativeAtomicMass/avogadrosNumber)
}
