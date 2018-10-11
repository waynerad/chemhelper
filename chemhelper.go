package main

import (
	"fmt"
	"strconv"
)

const avogadrosNumber = 6.0221420e23
const cm3perL = 1000

type element struct {
	name               string
	symbol             string
	atomicNumber       int
	relativeAtomicMass float64
	radioactive        bool
}

type isotope struct {
	atomicNumber int
	massNumber   int
	isotopicMass float64
	abundance    float64
}

// A formele is a "formula element".
// A formula like "2Hg + O2" has 2 formula elements, "2Hg" and "O2".
// Formulas don't take into account molecular divisions.
// They are used to obtain relative molecular weights and to make sure
// equations balance.
// If you need molecular distinctions, use molele/molecule.
type formele struct {
	coef          int
	elementNumber int
	numberOfAtoms int
}

type formula struct {
	components []formele
}

type iChemHelper interface {
	init()

	nameForSym(sbl string) string
	atomicNumberSym(sbl string) int
	atomicMassSym(sbl string) float64
	radioactiveSym(sbl string) bool

	makeFormele(coef string, elementSymbol string, atomCount string) formele
	parseFormula(fmla string) formula
	countFormulaAtoms(fmla formula) map[int]int
	determineMolecularMass(atomCounts map[int]int) float64
}

type chemHelper struct {
	elem         []element
	atomicNumber []int
	symbol       map[string]int
	isos         []isotope
	isolookup    map[int]map[int]int
}

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
		{6, 13, 13.003354, .01108}}
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

func isRuneSpace(ch rune) bool {
	return ch == 32
}

func isRunePlus(ch rune) bool {
	return ch == 43
}

func isRuneNumeric(ch rune) bool {
	return (ch >= 48) && (ch <= 58)
}

func isRuneLowercase(ch rune) bool {
	return ch > 90
}

// generic error-checker -- substitute with particular error checker when
// appropriate
func checkError(err error, message string) {
	if err == nil {
		return
	}
	panic(message)
}

func strToInt(str string) int {
	result, err := strconv.ParseInt(str, 10, 32)
	checkError(err, "cannot parse as int")
	return int(result)
}

// ------------------------------------------------------------------------
// Begin chemHelper interface implementation
// ------------------------------------------------------------------------

func (self *chemHelper) init() {
	self.elem = constructElementRelativeAtomicMassesTable()
	self.atomicNumber, self.symbol = makeRelativeAtomicMassesTableIndexes(self.elem)
	self.isos = constructIsotopeRelativeMassesTable()
	self.isolookup = makeIsotopicMassesTableIndexes(self.isos)
}

func (self *chemHelper) nameForSym(sbl string) string {
	return self.elem[self.symbol[sbl]].name
}

func (self *chemHelper) atomicNumberSym(sbl string) int {
	return self.elem[self.symbol[sbl]].atomicNumber
}

func (self *chemHelper) atomicMassSym(sbl string) float64 {
	return self.elem[self.symbol[sbl]].relativeAtomicMass
}

func (self *chemHelper) radioactiveSym(sbl string) bool {
	return self.elem[self.symbol[sbl]].radioactive
}

func (self *chemHelper) makeFormele(coef string, elementSymbol string, atomCount string) formele {
	if coef == "" {
		coef = "1"
	}
	if atomCount == "" {
		atomCount = "1"
	}
	var result formele
	result.coef = strToInt(coef)
	result.elementNumber = self.elem[self.symbol[elementSymbol]].atomicNumber
	result.numberOfAtoms = strToInt(atomCount)
	return result
}
func (self *chemHelper) parseFormula(fmla string) formula {
	var result formula
	// formulas look like: 2Hg + O2
	inCoef := true
	coefStr := ""
	symbolStr := ""
	countStr := ""
	for _, char := range fmla {
		if isRuneSpace(char) || isRunePlus(char) {
			if isRunePlus(char) {
				fmele := self.makeFormele(coefStr, symbolStr, countStr)
				result.components = append(result.components, fmele)
				inCoef = true
				coefStr = ""
				symbolStr = ""
			}
		} else {
			if isRuneNumeric(char) {
				if inCoef {
					coefStr += string(char)
				} else {
					countStr += string(char)
				}
			} else {
				if isRuneLowercase(char) {
					symbolStr += string(char)
				} else {
					if symbolStr != "" {
						// we need to push a previous symbol
						fmele := self.makeFormele(coefStr, symbolStr, countStr)
						result.components = append(result.components, fmele)
					}
					symbolStr = string(char)
					inCoef = false
					countStr = ""
				}
			}
		}
	}
	// handle final symbol at end
	if symbolStr != "" {
		fmele := self.makeFormele(coefStr, symbolStr, countStr)
		result.components = append(result.components, fmele)
	}
	return result
}
func (self *chemHelper) countFormulaAtoms(fmla formula) map[int]int {
	result := make(map[int]int)
	for _, fmele := range fmla.components {
		_, ok := result[fmele.elementNumber]
		if ok {
			result[fmele.elementNumber] += fmele.coef * fmele.numberOfAtoms
		} else {
			result[fmele.elementNumber] = fmele.coef * fmele.numberOfAtoms
		}
	}
	return result
}
func (self *chemHelper) determineMolecularMass(atomCounts map[int]int) float64 {
	var result float64
	result = 0.0
	for elemtNum, numberOfAtoms := range atomCounts {
		result += self.elem[self.atomicNumber[elemtNum]].relativeAtomicMass * float64(numberOfAtoms)
	}
	return result
}

// ------------------------------------------------------------------------
// End chemHelper interface implementation
// ------------------------------------------------------------------------

func testMain() {
	elem := constructElementRelativeAtomicMassesTable()
	testCheckForDuplicates(elem)
	atomicNumber, symbol := makeRelativeAtomicMassesTableIndexes(elem)
	for num := 0; num < 115; num++ {
		fmt.Println(num, elem[atomicNumber[num]])
	}
	fmt.Println("symbol", symbol)
	var chemh chemHelper
	chemh.init()
	fmt.Println("chemhelp", chemh)
	fmt.Println("H", chemh.nameForSym("H"))
	fmt.Println("O", chemh.atomicNumberSym("O"))
	fmt.Println("C", chemh.atomicMassSym("C"))
	fmt.Println("Si", chemh.radioactiveSym("Si"))
	testFormula := chemh.parseFormula("2Hg + O2")
	fmt.Println("2Hg + O2 formula", testFormula)
	testFormula = chemh.parseFormula("7H2SO4 + 2Hg + O2")
	fmt.Println("7H2SO4 + 2Hg + O2 formula", testFormula)
	atomCounts := chemh.countFormulaAtoms(testFormula)
	fmt.Println("atomCounts", atomCounts)
	molarMass := chemh.determineMolecularMass(atomCounts)
	fmt.Println("molarMass", molarMass)
}

func main() {
	var chemh chemHelper
	chemh.init()

	massOfOCompoundA := 0.22564
	fmt.Println("E1.1 Law of multiple proportions: B", 0.90255/massOfOCompoundA, "C", 1.3539/massOfOCompoundA, "D", 1.5795/massOfOCompoundA)

	fmt.Println("E1.2 Number of electrons, protons, and neutrons in 222Rn (Radon-222):", chemh.atomicNumberSym("Rn"), chemh.atomicNumberSym("Rn"), 222-chemh.atomicNumberSym("Rn"))

	// fmt.Println("E1.3 Relative atomic mass of carbon (calculated from isotopes)", (isos[isolookup[chemh.xyzSym("C"]].atomicNumber][12]].isotopicMass*isos[isolookup[chemh.xyzSym("C"]].atomicNumber][12]].abundance)+(isos[isolookup[chemh.xyzSym("C"]].atomicNumber][13]].isotopicMass*isos[isolookup[chemh.xyzSym("C"]].atomicNumber][13]].abundance))

	fmt.Println("E1.3 Relative atomic mass of carbon (from table)", chemh.atomicMassSym("C"))

	fmt.Println("Relative molecular mass of water:", (chemh.atomicMassSym("H")*2)+chemh.atomicMassSym("O"))

	fmt.Println("Mass of a single 12C atom:", chemh.atomicMassSym("C")/avogadrosNumber, "g")

	fmt.Println("Ratio of mass of sodium atom (Na) to carbon atom (C)", chemh.atomicMassSym("Na")/chemh.atomicMassSym("C"))

	fmt.Println("E1.4 Mass of one uranium atom:", chemh.atomicMassSym("U")/avogadrosNumber, "g")

	fmt.Println("Chemical amount of iron in 8.232g:", 8.232/chemh.atomicMassSym("Fe"), "mol")

	molarMassOfWater := (chemh.atomicMassSym("H") * 2) + chemh.atomicMassSym("O")
	fmt.Println("Amount of water needed for 0.2000 mol:", 0.2*molarMassOfWater, "g")

	molarMassOfNO2 := chemh.atomicMassSym("N") + (chemh.atomicMassSym("O") * 2)
	fmt.Println("E1.5a chemical amount of NO2 in 4.00 g of NO2:", 4.0/molarMassOfNO2, "mol")
	fmt.Println("E1.5b number of molecules of NO2:", (4.0/molarMassOfNO2)*avogadrosNumber)

	liquidBenzeneDensity := 0.8765 // in g/cm^3
	fmt.Println("Liquid benzene mass:", (0.2124*cm3perL)*liquidBenzeneDensity, "g")
	molecularMassBenzene := (chemh.atomicMassSym("C") * 6) + (chemh.atomicMassSym("H") * 6)
	fmt.Println("Molecular mass of benzene:", molecularMassBenzene)
	fmt.Println("Chemical amount of benzene:", ((0.2124*cm3perL)*liquidBenzeneDensity)/molecularMassBenzene, "mol")

	iceDensityNear0C := 0.92 // in g/cm^3
	fmt.Println("Molar volume of solid water (ice):", molarMassOfWater/iceDensityNear0C, "cm^3/mol")

	oxygenDensityRoomTempSeaLevel := 0.00130 // in g/cm^3
	fmt.Println("Molar volume of oxygen (O2) at room temperature:", ((chemh.atomicMassSym("O")*2)/oxygenDensityRoomTempSeaLevel)/cm3perL, "L/mol")

	fmt.Println("Volume per H2O molecule for solid water (ice):", (molarMassOfWater/avogadrosNumber)/iceDensityNear0C)

	fmt.Println("Molecular mass of water from formula:", chemh.determineMolecularMass(chemh.countFormulaAtoms(chemh.parseFormula("H2O"))))
	fmt.Println("Molecular mass of benzene from formula:", chemh.determineMolecularMass(chemh.countFormulaAtoms(chemh.parseFormula("C6H6"))))

}
