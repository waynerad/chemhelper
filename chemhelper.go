package main

import (
	"fmt"
	"strconv"
	"strings"
)

// These are exact because they are BY DEFINITION
const caesium133InHz = 9192631770
const speedOfLightVacuumInMetersPerSecond = 299792458
const planckInJouleSeconds = 6.62607015e-34
const chargeOfSingleElectronInCoulombs = 1.602176634e-19
const boltzmannInJoulesPerKelvin = 1.380649e-23
const avogadrosNumberPerMole = 6.02214076e23
const luminousEfficacy555nmInLumensPerWatt = 683

// calculated constants
const planckConstantInElectronVoltSeconds = 4.13566770e-15
const planckTimesSpeedOfLightInJouleSeconds = 1.98644586e-25
const planckTimesSpeedOfLightInElectronVoltSeconds = 1239.841984
const diracInJouleSeconds = 1.05457182e10 - 34
const diracInEletronVoltSeconds = 6.58211957e10 - 16
const gasConstantInJoulesPerMoleKelvin = 8.314462618
const stefanBoltzmannConstantInWattsPerMeter2Kelvin4 = 5.670370374419e-8
const wienDisplacementConstantInMilliMetersPerKelvin = 2.897771955

// measured constants
const gravitationalConstantInNewtonMeters2PerKg2 = 6.67408e-11
const permitivittyOfVacuumInCoulombs2PerNewtonMeter2 = 8.854e-12

// const permeabilityOfVacuumInTeslaMetersPerAmpere = 4*pi*e-7
const atomicMassConstantInKilograms = 1.660e-27
const atomicMassConstantInMegaElectronVolts = 931.4
const electronMassInKilograms = 9.109e-31
const electronMassInMegaElectronVolts = 0.5109
const electronMassInAtomicMassConstants = 0.0005485
const protonMassInKilograms = 1.672e-27
const protonMassInMegaElectronVolts = 938.2
const protonMassInAtomicMassUnits = 1.007
const neutronMassInKilograms = 1.674e-27
const neutronMassInMegaElectronVolts = 939.5
const neutronMassInInAtomicMassUnits = 1.008
const earthStandardGravityInMetersPerSecond2 = 9.80665
const hubbleConstantInKilometersPerSecondPerMegaparsec = 69.3
const hubbleConstantInPerSecond = 2.25e-18

// These are assorted conversion factors
const cm3perL = 1000
const m3perL = 10e-3
const coulombsPerEsu = 3.3356e-10
const pascalsPerAtm = 101325
const permittivityOfVacuum = 8.854e-12 // in C^2 J^-1 m^-1
const kelvinCelsiusDifference = 273.15
const boilingPointOfWaterAt1Atm = 99.97 // in degrees C (celsius)
const metersPerMile = 1609.344
const secondsPerHour = 3600
const secondsPerDay = 86400
const secondsPerYear = 31556923.487999998
const daysPerYear = 365.24217
const metersPerLightYear = 9460527659385452.904415084
const cmPerInch = 2.54
const mPerInch = 2.54e-3
const gramsPerPound = 453.59
const joulesPerKWh = 3.6e6
const dm3perGallon = 3.785
const m3perGallon = 3.785e-2
const feetPerFurlong = 600
const weeksPerFortnight = 2
const feetPerMile = 5280
const metersPerFoot = 3.28084
const daysPerWeek = 7

// all metric prefixes
const yotta = 1e24
const zetta = 1e21
const exa = 1e18
const peta = 1e15
const tera = 1e12
const giga = 1e9
const mega = 1e6
const kilo = 1e3
const hecto = 1e2
const deca = 1e1
const deci = 1e-1
const centi = 1e-2
const milli = 1e-3
const micro = 1e-6
const nano = 1e-9
const pico = 1e-12
const femto = 1e-15
const atto = 1e-18
const zepto = 1e-21
const yocto = 1e-24

const angstrom = 1e-10

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

// A formele is a "formula element": can be an atom, a grouping of atoms, or a molecule (with coefficient in front)

const nodeRoot = 1
const nodeMolecule = 2
const nodeGroup = 3
const nodeAtom = 4

// ids are the array indices
type formele struct {
	nodeType     int
	atomicNumber int
	coefficient  int
	parentNode   int
}

type formula struct {
	components []formele
}

type equation struct {
	leftSide  formula
	rightSide formula
	oneWay    bool
}

type iChemHelper interface {
	init()
	nameForSym(sbl string) string
	atomicNumberBySym(sbl string) int
	atomicMassBySym(sbl string) float64
	radioactiveBySym(sbl string) bool
	isotopeIsotopicMassBySym(sbl string, isotopeNumber int) float64
	isotopeAbundanceBySym(sbl string, isotopeNumber int) float64
	parseFormula(fmla string) formula
	parseFormulaStrict(fmla string) formula
	debugPrintFormula(fmla *formula)
	debugPrintlnFormula(fmla *formula)
	formulaToString(fmla *formula) string
	formulaToStringStrict(fmla *formula) string
	printFormula(fmla *formula)
	printFormulaStrict(fmla *formula)
	printlnFormula(fmla *formula)
	printlnFormulaStrict(fmla *formula)
	countFormulaAtoms(fmla formula) map[int]int
	compareFormulaAtoms(fmla1 map[int]int, fmla2 map[int]int) bool
	determineMolecularMass(atomCounts map[int]int) float64
	determineMolecularMassForFormula(fmla string) float64
	kelvinToCelsius(temp float64) float64
	celsiusToKelvin(temp float64) float64
	kelvinToFahrenheit(temp float64) float64
	fahrenheitToKelvin(temp float64) float64
	celsiusToFahrenheit(temp float64) float64
	fahrenheitToCelsius(temp float64) float64
	makeEquation(fmla1 formula, fmla2 formula) equation
	parseEquation(eqn string) equation
	printEquation(eqn *equation)
	printEquationStrict(eqn *equation)
	printlnEquation(eqn *equation)
	printlnEquationStrict(eqn *equation)
	balanceEquation(eqn *equation)
}

type chemHelper struct {
	elem         []element
	atomicNumber []int
	symbol       map[string]int
	isos         []isotope
	isolookup    map[int]map[int]int
	prefixes     map[string]float64
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
	relativeIsotopeMassesTable := []isotope{
		{5, 10, 10.013, 0.1961},
		{5, 11, 999.999, 0.8039},
		{6, 12, 12.0, 0.98892},
		{6, 13, 13.003354, .01108},
		{10, 20, 19.99212, 0.9000},
		{10, 21, 20.99316, 0.0027},
		{10, 22, 21.99132, 0.0973},
		{14, 28, 27.97693, 0.9221},
		{14, 29, 28.97649, 0.0470},
		{14, 30, 29.97376, 0.0309},
		{40, 90, 999.999, 99.99},
		{40, 91, 90.9056, 0.1127},
		{40, 92, 91.9050, 0.1717},
		{40, 94, 93.9063, 0.1733},
		{40, 96, 95.9083, 0.0278},
		{92, 238, 238.05078826, 0.992745}}
	return relativeIsotopeMassesTable
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

func isRuneUppercaseLetter(ch rune) bool {
	// This works because chemical symbols never use characters outside the ASCII range
	return (ch >= 65) && (ch <= 90)
}

func isRuneLowercaseLetter(ch rune) bool {
	// This works because chemical symbols never use characters outside the ASCII range
	return (ch >= 97) && (ch <= 122)
}

func isRuneOpenParen(ch rune) bool {
	return ch == 40
}

func isRuneCloseParen(ch rune) bool {
	return ch == 41
}

// generic error-checker -- substitute with particular error checker when
// appropriate
// func checkError(err error, message string) {
func checkError(err error) {
	if err == nil {
		return
	}
	panic(err.Error())
}

func strToInt(str string) int {
	result, err := strconv.ParseInt(str, 10, 32)
	checkError(err)
	return int(result)
}

func intToStr(ival int) string {
	return strconv.FormatInt(int64(ival), 10)
}

// returns new stack pointer and point to new stack itself if it had to be expanded
func pushNode(currentNode int, nodeStack []int, stackSP int) ([]int, int) {
	if stackSP == len(nodeStack) {
		nodeStack = append(nodeStack, currentNode)
		stackSP = len(nodeStack)
	} else {
		nodeStack[stackSP] = currentNode
		stackSP++
	}
	return nodeStack, stackSP
}

// returns ID of new node
func addNode(result *formula, newNode formele) int {
	result.components = append(result.components, newNode)
	return len(result.components) - 1
}

// returns id of new node and new stack pointer
func addAndPushNode(result *formula, newNode formele, currentNode int, nodeStack []int, stackSP int) (int, []int, int) {
	newId := addNode(result, newNode)
	nodeStack, stackSP = pushNode(currentNode, nodeStack, stackSP)
	return newId, nodeStack, stackSP
}

// returns ID of node and new stack pointer
func popNode(nodeStack []int, stackSP int) (int, int) {
	stackSP--
	return nodeStack[stackSP], stackSP
}

func handlePostParenCount(strict bool, position int, result *formula, postParenCountNode int, postParenCountStr string) (bool, bool) {
	if postParenCountStr == "" {
		if strict {
			panic("Cannot parse formula: Cannot close parenthesis in strict mode without specifying group count: position: " + intToStr(position))
		} else {
			postParenCountStr = "1"
		}
	}
	postParenCountNum := strToInt(postParenCountStr)
	result.components[postParenCountNode].coefficient = postParenCountNum
	inPostParenCount := false
	inNothing := true
	return inPostParenCount, inNothing
}

func addAtom(self *chemHelper, strict bool, position int, result *formula, currentNode int, symbolStr string, atomCountStr string) {
	var atomCountNum int
	if atomCountStr == "" {
		if strict {
			panic("Cannot parse formula: atom count cannot be missing in strict mode: position: " + intToStr(position))
		} else {
			atomCountNum = 1
		}
	} else {
		atomCountNum = strToInt(atomCountStr)
	}
	newNode := formele{nodeAtom, self.elem[self.symbol[symbolStr]].atomicNumber, atomCountNum, currentNode}
	_ = addNode(result, newNode)
}

func referenceaddAtom(self *chemHelper, strict bool, position int, result *formula, currentNode int, symbolStr string, atomCountStr string) {
	var atomCountNum int
	if atomCountStr == "" {
		if strict {
			panic("Cannot parse formula: atom count cannot be missing in strict mode: position: " + intToStr(position))
		} else {
			atomCountNum = 1
		}
	} else {
		atomCountNum = strToInt(atomCountStr)
	}
	newNode := formele{nodeAtom, self.elem[self.symbol[symbolStr]].atomicNumber, atomCountNum, currentNode}
	_ = addNode(result, newNode)
}

func moleculePushNode(coefNum int, result *formula, currentNode int, nodeStack []int, stackSP int) (int, []int, int) {
	newNode := formele{nodeMolecule, 0, coefNum, currentNode}
	currentNode, nodeStack, stackSP = addAndPushNode(result, newNode, currentNode, nodeStack, stackSP)
	return currentNode, nodeStack, stackSP
}

func pushImplicitMoleculeGroupParen(coefNum int, result *formula, currentNode int, nodeStack []int, stackSP int) (int, []int, int) {
	// molecule push
	if currentNode == 0 {
		currentNode, nodeStack, stackSP = moleculePushNode(coefNum, result, currentNode, nodeStack, stackSP)
	}
	newNode := formele{nodeGroup, 0, 0, currentNode} // We'll fill in coef at the end
	currentNode, nodeStack, stackSP = addAndPushNode(result, newNode, currentNode, nodeStack, stackSP)
	return currentNode, nodeStack, stackSP
}

func needToPushImplicitGroupParen(result *formula, currentNode int) bool {
	if currentNode == 0 {
		return true
	}
	if result.components[currentNode].nodeType == nodeMolecule {
		return true
	}
	return false
}

func handleCloseParen(self *chemHelper, strict bool, position int, result *formula, nodeStack []int, stackSP int, currentNode int, inCoef bool, inSymbol bool, symbolStr string, inAtomCount bool, atomCountStr string, inPostParenCount bool, postParenCountStr string, postParenCountNode int) (bool, bool, int, int, int, bool, string, bool) {
	var inNothing bool
	if inCoef {
		panic("Cannot parse formula: Close parenthesis encountered in a molecule coefficient: position:" + intToStr(position))
	}
	if inAtomCount {
		if (!strict) && (currentNode == 0) {
			// shouldn't this be impossible??
			panic("Cannot parse formula: Close parenthesis encountered at top level: position:" + intToStr(position))
		}
		addAtom(self, strict, position, result, currentNode, symbolStr, atomCountStr)
		inAtomCount = false
	}
	if !strict {
		if inSymbol {
			if currentNode == 0 {
				// shouldn't this be impossible??
				panic("Cannot parse formula: Close parenthesis encountered at top level: position:" + intToStr(position))
			}
			addAtom(self, strict, position, result, currentNode, symbolStr, "1")
			inSymbol = false
		}
	}
	if inPostParenCount {
		// we have to check for this because multiple nested parenthesis are possible
		inPostParenCount, inNothing = handlePostParenCount(strict, position, result, postParenCountNode, postParenCountStr)
	}
	postParenCountNode = currentNode
	currentNode, stackSP = popNode(nodeStack, stackSP)
	inPostParenCount = true
	postParenCountStr = ""
	inNothing = false // after we close the parenthesis, we are not in a coefficient, we are not in a symbol, we are not in an atom count, and we are not in a post-parenthesis count (any more), so... we are in nothing
	return inSymbol, inAtomCount, postParenCountNode, currentNode, stackSP, inPostParenCount, postParenCountStr, inNothing
}

func referencehandleCloseParen(self *chemHelper, strict bool, position int, result *formula, nodeStack []int, stackSP int, currentNode int, inSymbol bool, symbolStr string, inAtomCount bool, atomCountStr string, inPostParenCount bool, postParenCountStr string, postParenCountNode int) (bool, bool, int, int, int, bool, string) {
	if inAtomCount {
		if (!strict) && (currentNode == 0) {
			// shouldn't this be impossible??
			panic("Close parenthesis encountered at top level: position:" + intToStr(position))
		}
		addAtom(self, strict, position, result, currentNode, symbolStr, atomCountStr)
		inAtomCount = false
	}
	if !strict {
		if inSymbol {
			if currentNode == 0 {
				// shouldn't this be impossible??
				panic("Close parenthesis encountered at top level: position:" + intToStr(position))
			}
			addAtom(self, strict, position, result, currentNode, symbolStr, "1")
			inSymbol = false
		}
	}
	if inPostParenCount {
		// we have to check for this because multiple nested parenthesis are possible
		// underscore is hack
		_, inPostParenCount = handlePostParenCount(strict, position, result, postParenCountNode, postParenCountStr)
	}
	postParenCountNode = currentNode
	currentNode, stackSP = popNode(nodeStack, stackSP)
	inPostParenCount = true
	postParenCountStr = ""
	return inAtomCount, inSymbol, postParenCountNode, currentNode, stackSP, inPostParenCount, postParenCountStr
}

func assertInOneOf(inCoef bool, inSymbol bool, inAtomCount bool, inPostParenCount bool, inNothing bool) {
	count := 0
	if inCoef {
		count++
	}
	if inSymbol {
		count++
	}
	if inAtomCount {
		count++
	}
	if inPostParenCount {
		count++
	}
	if inNothing {
		count++
	}
	if count != 1 {
		panic("Cannot parse formula: Assert(in one of) failed")
	}
}

func parseFormulaInternal(self *chemHelper, strict bool, fmla string) formula {
	// strict mode requires all coefficients and parenthesis to be explicitly specified
	// formulas look like: 2Hg + O2
	// formulas can have subgroups:
	var result formula
	var nodeStack []int
	nodeStack = make([]int, 0)
	stackSP := 0
	rootNode := formele{nodeRoot, 0, 0, -1}
	result.components = append(result.components, rootNode)
	currentNode := 0 // 0 == root node
	var coefNum int
	coefStr := ""
	inCoef := true
	symbolStr := ""
	inSymbol := false
	atomCountStr := ""
	inAtomCount := false
	postParenCountStr := ""
	inPostParenCount := false
	var postParenCountNode int
	position := 0
	inNothing := false
	inImplicitMoleculeGroup := false
	for _, char := range fmla {
		assertInOneOf(inCoef, inSymbol, inAtomCount, inPostParenCount, inNothing)
		position++
		valid := false
		if isRuneNumeric(char) {
			valid = true
			if inCoef {
				coefStr += string(char)
			} else {
				if inSymbol {
					inSymbol = false
					inAtomCount = true
					atomCountStr = string(char)
				} else {
					if inAtomCount {
						atomCountStr += string(char)
					} else {
						if inPostParenCount {
							postParenCountStr += string(char)
						} else {
						}
					}
				}
			}
		} else {
			if inPostParenCount {
				inPostParenCount, inNothing = handlePostParenCount(strict, position, &result, postParenCountNode, postParenCountStr)
			}
		}
		if isRuneSpace(char) || isRunePlus(char) {
			if !strict {
				if inImplicitMoleculeGroup {
					inSymbol, inAtomCount, postParenCountNode, currentNode, stackSP, inPostParenCount, postParenCountStr, inNothing = handleCloseParen(self, strict, position, &result, nodeStack, stackSP, currentNode, inCoef, inSymbol, symbolStr, inAtomCount, atomCountStr, inPostParenCount, postParenCountStr, postParenCountNode)
					result.components[postParenCountNode].coefficient = 1
					inPostParenCount = false
					inImplicitMoleculeGroup = false
				}
			}
			if isRuneSpace(char) {
				valid = true
				if !inCoef {
					inNothing = true
				}
			}
			if isRunePlus(char) {
				valid = true
				currentNode, stackSP = popNode(nodeStack, stackSP)
				if currentNode != 0 {
					panic("Cannot parse formula: Plus character detected inside a molecule. Plus can only be used between molecules. Position: " + intToStr(position))
				}
				inCoef = true
				coefStr = ""
				inNothing = false
			}
		}
		if isRuneUppercaseLetter(char) {
			valid = true
			if inCoef || inNothing {
				if inCoef {
					if strict {
						panic("Cannot parse formula: Cannot start symbol without parenthesis in strict mode: position: " + intToStr(position))
					} else {
						if inCoef {
							if coefStr == "" {
								if strict {
									panic("Cannot parse formula: Cannot parse formula: coefficient missing: position: " + intToStr(position))
								} else {
									coefNum = 1
								}
							} else {
								coefNum = strToInt(coefStr)
							}
							inCoef = false
						}
					}
				}
				if inNothing {
					coefNum = 1
				}
				if needToPushImplicitGroupParen(&result, currentNode) {
					currentNode, nodeStack, stackSP = pushImplicitMoleculeGroupParen(coefNum, &result, currentNode, nodeStack, stackSP)
					inImplicitMoleculeGroup = true
				}
			}
			if inSymbol {
				if strict {
					panic("Cannot parse formula: Cannot switch symbols within a symbol in strict mode: position: " + intToStr(position))
				} else {
					addAtom(self, strict, position, &result, currentNode, symbolStr, "1")
				}
			}
			if inAtomCount {
				if !strict {
					if needToPushImplicitGroupParen(&result, currentNode) {
						currentNode, nodeStack, stackSP = pushImplicitMoleculeGroupParen(coefNum, &result, currentNode, nodeStack, stackSP)
						inImplicitMoleculeGroup = true
					}
				}
				addAtom(self, strict, position, &result, currentNode, symbolStr, atomCountStr)
				inAtomCount = false
			}
			symbolStr = string(char)
			inSymbol = true
			inNothing = false
		}
		if isRuneLowercaseLetter(char) {
			valid = true
			if inCoef {
				if strict {
					panic("Cannot parse formula: Cannot start symbol without parenthesis: position: " + intToStr(position))
				}
			}
			if !inSymbol {
				panic("Cannot parse formula: Cannot start symbol with a lowercase letter: position: " + intToStr(position))
			}
			symbolStr += string(char)
		}
		if isRuneOpenParen(char) {
			valid = true
			// molecule push
			if inCoef {
				if coefStr == "" {
					if strict {
						panic("Cannot parse formula: coefficient missing: position: " + intToStr(position))
					} else {
						coefNum = 1
					}
				} else {
					coefNum = strToInt(coefStr)
				}
				inCoef = false
			}
			if currentNode == 0 {
				currentNode, nodeStack, stackSP = moleculePushNode(coefNum, &result, currentNode, nodeStack, stackSP)
				inImplicitMoleculeGroup = false
			}
			if inAtomCount {
				addAtom(self, strict, position, &result, currentNode, symbolStr, atomCountStr)
				inAtomCount = false
			}
			if inSymbol {
				if strict {
					panic("Cannot parse formula: Atom counts cannot be omitted in strict mode: position: " + intToStr(position))
				} else {
					addAtom(self, strict, position, &result, currentNode, symbolStr, "1")
					inSymbol = false
				}
			}
			if !strict {
				if inImplicitMoleculeGroup {
					inPostParenCount, inNothing = handlePostParenCount(strict, position, &result, currentNode, "1")
					currentNode, stackSP = popNode(nodeStack, stackSP)
					inImplicitMoleculeGroup = false
				}
			}
			newNode := formele{nodeGroup, 0, 0, currentNode} // We'll fill in coef at the end
			currentNode, nodeStack, stackSP = addAndPushNode(&result, newNode, currentNode, nodeStack, stackSP)
			inNothing = true
		}
		if isRuneCloseParen(char) {
			valid = true
			inSymbol, inAtomCount, postParenCountNode, currentNode, stackSP, inPostParenCount, postParenCountStr, inNothing = handleCloseParen(self, strict, position, &result, nodeStack, stackSP, currentNode, inCoef, inSymbol, symbolStr, inAtomCount, atomCountStr, inPostParenCount, postParenCountStr, postParenCountNode)
		}
		if !valid {
			panic("Cannot parse formula: invalid character: position: " + intToStr(position))
		}
	}
	if !strict {
		if inImplicitMoleculeGroup {
			inSymbol, inAtomCount, postParenCountNode, currentNode, stackSP, inPostParenCount, postParenCountStr, inNothing = handleCloseParen(self, strict, position, &result, nodeStack, stackSP, currentNode, inCoef, inSymbol, symbolStr, inAtomCount, atomCountStr, inPostParenCount, postParenCountStr, postParenCountNode)
			result.components[postParenCountNode].coefficient = 1
			inPostParenCount = false
			inImplicitMoleculeGroup = false
		}
	}
	if inPostParenCount {
		inPostParenCount, inNothing = handlePostParenCount(strict, position, &result, postParenCountNode, postParenCountStr)
	}
	return result
}

func referenceparseFormulaInternal(self *chemHelper, strict bool, fmla string) formula {
	// strict mode requires all coefficients and parenthesis to be explicitly specified
	// formulas look like: 2Hg + O2
	// formulas can have subgroups: Ca(NO3)2
	var result formula
	var nodeStack []int
	nodeStack = make([]int, 0)
	stackSP := 0
	rootNode := formele{nodeRoot, 0, 0, -1}
	result.components = append(result.components, rootNode)
	currentNode := 0 // 0 == root node
	var coefNum int
	coefStr := ""
	inCoef := true
	symbolStr := ""
	inSymbol := false
	atomCountStr := ""
	inAtomCount := false
	// var postParenCountNum int
	postParenCountStr := ""
	inPostParenCount := false
	var postParenCountNode int
	position := 0
	inImplicitMoleculeGroup := false
	for _, char := range fmla {
		position++
		valid := false
		if isRuneNumeric(char) {
			valid = true
			if inCoef {
				coefStr += string(char)
			} else {
				if inSymbol {
					inSymbol = false
					inAtomCount = true
					atomCountStr = string(char)
				} else {
					if inAtomCount {
						atomCountStr += string(char)
					} else {
						if inPostParenCount {
							postParenCountStr += string(char)
						} else {
							panic("dont know what to do with this number: position: " + intToStr(position))
						}
					}
				}
			}
		} else {
			if inPostParenCount {
				inCoef, inPostParenCount = handlePostParenCount(strict, position, &result, postParenCountNode, postParenCountStr)
			}
		}
		if isRuneSpace(char) || isRunePlus(char) {
			if !strict {
				if inImplicitMoleculeGroup {
					// hack
					// inAtomCount, inSymbol, postParenCountNode, currentNode, stackSP, inPostParenCount, postParenCountStr = handleCloseParen(self, strict, position, &result, nodeStack, stackSP, currentNode, inSymbol, symbolStr, inAtomCount, atomCountStr, inPostParenCount, postParenCountStr, postParenCountNode)
					inSymbol, inAtomCount, postParenCountNode, currentNode, stackSP, inPostParenCount, postParenCountStr, _ = handleCloseParen(self, strict, position, &result, nodeStack, stackSP, currentNode, inCoef, inSymbol, symbolStr, inAtomCount, atomCountStr, inPostParenCount, postParenCountStr, postParenCountNode)
					result.components[postParenCountNode].coefficient = 1
					inPostParenCount = false
					inImplicitMoleculeGroup = false
				}
			}
			if isRuneSpace(char) {
				valid = true
			}
			if isRunePlus(char) {
				valid = true
				currentNode, stackSP = popNode(nodeStack, stackSP)
				if currentNode != 0 {
					panic("Plus character detected inside a molecule. Plus can only be used between molecules. Position: " + intToStr(position))
				}
				inCoef = true
				coefStr = ""
			}
		}
		if isRuneUppercaseLetter(char) {
			valid = true
			if inCoef {
				if strict {
					panic("Cannot start symbol without parenthesis in strict mode: position: " + intToStr(position))
				} else {
					if coefStr == "" {
						if strict {
							panic("Cannot parse formula: coefficient missing: position: " + intToStr(position))
						} else {
							coefNum = 1
						}
					} else {
						coefNum = strToInt(coefStr)
					}
					inCoef = false
					if needToPushImplicitGroupParen(&result, currentNode) {
						currentNode, nodeStack, stackSP = pushImplicitMoleculeGroupParen(coefNum, &result, currentNode, nodeStack, stackSP)
						inImplicitMoleculeGroup = true
					}
				}
			}
			if inSymbol {
				if strict {
					panic("Cannot switch symbols within a symbol in strict mode: position: " + intToStr(position))
				} else {
					addAtom(self, strict, position, &result, currentNode, symbolStr, "1")
				}
			}
			if inAtomCount {
				if !strict {
					if needToPushImplicitGroupParen(&result, currentNode) {
						currentNode, nodeStack, stackSP = pushImplicitMoleculeGroupParen(coefNum, &result, currentNode, nodeStack, stackSP)
						inImplicitMoleculeGroup = true
					}
				}
				addAtom(self, strict, position, &result, currentNode, symbolStr, atomCountStr)
				inAtomCount = false
			}
			symbolStr = string(char)
			inSymbol = true
		}
		if isRuneLowercaseLetter(char) {
			valid = true
			if inCoef {
				if strict {
					panic("Cannot start symbol without parenthesis: position: " + intToStr(position))
				}
			}
			if !inSymbol {
				panic("Cannot start symbol with a lowercase letter: position: " + intToStr(position))
			}
			symbolStr += string(char)
		}
		if isRuneOpenParen(char) {
			valid = true
			// molecule push
			if inCoef {
				if coefStr == "" {
					if strict {
						panic("Cannot parse formula: coefficient missing: position: " + intToStr(position))
					} else {
						coefNum = 1
					}
				} else {
					coefNum = strToInt(coefStr)
				}
				inCoef = false
			}
			if currentNode == 0 {
				currentNode, nodeStack, stackSP = moleculePushNode(coefNum, &result, currentNode, nodeStack, stackSP)
				inImplicitMoleculeGroup = false
			}
			if inAtomCount {
				addAtom(self, strict, position, &result, currentNode, symbolStr, atomCountStr)
				inAtomCount = false
			}
			if inSymbol {
				if strict {
					panic("Atom counts cannot be omitted in strict mode: position: " + intToStr(position))
				} else {
					addAtom(self, strict, position, &result, currentNode, symbolStr, "1")
					inSymbol = false
				}
			}
			if !strict {
				if inImplicitMoleculeGroup {
					// if we're in an implicit ground and going into an explicit group, we have to pop out of the implicit group
					inCoef, inPostParenCount = handlePostParenCount(strict, position, &result, currentNode, "1")
					currentNode, stackSP = popNode(nodeStack, stackSP)
					inImplicitMoleculeGroup = false
				}
			}
			newNode := formele{nodeGroup, 0, 0, currentNode} // We'll fill in coef at the end
			currentNode, nodeStack, stackSP = addAndPushNode(&result, newNode, currentNode, nodeStack, stackSP)
		}
		if isRuneCloseParen(char) {
			valid = true
			// hack
			// inAtomCount, inSymbol, postParenCountNode, currentNode, stackSP, inPostParenCount, postParenCountStr = handleCloseParen(self, strict, position, &result, nodeStack, stackSP, currentNode, inSymbol, symbolStr, inAtomCount, atomCountStr, inPostParenCount, postParenCountStr, postParenCountNode)
			inSymbol, inAtomCount, postParenCountNode, currentNode, stackSP, inPostParenCount, postParenCountStr, _ = handleCloseParen(self, strict, position, &result, nodeStack, stackSP, currentNode, inCoef, inSymbol, symbolStr, inAtomCount, atomCountStr, inPostParenCount, postParenCountStr, postParenCountNode)
		}
		if !valid {
			panic("Invalid character: position: " + intToStr(position))
		}
	}
	if !strict {
		if inImplicitMoleculeGroup {
			// hack
			// inAtomCount, inSymbol, postParenCountNode, currentNode, stackSP, inPostParenCount, postParenCountStr = handleCloseParen(self, strict, position, &result, nodeStack, stackSP, currentNode, inSymbol, symbolStr, inAtomCount, atomCountStr, inPostParenCount, postParenCountStr, postParenCountNode)
			inSymbol, inAtomCount, postParenCountNode, currentNode, stackSP, inPostParenCount, postParenCountStr, _ = handleCloseParen(self, strict, position, &result, nodeStack, stackSP, currentNode, inCoef, inSymbol, symbolStr, inAtomCount, atomCountStr, inPostParenCount, postParenCountStr, postParenCountNode)
			result.components[postParenCountNode].coefficient = 1
			inPostParenCount = false
			inImplicitMoleculeGroup = false
		}
	}
	if inPostParenCount {
		inCoef, inPostParenCount = handlePostParenCount(strict, position, &result, postParenCountNode, postParenCountStr)
	}
	return result
}

func determineIfGroupNodeDispayableInNormalMode(fmla *formula, nodeNum int) bool {
	if fmla.components[nodeNum].coefficient != 1 {
		return true
	}
	parent := fmla.components[nodeNum].parentNode
	if parent < 0 {
		// should never happen?
		return true
	}
	if fmla.components[parent].nodeType == nodeMolecule {
		return false
	}
	return true
}

func strizeOneNodeBefore(self *chemHelper, strict bool, fmla *formula, nodeNum int) string {
	switch fmla.components[nodeNum].nodeType {
	case nodeMolecule:
		if strict {
			return " + " + intToStr(fmla.components[nodeNum].coefficient)
		} else {
			if fmla.components[nodeNum].coefficient == 1 {
				return " + "
			} else {
				return " + " + intToStr(fmla.components[nodeNum].coefficient)
			}
		}
	case nodeGroup:
		if strict {
			return "("
		} else {
			if determineIfGroupNodeDispayableInNormalMode(fmla, nodeNum) {
				return "("
			}
		}
	case nodeAtom:
		if strict {
			return self.elem[self.atomicNumber[fmla.components[nodeNum].atomicNumber]].symbol + intToStr(fmla.components[nodeNum].coefficient)
		} else {
			if fmla.components[nodeNum].coefficient == 1 {
				return self.elem[self.atomicNumber[fmla.components[nodeNum].atomicNumber]].symbol
			} else {
				return self.elem[self.atomicNumber[fmla.components[nodeNum].atomicNumber]].symbol + intToStr(fmla.components[nodeNum].coefficient)
			}
		}
	}
	return ""
}

func referencestrizeOneNodeBefore(self *chemHelper, strict bool, fmla *formula, nodeNum int) string {
	switch fmla.components[nodeNum].nodeType {
	case nodeMolecule:
		if strict {
			return " + " + intToStr(fmla.components[nodeNum].coefficient)
		} else {
			if fmla.components[nodeNum].coefficient == 1 {
				return " + "
			} else {
				return " + " + intToStr(fmla.components[nodeNum].coefficient)
			}
		}
	case nodeGroup:
		if strict {
			return "("
		} else {
			if determineIfGroupNodeDispayableInNormalMode(fmla, nodeNum) {
				return "("
			}
		}
	case nodeAtom:
		if strict {
			return self.elem[self.atomicNumber[fmla.components[nodeNum].atomicNumber]].symbol + intToStr(fmla.components[nodeNum].coefficient)
		} else {
			if fmla.components[nodeNum].coefficient == 1 {
				return self.elem[self.atomicNumber[fmla.components[nodeNum].atomicNumber]].symbol
			} else {
				return self.elem[self.atomicNumber[fmla.components[nodeNum].atomicNumber]].symbol + intToStr(fmla.components[nodeNum].coefficient)
			}
		}
	}
	return ""
}

func strizeOneNodeAfter(strict bool, fmla *formula, nodeNum int) string {
	switch fmla.components[nodeNum].nodeType {
	case nodeGroup:
		if strict {
			return ")" + intToStr(fmla.components[nodeNum].coefficient)
		} else {
			if determineIfGroupNodeDispayableInNormalMode(fmla, nodeNum) {
				if fmla.components[nodeNum].coefficient == 1 {
					return ")"
				} else {
					return ")" + intToStr(fmla.components[nodeNum].coefficient)
				}
			}
		}
	}
	return ""
}

func recurseStrizeFormulaNode(self *chemHelper, strict bool, fmla *formula, nodeNum int) string {
	result := strizeOneNodeBefore(self, strict, fmla, nodeNum)
	lcomp := len(fmla.components)
	for ii := 0; ii < lcomp; ii++ {
		if fmla.components[ii].parentNode == nodeNum {
			result += recurseStrizeFormulaNode(self, strict, fmla, ii)
		}
	}
	result += strizeOneNodeAfter(strict, fmla, nodeNum)
	return result
}

func referencerecurseStrizeFormulaNode(self *chemHelper, strict bool, fmla *formula, nodeNum int) string {
	result := strizeOneNodeBefore(self, strict, fmla, nodeNum)
	lcomp := len(fmla.components)
	for ii := 0; ii < lcomp; ii++ {
		if fmla.components[ii].parentNode == nodeNum {
			result += recurseStrizeFormulaNode(self, strict, fmla, ii)
		}
	}
	result += strizeOneNodeAfter(strict, fmla, nodeNum)
	return result
}

func printIndent(indentlvl int) {
	for ii := 0; ii < (indentlvl * 4); ii++ {
		fmt.Print(" ")
	}
}

func debugPrintOneNode(self *chemHelper, fmla *formula, nodeNum int, indentlvl int) {
	printIndent(indentlvl)
	// uncomment the next line to get node ID numbers for debugging
	// fmt.Print("Node " + intToStr(nodeNum) + ": ")
	switch fmla.components[nodeNum].nodeType {
	case nodeRoot:
		fmt.Println("Root node")
	case nodeMolecule:
		fmt.Println("Molecule")
		printIndent(indentlvl + 1)
		fmt.Println("Number of molecules:", fmla.components[nodeNum].coefficient)
	case nodeGroup:
		fmt.Println("Group")
		printIndent(indentlvl + 1)
		fmt.Println("Number of group repetitions:", fmla.components[nodeNum].coefficient)
	case nodeAtom:
		fmt.Println("Atom")
		printIndent(indentlvl + 1)
		fmt.Println("Atom element number:", fmla.components[nodeNum].atomicNumber)
		printIndent(indentlvl + 1)
		fmt.Println("Atom symbol:", self.elem[self.atomicNumber[fmla.components[nodeNum].atomicNumber]].symbol)
		printIndent(indentlvl + 1)
		fmt.Println("Atom name:", self.elem[self.atomicNumber[fmla.components[nodeNum].atomicNumber]].name)
		printIndent(indentlvl + 1)
		fmt.Println("Number of atoms:", fmla.components[nodeNum].coefficient)
	}
}

func referencedebugPrintOneNode(self *chemHelper, fmla *formula, nodeNum int, indentlvl int) {
	printIndent(indentlvl)
	fmt.Print("Node " + intToStr(nodeNum) + ": ")
	switch fmla.components[nodeNum].nodeType {
	case nodeRoot:
		fmt.Println("Root node")
	case nodeMolecule:
		fmt.Println("Molecule")
		printIndent(indentlvl + 1)
		fmt.Println("Number of molecules:", fmla.components[nodeNum].coefficient)
	case nodeGroup:
		fmt.Println("Group")
		printIndent(indentlvl + 1)
		fmt.Println("Number of group repetitions:", fmla.components[nodeNum].coefficient)
	case nodeAtom:
		fmt.Println("Atom")
		printIndent(indentlvl + 1)
		fmt.Println("Atom element number:", fmla.components[nodeNum].atomicNumber)
		printIndent(indentlvl + 1)
		fmt.Println("Atom symbol:", self.elem[self.atomicNumber[fmla.components[nodeNum].atomicNumber]].symbol)
		printIndent(indentlvl + 1)
		fmt.Println("Atom name:", self.elem[self.atomicNumber[fmla.components[nodeNum].atomicNumber]].name)
		printIndent(indentlvl + 1)
		fmt.Println("Number of atoms:", fmla.components[nodeNum].coefficient)
	}
}

func recursePrettyPrintFormulaNode(self *chemHelper, fmla *formula, nodeNum int, indentlvl int) {
	printIndent(indentlvl)
	fmt.Println("[")
	debugPrintOneNode(self, fmla, nodeNum, indentlvl)
	lcomp := len(fmla.components)
	for ii := 0; ii < lcomp; ii++ {
		if fmla.components[ii].parentNode == nodeNum {
			recursePrettyPrintFormulaNode(self, fmla, ii, indentlvl+1)
		}
	}
	printIndent(indentlvl)
	fmt.Println("]")
}

func referencerecursePrettyPrintFormulaNode(self *chemHelper, fmla *formula, nodeNum int, indentlvl int) {
	printIndent(indentlvl)
	fmt.Println("[")
	debugPrintOneNode(self, fmla, nodeNum, indentlvl)
	lcomp := len(fmla.components)
	for ii := 0; ii < lcomp; ii++ {
		if fmla.components[ii].parentNode == nodeNum {
			recursePrettyPrintFormulaNode(self, fmla, ii, indentlvl+1)
		}
	}
	printIndent(indentlvl)
	fmt.Println("]")
}

func testParsingFormulaCheckMatch(fmA *formula, fmB *formula) {
	if len(fmA.components) != len(fmB.components) {
		panic("Parse Formula Test: Check Match: formulas have unequal lengths: " + intToStr(len(fmA.components)) + " vs " + intToStr(len(fmB.components)))
	}
	for ii := 0; ii < len(fmA.components); ii++ {
		if fmA.components[ii].nodeType != fmB.components[ii].nodeType {
			panic("Parse Formula Test: Node types don't match. Entry: " + intToStr(ii))
		}
		if fmA.components[ii].atomicNumber != fmB.components[ii].atomicNumber {
			panic("Parse Formula Test: Atomic numbers don't match. Entry: " + intToStr(ii))
		}
		if fmA.components[ii].coefficient != fmB.components[ii].coefficient {
			panic("Parse Formula Test: Coefficients don't match. Entry: " + intToStr(ii))
		}
		if fmA.components[ii].parentNode != fmB.components[ii].parentNode {
			panic("Parse Formula Test: Parent nodes don't match. Entry: " + intToStr(ii))
		}
	}
	fmt.Println("Check match passed.")
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
	idx, ok := self.symbol[sbl]
	if !ok {
		panic("nameForSym: symbol does not exist: " + sbl)
	}
	return self.elem[idx].name
}

func (self *chemHelper) atomicNumberBySym(sbl string) int {
	idx, ok := self.symbol[sbl]
	if !ok {
		panic("atomicNumberBySym: symbol does not exist: " + sbl)
	}
	return self.elem[idx].atomicNumber
}
func (self *chemHelper) atomicMassBySym(sbl string) float64 {
	idx, ok := self.symbol[sbl]
	if !ok {
		panic("atomicMassBySym: symbol does not exist: " + sbl)
	}
	return self.elem[idx].relativeAtomicMass
}
func (self *chemHelper) radioactiveBySym(sbl string) bool {
	idx, ok := self.symbol[sbl]
	if !ok {
		panic("radioactiveBySym: symbol does not exist: " + sbl)
	}
	return self.elem[idx].radioactive
}

func (self *chemHelper) isotopeIsotopicMassBySym(sbl string, isotopeNumber int) float64 {
	// this one-liner would work if we didn't have to check for missing entries in tables!
	// return self.isos[self.isolookup[self.elem[self.symbol[sbl]].atomicNumber][isotopeNumber]].isotopicMass
	elemIdx, ok := self.symbol[sbl]
	if !ok {
		panic("isotopeIsotopicMassBySym: symbol does not exist: " + sbl)
	}
	atomnumber := self.elem[elemIdx].atomicNumber
	isoForElem, ok := self.isolookup[atomnumber]
	if !ok {
		panic("isotopeIsotopicMassBySym: isotope table does not have isotopes for element: " + sbl)
	}
	isoIdx, ok := isoForElem[isotopeNumber]
	if !ok {
		panic("isotopeIsotopicMassBySym: isotope table has isotopes for element but not the isotope requested: element: " + sbl + ", isotope: " + intToStr(isotopeNumber))
	}
	return self.isos[isoIdx].isotopicMass
}

func (self *chemHelper) isotopeAbundanceBySym(sbl string, isotopeNumber int) float64 {
	// this one-liner would work if we didn't have to check for missing entries in tables!
	// return self.isos[self.isolookup[self.elem[self.symbol[sbl]].atomicNumber][isotopeNumber]].abundance
	elemIdx, ok := self.symbol[sbl]
	if !ok {
		panic("isotopeIsotopicMassBySym: symbol does not exist: " + sbl)
	}
	atomnumber := self.elem[elemIdx].atomicNumber
	isoForElem, ok := self.isolookup[atomnumber]
	if !ok {
		panic("isotopeIsotopicMassBySym: isotope table does not have isotopes for element: " + sbl)
	}
	isoIdx, ok := isoForElem[isotopeNumber]
	if !ok {
		panic("isotopeIsotopicMassBySym: isotope table has isotopes for element but not the isotope requested: element: " + sbl + ", isotope: " + intToStr(isotopeNumber))
	}
	return self.isos[isoIdx].abundance
}

func (self *chemHelper) parseFormula(fmla string) formula {
	return parseFormulaInternal(self, false, fmla)
}

func (self *chemHelper) parseFormulaStrict(fmla string) formula {
	return parseFormulaInternal(self, true, fmla)
}

func (self *chemHelper) debugPrintFormula(fmla *formula) {
	indentlvl := 0
	recursePrettyPrintFormulaNode(self, fmla, 0, indentlvl)
}

func (self *chemHelper) debugPrintlnFormula(fmla *formula) {
	indentlvl := 0
	recursePrettyPrintFormulaNode(self, fmla, 0, indentlvl)
	fmt.Println("")
}

func (self *chemHelper) formulaToString(fmla *formula) string {
	return recurseStrizeFormulaNode(self, false, fmla, 0)[3:]
}

func (self *chemHelper) formulaToStringStrict(fmla *formula) string {
	return recurseStrizeFormulaNode(self, true, fmla, 0)[3:]
}

func (self *chemHelper) printFormula(fmla *formula) {
	fmt.Print(self.formulaToString(fmla))
}

func (self *chemHelper) printFormulaStrict(fmla *formula) {
	fmt.Print(self.formulaToStringStrict(fmla))
}

func (self *chemHelper) printlnFormula(fmla *formula) {
	fmt.Println(self.formulaToString(fmla))
}

func (self *chemHelper) printlnFormulaStrict(fmla *formula) {
	fmt.Println(self.formulaToStringStrict(fmla))
}

func (self *chemHelper) countFormulaAtoms(fmla formula) map[int]int {
	var result map[int]int
	result = make(map[int]int)
	lform := len(fmla.components)
	for ii := 0; ii < lform; ii++ {
		if fmla.components[ii].nodeType == nodeAtom {
			atomicNumber := fmla.components[ii].atomicNumber
			numAtoms := fmla.components[ii].coefficient
			parent := fmla.components[ii].parentNode
			for parent != 0 {
				numAtoms *= fmla.components[parent].coefficient
				parent = fmla.components[parent].parentNode
			}
			_, ok := result[atomicNumber]
			if ok {
				result[atomicNumber] += numAtoms
			} else {
				result[atomicNumber] = numAtoms
			}
		}
	}
	return result
}

func (self *chemHelper) compareFormulaAtoms(fmla1 map[int]int, fmla2 map[int]int) bool {
	// this is a dumb and blunt function that could be optimized
	for atomNum, numAtoms1 := range fmla1 {
		numAtoms2, ok := fmla2[atomNum]
		if !ok {
			return false
		}
		if numAtoms2 != numAtoms1 {
			return false
		}
	}
	for atomNum, numAtoms2 := range fmla2 {
		numAtoms1, ok := fmla1[atomNum]
		if !ok {
			return false
		}
		if numAtoms1 != numAtoms2 {
			return false
		}
	}
	return true
}

func (self *chemHelper) determineMolecularMass(atomCounts map[int]int) float64 {
	var result float64
	result = 0.0
	for elemtNum, numberOfAtoms := range atomCounts {
		result += self.elem[self.atomicNumber[elemtNum]].relativeAtomicMass * float64(numberOfAtoms)
	}
	return result
}

func (self *chemHelper) determineMolecularMassForFormula(fmla string) float64 {
	return self.determineMolecularMass(self.countFormulaAtoms(self.parseFormula(fmla)))
}

func (self *chemHelper) kelvinToCelsius(temp float64) float64 {
	return temp - kelvinCelsiusDifference
}

func (self *chemHelper) celsiusToKelvin(temp float64) float64 {
	return temp + kelvinCelsiusDifference
}

func (self *chemHelper) kelvinToFahrenheit(temp float64) float64 {
	return ((temp - kelvinCelsiusDifference) * 1.8) + 32
}

func (self *chemHelper) fahrenheitToKelvin(temp float64) float64 {
	return ((temp - 32) / 1.8) + kelvinCelsiusDifference
}

func (self *chemHelper) celsiusToFahrenheit(temp float64) float64 {
	return (temp * 1.8) + 32
}

func (self *chemHelper) fahrenheitToCelsius(temp float64) float64 {
	return (temp - 32) / 1.8
}

func (self *chemHelper) makeEquation(leftSide formula, rightSide formula, oneWay bool) equation {
	var result equation
	result.leftSide = leftSide
	result.rightSide = rightSide
	result.oneWay = oneWay
	return result
}

func (self *chemHelper) parseEquation(eqn string) equation {
	var oneWay bool
	var leftSide formula
	var rightSide formula

	ii := strings.Index(eqn, "<->")
	if ii < 0 {
		ii := strings.Index(eqn, "->")
		if ii < 0 {
			panic("cannot find equation arrow")
		}
		oneWay = true
		leftSide = self.parseFormula(eqn[:ii])
		rightSide = self.parseFormula(eqn[ii+2:])
	} else {
		oneWay = false
		leftSide = self.parseFormula(eqn[:ii])
		rightSide = self.parseFormula(eqn[ii+3:])
	}
	return self.makeEquation(leftSide, rightSide, oneWay)
}

func (self *chemHelper) printEquation(eqn *equation) {
	self.printFormula(&eqn.leftSide)
	fmt.Print(" ")
	if !eqn.oneWay {
		fmt.Print("<")
	}
	fmt.Print("-> ")
	self.printFormula(&eqn.rightSide)
}

func (self *chemHelper) printEquationStrict(eqn *equation) {
	self.printFormulaStrict(&eqn.leftSide)
	fmt.Print(" ")
	if !eqn.oneWay {
		fmt.Print("<")
	}
	fmt.Print("-> ")
	self.printFormulaStrict(&eqn.rightSide)
}

func (self *chemHelper) printlnEquation(eqn *equation) {
	self.printEquation(eqn)
	fmt.Println("")
}

func (self *chemHelper) printlnEquationStrict(eqn *equation) {
	self.printEquationStrict(eqn)
	fmt.Println("")
}

func (self *chemHelper) balanceEquation(eqn *equation) {
	type coefIdx struct {
		side int
		node int
	}
	coefs := make([]coefIdx, 0)
	// find our coefficient nodes and initialize
	for node := 0; node < len(eqn.leftSide.components); node++ {
		if eqn.leftSide.components[node].nodeType == nodeMolecule {
			coefs = append(coefs, coefIdx{0, node})
			eqn.leftSide.components[node].coefficient = 1
		}
	}
	for node := 0; node < len(eqn.rightSide.components); node++ {
		if eqn.rightSide.components[node].nodeType == nodeMolecule {
			coefs = append(coefs, coefIdx{1, node})
			eqn.rightSide.components[node].coefficient = 1
		}
	}
	// Ok, now we find the coefficients that balance the equation!!!
	// This is a very inefficient brute force algorithm
	// However, since chemical equation coefficients are almost always
	// small (single-digit numbers), it's still instant from the
	// user's point of view!
	max := 2
	pos := 0
	var overMax bool
	for {
		// debug
		// fmt.Print("Current coefficients:")
		// for ii := 0; ii < len(coefs); ii++ {
		// 	fmt.Print(" ")
		// 	if coefs[ii].side == 0 {
		// 		fmt.Print(eqn.leftSide.components[coefs[ii].node].coefficient)
		// 	} else {
		// 		fmt.Print(eqn.rightSide.components[coefs[ii].node].coefficient)
		// 	}
		// }
		// fmt.Println("")
		// compare
		// these could be optimized to eliminate the repeat memory allocation
		leftAtoms := self.countFormulaAtoms(eqn.leftSide)
		rightAtoms := self.countFormulaAtoms(eqn.rightSide)
		match := self.compareFormulaAtoms(leftAtoms, rightAtoms)
		if match {
			return
		}
		// increment
		pos = 0
		keepGoing := true
		for keepGoing {
			if coefs[pos].side == 0 {
				eqn.leftSide.components[coefs[pos].node].coefficient++
				overMax = eqn.leftSide.components[coefs[pos].node].coefficient > max
			} else {
				eqn.rightSide.components[coefs[pos].node].coefficient++
				overMax = eqn.rightSide.components[coefs[pos].node].coefficient > max
			}
			keepGoing = false
			if overMax {
				keepGoing = true
				for ii := 0; ii <= pos; ii++ {
					if coefs[ii].side == 0 {
						eqn.leftSide.components[coefs[ii].node].coefficient = 1
					} else {
						eqn.rightSide.components[coefs[ii].node].coefficient = 1
					}
				}
				pos++
				if pos == len(coefs) {
					if max == 128 {
						panic("endless loop in coefficient search for balancing equation")
					}
					max <<= 1
					pos = 0
					keepGoing = false
				}
			}
		}
	}
	}

// ------------------------------------------------------------------------
// End chemHelper interface implementation
// ------------------------------------------------------------------------

func testParsingFormulaMain() {
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
	fmt.Println("O", chemh.atomicNumberBySym("O"))
	fmt.Println("C", chemh.atomicMassBySym("C"))
	fmt.Println("Si", chemh.radioactiveBySym("Si"))
	// this should cause error
	// fmt.Println("Xq", chemh.radioactiveBySym("Xq"))
	testFormula := chemh.parseFormula("2Hg + O2")
	fmt.Println("2Hg + O2 formula", testFormula)
	testFormula = chemh.parseFormula("7H2SO4 + 2Hg + O2")
	fmt.Println("7H2SO4 + 2Hg + O2 formula", testFormula)
	atomCounts := chemh.countFormulaAtoms(testFormula)
	fmt.Println("atomCounts", atomCounts)
	molarMass := chemh.determineMolecularMass(atomCounts)
	fmt.Println("molarMass", molarMass)
	// change parameters to test error messages
	fmt.Println("238U mass:", chemh.isotopeIsotopicMassBySym("U", 238))
	fmt.Println("238U abundance:", chemh.isotopeAbundanceBySym("U", 238))
}

func testParsingFormulaVerify1(fmA *formula) {
	var result formula
	result.components = append(result.components, formele{1, 0, 0, -1})
	result.components = append(result.components, formele{2, 0, 89, 0})
	result.components = append(result.components, formele{3, 0, 75, 1})
	result.components = append(result.components, formele{4, 10, 17, 2})
	result.components = append(result.components, formele{4, 17, 29, 2})
	result.components = append(result.components, formele{2, 0, 47, 0})
	result.components = append(result.components, formele{3, 0, 96, 5})
	result.components = append(result.components, formele{4, 1, 2, 6})
	result.components = append(result.components, formele{3, 0, 3, 6})
	result.components = append(result.components, formele{4, 16, 1, 8})
	result.components = append(result.components, formele{4, 8, 4, 8})
	testParsingFormulaCheckMatch(fmA, &result)
}

func testParsingFormulaVerify2(fmA *formula) {
	var result formula
	result.components = append(result.components, formele{1, 0, 0, -1})
	result.components = append(result.components, formele{2, 0, 1, 0})
	result.components = append(result.components, formele{3, 0, 1, 1})
	result.components = append(result.components, formele{4, 10, 1, 2})
	result.components = append(result.components, formele{4, 17, 1, 2})
	result.components = append(result.components, formele{2, 0, 1, 0})
	result.components = append(result.components, formele{3, 0, 1, 5})
	result.components = append(result.components, formele{4, 1, 2, 6})
	result.components = append(result.components, formele{3, 0, 3, 6})
	result.components = append(result.components, formele{4, 16, 1, 8})
	result.components = append(result.components, formele{4, 8, 4, 8})
	testParsingFormulaCheckMatch(fmA, &result)
}

func testParsingFormulaVerify3(fmA *formula) {
	var result formula
	result.components = append(result.components, formele{1, 0, 0, -1})
	result.components = append(result.components, formele{2, 0, 14, 0})
	result.components = append(result.components, formele{3, 0, 1, 1})
	result.components = append(result.components, formele{4, 10, 1, 2})
	result.components = append(result.components, formele{4, 17, 1, 2})
	result.components = append(result.components, formele{2, 0, 15, 0})
	result.components = append(result.components, formele{3, 0, 1, 5})
	result.components = append(result.components, formele{4, 1, 2, 6})
	result.components = append(result.components, formele{3, 0, 3, 5})
	result.components = append(result.components, formele{4, 16, 1, 8})
	result.components = append(result.components, formele{4, 8, 4, 8})
	testParsingFormulaCheckMatch(fmA, &result)
}

func testParsingFormulaVerify4(fmA *formula) {
	var result formula
	result.components = append(result.components, formele{1, 0, 0, -1})
	result.components = append(result.components, formele{2, 0, 1, 0})
	result.components = append(result.components, formele{3, 0, 1, 1})
	result.components = append(result.components, formele{4, 10, 1, 2})
	result.components = append(result.components, formele{4, 17, 1, 2})
	result.components = append(result.components, formele{2, 0, 1, 0})
	result.components = append(result.components, formele{3, 0, 1, 5})
	result.components = append(result.components, formele{4, 1, 2, 6})
	result.components = append(result.components, formele{3, 0, 3, 5})
	result.components = append(result.components, formele{4, 16, 1, 8})
	result.components = append(result.components, formele{4, 8, 4, 8})
	testParsingFormulaCheckMatch(fmA, &result)
}

func testParsingFormulaVerify5(fmA *formula) {
	var result formula
	result.components = append(result.components, formele{1, 0, 0, -1})
	result.components = append(result.components, formele{2, 0, 1, 0})
	result.components = append(result.components, formele{3, 0, 1, 1})
	result.components = append(result.components, formele{4, 10, 27, 2})
	result.components = append(result.components, formele{3, 0, 35, 1})
	result.components = append(result.components, formele{4, 17, 1, 4})
	result.components = append(result.components, formele{4, 10, 49, 4})
	result.components = append(result.components, formele{2, 0, 1, 0})
	result.components = append(result.components, formele{3, 0, 1, 7})
	result.components = append(result.components, formele{4, 1, 1, 8})
	result.components = append(result.components, formele{4, 8, 1, 8})
	testParsingFormulaCheckMatch(fmA, &result)
}

func testParsingFormulaVerify6(fmA *formula) {
	var result formula
	result.components = append(result.components, formele{1, 0, 0, -1})
	result.components = append(result.components, formele{2, 0, 1, 0})
	result.components = append(result.components, formele{3, 0, 2, 1})
	result.components = append(result.components, formele{4, 7, 1, 2})
	result.components = append(result.components, formele{4, 1, 4, 2})
	result.components = append(result.components, formele{3, 0, 1, 1})
	result.components = append(result.components, formele{4, 16, 1, 5})
	result.components = append(result.components, formele{4, 8, 4, 5})
	testParsingFormulaCheckMatch(fmA, &result)
}

func testParsingFormulaVerify7(fmA *formula) {
	var result formula
	result.components = append(result.components, formele{1, 0, 0, -1})
	result.components = append(result.components, formele{2, 0, 1, 0})
	result.components = append(result.components, formele{3, 0, 1, 1})
	result.components = append(result.components, formele{4, 20, 1, 2})
	result.components = append(result.components, formele{3, 0, 2, 1})
	result.components = append(result.components, formele{4, 7, 1, 4})
	result.components = append(result.components, formele{4, 8, 3, 4})
	testParsingFormulaCheckMatch(fmA, &result)
}

func testParsingFormula1() {
	var chemh chemHelper
	chemh.init()
	fmt.Println("BEGIN TEST 1")
	formulaStr := "89(Ne17Cl29)75 + 47(H2(S1O4)3)96"
	fmla := chemh.parseFormulaStrict(formulaStr)
	fmt.Print(formulaStr, " ")
	testParsingFormulaVerify1(&fmla)
	fmt.Println("fmla", fmla)
	fmt.Println(formulaStr)
	chemh.debugPrintlnFormula(&fmla)
	fmt.Print(formulaStr, ": (strict): ")
	chemh.printlnFormulaStrict(&fmla)
	fmt.Print(formulaStr, ": (normal): ")
	chemh.printlnFormula(&fmla)
	fmt.Println("END OF TEST 1")
}

func testParsingFormula2() {
	var chemh chemHelper
	chemh.init()
	fmt.Println("")
	fmt.Println("BEGIN TEST 2")
	formulaStr := "(NeCl) + (H2(S1O4)3)"
	fmla := chemh.parseFormula(formulaStr)
	fmt.Print(formulaStr, " ")
	testParsingFormulaVerify2(&fmla)
	fmt.Println("fmla", fmla)
	fmt.Println(formulaStr)
	chemh.debugPrintlnFormula(&fmla)
	fmt.Print(formulaStr, ": (strict): ")
	chemh.printlnFormulaStrict(&fmla)
	fmt.Print(formulaStr, ": (normal): ")
	chemh.printlnFormula(&fmla)
	fmt.Println("END OF TEST 2")
}

func testParsingFormula3() {
	var chemh chemHelper
	chemh.init()
	fmt.Println("")
	fmt.Println("BEGIN TEST 3")
	formulaStr := "14NeCl + 15H2(S1O4)3"
	fmla := chemh.parseFormula(formulaStr)
	fmt.Print(formulaStr, " ")
	testParsingFormulaVerify3(&fmla)
	fmt.Println("fmla", fmla)
	fmt.Println(formulaStr)
	chemh.debugPrintlnFormula(&fmla)
	fmt.Print(formulaStr, ": (strict): ")
	chemh.printlnFormulaStrict(&fmla)
	fmt.Print(formulaStr, ": (normal): ")
	chemh.printlnFormula(&fmla)
	fmt.Println("END OF TEST 3")
}

func testParsingFormula4() {
	var chemh chemHelper
	chemh.init()
	fmt.Println("")
	fmt.Println("BEGIN TEST 4")
	formulaStr := "NeCl + H2(SO4)3"
	fmla := chemh.parseFormula(formulaStr)
	fmt.Print(formulaStr, " ")
	testParsingFormulaVerify4(&fmla)
	fmt.Println("fmla", fmla)
	fmt.Println(formulaStr)
	chemh.debugPrintlnFormula(&fmla)
	fmt.Print(formulaStr, ": (strict): ")
	chemh.printlnFormulaStrict(&fmla)
	fmt.Print(formulaStr, ": (normal): ")
	chemh.printlnFormula(&fmla)
	fmt.Println("END OF TEST 4")
}

func testParsingFormula5() {
	var chemh chemHelper
	chemh.init()
	fmt.Println("")
	fmt.Println("BEGIN TEST 5")
	formulaStr := "Ne27(ClNe49)35 + HO"
	fmla := chemh.parseFormula(formulaStr)
	fmt.Print(formulaStr, " ")
	fmt.Println("fmla", fmla)
	testParsingFormulaVerify5(&fmla)
	fmt.Println("fmla", fmla)
	fmt.Println(formulaStr)
	chemh.debugPrintlnFormula(&fmla)
	fmt.Print(formulaStr, ": (strict): ")
	chemh.printlnFormulaStrict(&fmla)
	fmt.Print(formulaStr, ": (normal): ")
	chemh.printlnFormula(&fmla)
	fmt.Println("END OF TEST 5")
}

func testParsingFormula6() {
	var chemh chemHelper
	chemh.init()
	fmt.Println("")
	fmt.Println("BEGIN TEST 6")
	formulaStr := "(NH4)2SO4"
	// formulaStr := "(NH4)2(SO4)"
	fmla := chemh.parseFormula(formulaStr)
	fmt.Println("fmla", fmla)
	chemh.debugPrintlnFormula(&fmla)
	fmt.Print(formulaStr, " ")
	testParsingFormulaVerify6(&fmla)
	fmt.Println("fmla", fmla)
	fmt.Println(formulaStr)
	chemh.debugPrintlnFormula(&fmla)
	fmt.Print(formulaStr, ": (strict): ")
	chemh.printlnFormulaStrict(&fmla)
	fmt.Print(formulaStr, ": (normal): ")
	chemh.printlnFormula(&fmla)
	fmt.Println("END OF TEST 6")
}

func testParsingFormula7() {
	var chemh chemHelper
	chemh.init()
	fmt.Println("")
	fmt.Println("BEGIN TEST 7")
	formulaStr := "Ca(NO3)2"
	fmla := chemh.parseFormula(formulaStr)
	fmt.Print(formulaStr, " ")
	testParsingFormulaVerify7(&fmla)
	fmt.Println("fmla", fmla)
	fmt.Println(formulaStr)
	chemh.debugPrintlnFormula(&fmla)
	fmt.Print(formulaStr, ": (strict): ")
	chemh.printlnFormulaStrict(&fmla)
	fmt.Print(formulaStr, ": (normal): ")
	chemh.printlnFormula(&fmla)
	fmt.Println("END OF TEST 7")
}

func testMain() {
	testParsingFormulaMain()
	testParsingFormula1()
	testParsingFormula2()
	testParsingFormula3()
	testParsingFormula4()
	testParsingFormula5()
	testParsingFormula6()
	testParsingFormula7()
}

func appendix() {
	fmt.Println("B.1a 65.2 nanograms = 6.52e1 ng =", 6.52e1*nano, "g =", (6.52e1*nano)/kilo, "kg")
	fmt.Println("B.1b 88 picoseconds = 8.8e1 ps =", 8.8e1*pico, "s")
	fmt.Println("B.1c 5.4 terawatts =", 5.4*tera, "W =", 5.4*tera, "kg m^2 s^-3")
	fmt.Println("B.1d 17 kilovolts = 1.7e1 kV = ", 1.7e1*kilo, "V =", 1.7e1*kilo, "kg m^2 s^-3 A^-1")

	fmt.Println("B.2a 66 microK = 6.6e1 microK = ", 6.6e1*micro, "K")
	fmt.Println("B.2b 15.9 MJ = 1.59e1 MJ = ", 1.59e1*mega, "J =", 1.59e1*mega, "kg m^2 s^-2")
	fmt.Println("B.2c 0.13 mg = 1.3e-1 mg = ", 1.3e-1*milli, "g =", (1.3e-1*milli)/kilo, "kg")
	fmt.Println("B.2d 62 GPa = 6.2e1 GPa = ", 6.2e1*giga, "Pa =", (6.2e1*giga)/pascalsPerAtm, "Atm")

	var chemh chemHelper
	fmt.Println("B.3a 9001 F =", chemh.fahrenheitToCelsius(9001), "C")
	fmt.Println("B.3b 98.6 F =", chemh.fahrenheitToCelsius(98.6), "C")
	fmt.Println("B.3c 20 F above boiling point of water at 1 Atm =", chemh.fahrenheitToCelsius(chemh.celsiusToFahrenheit(boilingPointOfWaterAt1Atm)+20), "C")
	fmt.Println("B.3d -40 F =", chemh.fahrenheitToCelsius(-40), "C")

	fmt.Println("B.4a 5,000 C = ", chemh.celsiusToFahrenheit(5000), "F")
	fmt.Println("B.4b 40.0 C = ", chemh.celsiusToFahrenheit(40.0), "F")
	fmt.Println("B.4c 212 C = ", chemh.celsiusToFahrenheit(212), "F")
	fmt.Println("B.4d -40 C = ", chemh.celsiusToFahrenheit(-40), "F")

	fmt.Println("B.5a 9001 F =", chemh.fahrenheitToKelvin(9001), "K")
	fmt.Println("B.5b 98.6 F =", chemh.fahrenheitToKelvin(98.6), "K")
	fmt.Println("B.5c 20 F above boiling point of water at 1 Atm =", chemh.fahrenheitToKelvin(chemh.celsiusToFahrenheit(boilingPointOfWaterAt1Atm)+20), "K")
	fmt.Println("B.5d -40 F =", chemh.fahrenheitToKelvin(-40), "K")

	fmt.Println("B.6a 5,000 C = ", chemh.celsiusToKelvin(5000), "K")
	fmt.Println("B.6b 40.0 C = ", chemh.celsiusToKelvin(40.0), "K")
	fmt.Println("B.6c 212 C = ", chemh.celsiusToKelvin(212), "K")
	fmt.Println("B.6d -40 C = ", chemh.celsiusToKelvin(-40), "K")

	fmt.Println("B.7a 55.0 miles per hour =", (55.0*metersPerMile)/secondsPerHour, "m/s")
	fmt.Println("B.7b 1.15 g cm^-3 =", (1.15/kilo)/(centi*centi*centi), "kg m^-3")
	fmt.Println("B.7c 1.6e-19 C A (angstrom) =", 1.6e-19/(angstrom*angstrom), "s m WRONG")
	fmt.Println("B.7d 0.15 mol L^-1 =", 1.5e-1/m3perL, "mol m^-3")
	fmt.Println("B.7e 5.7e3 L atm day^-1 =", (5.7e3*m3perL*pascalsPerAtm)/secondsPerDay, "m^3 Pa s^-1")
	fmt.Println("B.8a 67.3 atm =", 67.3*pascalsPerAtm, "pascals =", 67.3*pascalsPerAtm, "kg m^-1s^-2 (= N m^2)")
	fmt.Println("B.8b 1.0 * 10^4 V cm^-1 =", 1e4, "kg m^2 s^-3 A^-1 =", (1e4 / centi), "kg m^1 s^-3 A^-1 (= N s^-1 A^-1")
	fmt.Println("B.8c 7.4 Ang year^-1 =", (7.4*angstrom)*secondsPerYear, "m s^-1")
	fmt.Println("B.8d 22.4 L mol^-1 =", 22.4*m3perL, "m^3 mol^-1")
	fmt.Println("B.8e 14.7 lb inch^-2 =", ((14.7*gramsPerPound)/kilo)/(mPerInch*mPerInch), "kg m^-2")

	fmt.Println("B.9 15.3 kWh =", 15.3*kilo*secondsPerHour, "kg m^2 s^-2 (= J)")
	fmt.Println("B.10 30 mpg (miles per gallon) =", (30.0*metersPerMile)/dm3perGallon, "m dm-3, conversion factor =", metersPerMile/dm3perGallon, ", 30 mpg =", 30.0*metersPerMile/m3perGallon, "m m^-3")
	fmt.Println("B.11 404 in^3 =", 404.0*(cmPerInch*cmPerInch*cmPerInch), "cm^3 and", 404.0*(cmPerInch*cmPerInch*cmPerInch)/cm3perL, "L")
	fmt.Println("B.12a 3.00 x 10^8 m s^-1 =", 3e8/metersPerMile, "miles per second")
	fmt.Println("B.12b 3.00 x 10^8 m s^-1 =", 3e8*(secondsPerDay*daysPerWeek*weeksPerFortnight)/(metersPerFoot*feetPerFurlong), "furlongs per fortnight")
}

func main() {
	appendix()
	var chemh chemHelper
	chemh.init()

	massOfOCompoundA := 0.22564
	fmt.Println("E1.1 Law of multiple proportions: Mass of O combined with 1.0000 g of Cl: Compound A: 0.22564 g, Compound B: 0.90255 g, Compound C: 1.3539 g, Compound D: 1.5795 g")
	fmt.Println("E1.1a Ratios to Compound A: Compound B", 0.90255/massOfOCompoundA, "Compound C", 1.3539/massOfOCompoundA, "Compound D", 1.5795/massOfOCompoundA)
	fmt.Println("E1.1b If Compound A is a multiple of Cl2O, then determine formulas for compounds B, C, and D: They're going to be Cl2O4 (and multiples: ClO2, Cl3O6, etc), Cl2O6 (and multiples: ClO3, Cl3O9, etc), and Cl2O7 (and multiples: Cl4O14, etc).")

	fmt.Println("E1.2 Number of electrons, protons, and neutrons in 222Rn (Radon-222):", chemh.atomicNumberBySym("Rn"), chemh.atomicNumberBySym("Rn"), 222-chemh.atomicNumberBySym("Rn"))

	// fmt.Println("E1.3 Relative atomic mass of carbon (calculated from isotopes)", (isos[isolookup[chemh.xyzSym("C"]].atomicNumber][12]].isotopicMass*isos[isolookup[chemh.xyzSym("C"]].atomicNumber][12]].abundance)+(isos[isolookup[chemh.xyzSym("C"]].atomicNumber][13]].isotopicMass*isos[isolookup[chemh.xyzSym("C"]].atomicNumber][13]].abundance))

	fmt.Println("E1.3 Relative atomic mass of carbon calculated from isotobes: 12C has an isotopic mass of 12.0000000 and abundance of 0.98892, and 13C has an isotopic mass of 13.003354 and abundance of 0.01108.")
	fmt.Println("    Multiplying out we get:", (chemh.isotopeIsotopicMassBySym("C", 12)*chemh.isotopeAbundanceBySym("C", 12))+(chemh.isotopeIsotopicMassBySym("C", 13)*chemh.isotopeAbundanceBySym("C", 13)))

	fmt.Println("E1.3 Relative atomic mass of carbon from table:", chemh.atomicMassBySym("C"))

	fmt.Println("Relative molecular mass of water:", (chemh.atomicMassBySym("H")*2)+chemh.atomicMassBySym("O"))

	fmt.Println("Mass of a single C atom:", chemh.atomicMassBySym("C")/avogadrosNumberPerMole, "g")
	fmt.Println("Mass of a single carbon-12 (12C) atom:", chemh.isotopeIsotopicMassBySym("C", 12)/avogadrosNumberPerMole, "g")

	fmt.Println("Ratio of mass of sodium atom (Na) to carbon atom (C)", chemh.atomicMassBySym("Na")/chemh.atomicMassBySym("C"))

	fmt.Println("E1.4 Mass of one uranium atom:", chemh.atomicMassBySym("U")/avogadrosNumberPerMole, "g")

	fmt.Println("    Mass of one uranium-238 isotope (238U) atom:", chemh.isotopeIsotopicMassBySym("U", 238)/avogadrosNumberPerMole, "g")

	fmt.Println("Chemical amount of iron in 8.232g:", 8.232/chemh.atomicMassBySym("Fe"), "mol")

	molarMassOfWater := (chemh.atomicMassBySym("H") * 2) + chemh.atomicMassBySym("O")
	fmt.Println("Amount of water needed for 0.2000 mol:", 0.2*molarMassOfWater, "g")

	molarMassOfNO2 := chemh.atomicMassBySym("N") + (chemh.atomicMassBySym("O") * 2)
	fmt.Println("E1.5a chemical amount of NO2 in 4.00 g of NO2:", 4.0/molarMassOfNO2, "mol")
	fmt.Println("E1.5b number of molecules of NO2 in 4.00 g of NO2:", (4.0/molarMassOfNO2)*avogadrosNumberPerMole)

	liquidBenzeneDensity := 0.8765 // in g/cm^3
	fmt.Println("Mass of 0.2124 L of liquid benzene:", (0.2124*cm3perL)*liquidBenzeneDensity, "g")
	molecularMassBenzene := (chemh.atomicMassBySym("C") * 6) + (chemh.atomicMassBySym("H") * 6)
	fmt.Println("Molecular mass of benzene:", molecularMassBenzene)
	fmt.Println("Chemical amount of benzene in 0.2124 L of benzene:", ((0.2124*cm3perL)*liquidBenzeneDensity)/molecularMassBenzene, "mol")

	iceDensityNear0C := 0.92 // in g/cm^3
	fmt.Println("Ice density near 0 C:", iceDensityNear0C, "g/cm^3")
	fmt.Println("Molar volume of solid water (ice):", molarMassOfWater/iceDensityNear0C, "cm^3/mol")

	oxygenDensityRoomTempSeaLevel := 0.00130 // in g/cm^3
	fmt.Println("Molar volume of oxygen (O2) at room temperature and pressure:", ((chemh.atomicMassBySym("O")*2)/oxygenDensityRoomTempSeaLevel)/cm3perL, "L/mol")

	fmt.Println("Volume per H2O molecule for solid water (ice):", (molarMassOfWater/avogadrosNumberPerMole)/iceDensityNear0C)

	fmt.Println("Molecular mass of water from formula:", chemh.determineMolecularMass(chemh.countFormulaAtoms(chemh.parseFormula("H2O"))))
	fmt.Println("Molecular mass of benzene from formula:", chemh.determineMolecularMass(chemh.countFormulaAtoms(chemh.parseFormula("C6H6"))))

	fmt.Println("P5. A sample of ascorbic acid (vitamin C) has 30.0 g of carbon, 40.0 g of oxygen, while aonther has 12.7 g of carbon. The mass of oxygen is:", 12.7*(40.0/30.0), "g")

	fmt.Println("P6. A sample compound has 25.0 g of hafnium and 31.5 g of tellurium. Another sample has 0.125 g of hafnium. The amount of tellurium in the second sample is:", 0.125*(31.5/25.0), "g")
	fmt.Println("P7. Compound 1 has 33.28% nitrogen (N) and 66.72% silicon (Si) by mass. Compound 2 has 39.94% N and 60.06% Si by mass.")
	fmt.Println("P7a. Amount of silicon that combines with 1.0000 g of nitrogen for Compound 1:", 66.72/33.28)
	fmt.Println("     Amount of silicon that combines with 1.0000 g of nitrogen for Compound 2:", 60.06/39.94)
	fmt.Println("P7b. Show these compounds satisfy the law of multiple proportions. The radio is:", ((66.72/33.28)/(60.06/39.94))*3, "/ 3")
	fmt.Println("     If the second has the formula Si3N4, the formula of the first compound is SiN (or a multiple)")

	fmt.Println("P8. Compound 1 has 86.979% iodine (I) and 13.021 % flourine (F) by mass. Compound 2 has 69.007% I and 30.993% F. Compound 3 has 57.191% I and 42.809% F. Compound 4 has 48.829% I and 51.171% F by mass. Compute mass of flourine that combines with 1.0000 g of iodine:")
	fmt.Println("     Compound 1:", 13.021/86.979)
	fmt.Println("     Compound 2:", 30.993/69.007)
	fmt.Println("     Compound 3:", 42.809/57.191)
	fmt.Println("     Compound 4:", 51.171/48.829)
	fmt.Println("P8b. Figure out the whole-number ratios for the previous four answers and show the compounds satisfy the law of multiple proportions")
	compound1 := 13.021 / 86.979
	fmt.Println("P8b. Compound 2:", (30.993/69.007)/compound1)
	fmt.Println("P8b. Compound 3:", (42.809/57.191)/compound1)
	fmt.Println("P8b. Compound 4:", (51.171/48.829)/compound1)

	fmt.Println("P9. Compound 1 is 76.10% vanadium (V) and 23.90% oxygen (O) by mass. Compound 2 is 67.98% V and 32.02% O. Compound 3 is 61.42% V and 38.58% O. Compound 4 is 56.02% V and 43.98% O. The relative numbers of atoms of oxygen in the compounds for a given mass of vanadium:")
	fmt.Println("    Compound 1:", ((23.90/chemh.atomicMassBySym("O"))/(76.10/chemh.atomicMassBySym("V")))*2, "/2")
	fmt.Println("    Compound 2:", ((32.02/chemh.atomicMassBySym("O"))/(67.98/chemh.atomicMassBySym("V")))*2, "/2")
	fmt.Println("    Compound 3:", ((38.58/chemh.atomicMassBySym("O"))/(61.42/chemh.atomicMassBySym("V")))*2, "/2")
	fmt.Println("    Compound 4:", ((43.98/chemh.atomicMassBySym("O"))/(56.02/chemh.atomicMassBySym("V")))*2, "/2")

	fmt.Println("P10. Compound 1 is 72.17% tungsten (W) and 27.83% chlorine (Cl) by mass. Compound 2 is 56.45% W and 43.55% Cl. Compound 3 is 50.91% W and 49.09% Cl. Compound 4 is 46.36% W and 53.64% Cl. Relative numbers of chlorine atoms for a given mass of tungsten:")
	fmt.Println("     Compound 1:", (27.83/chemh.atomicMassBySym("Cl"))/(72.17/chemh.atomicMassBySym("W")))
	fmt.Println("     Compound 2:", (43.55/chemh.atomicMassBySym("Cl"))/(56.45/chemh.atomicMassBySym("W")))
	fmt.Println("     Compound 3:", (49.09/chemh.atomicMassBySym("Cl"))/(50.91/chemh.atomicMassBySym("W")))
	fmt.Println("     Compound 4:", (53.64/chemh.atomicMassBySym("Cl"))/(46.36/chemh.atomicMassBySym("W")))
	fmt.Println("     If a molecule of each compound contains only one tungsten atom, the formulas for the four compounds are: Compound 1: WCl2, Compound 2: WCl4, Compound 3: WCl5, Compound 4: WCl6")

	fmt.Println("P11.  14.4 mL of hydrogen and 14.4 mL of oxygen from electrolysis of a liquid containing only hydrogen and oxygen.")
	fmt.Println("P11a. Chemical formula: HO because it's volume, not mass (14.4 mL)")
	fmt.Println("P11b. Reason more than one formula is possible: because it could be H2O2, H3O3, H4O4, etc")

	fmt.Println("P12. Liquid N2H4 decomposed into N2 and H2. N2 occupies 13.7 mL (room temperature and pressure). Volume of hydrogen under the same conditions:", 13.7*(4/2), "mL")

	fmt.Println("P13. Nitrogen dioxide (NO2) formed from dinitrogen oxide (N2O) and oxygen (O2) mixed (in the presence of a certain catalyst).")
	fmt.Println("     Formula: 2 N2O + 3 O2 -> 4 NO2")
	fmt.Println("     Volumes of N2O and oxygen needed to produce 4.0 L of NO2 if temperature and pressure held constant: (calculated by hand): 2.0 mL of N2O and 3.0 mL of O2 produce 4.0 L of NO2")
	fmt.Print("     Calculated again from equation using automatic equation balancer: ")
	nitrogenDioxideEqu := chemh.parseEquation("N2O + O2 -> NO2")
	chemh.balanceEquation(&nitrogenDioxideEqu)
	chemh.printlnEquation(&nitrogenDioxideEqu)

	fmt.Println("P14. Gaseous methanol (CH3OH) reacts with oxygen (O2) to produce water vapor and carbon dioxide. Volumes of water and carbon dioxide produced from 2.0 L of methanol if temperature and pressure held constant:")
	fmt.Println("     Equation from description: CH3OH + O2 -> H2O + CO2")
	fmt.Print("     balanced is: ")
	methanolEqu := chemh.parseEquation("CH3OH + O2 -> H2O + CO2")
	chemh.balanceEquation(&methanolEqu)
	chemh.printEquation(&methanolEqu)
	fmt.Println("")
	// balanced is: 2(CH3OH) + 3(O2) <-> 4(H2O) + 2(CO2)
	fmt.Println("     Therefore: the volume of water produced is 4.0 L and the volume of carbon dioxide is 2.0 L. (3.0 L of oxygen is also consumed in the process.)")

	fmt.Println("P15.  The isotope of plutonium used for nuclear fission is plutonium-239.")
	fmt.Println("P15a. Plutonium-239 ratio of neutrons to protons in the nucleus is", float64(239-chemh.atomicNumberBySym("Pu"))/float64(chemh.atomicNumberBySym("Pu")))
	fmt.Println("P15b. Plutonium-239 has", chemh.atomicNumberBySym("Pu"), "electrons")

	fmt.Println("P16.  The missing element in the first six periods, finally discovered in 1947 in the fission products from uranium, was called promethium.")
	fmt.Println("P16a. Promethium-145: ratio of the number of neutrons in a 145Pm nucleus to the number of protons:", float64(145-chemh.atomicNumberBySym("Pm"))/float64(chemh.atomicNumberBySym("Pm")))
	fmt.Println("P16b. Promethium-145: number of electrons:", chemh.atomicNumberBySym("Pm"))

	fmt.Println("P17. Americium-241 is used in smoke detectors. Number of protons:", chemh.atomicNumberBySym("Am"), ", number of neutrons:", 241-chemh.atomicNumberBySym("Am"), ", number of electrons:", chemh.atomicNumberBySym("Am"))

	fmt.Println("P18. In 1982 a single atom of 266 over 109 Mt (meitnerium-266) was reported. The atomic composition of this isotope is: number of protons:", chemh.atomicNumberBySym("Mt"), ", number of neutrons:", 266-chemh.atomicNumberBySym("Mt"), ", number of electrons:", chemh.atomicNumberBySym("Mt"))

	fmt.Println("P19. The natural abundances and isotopic masses of silicon (Si) relative to 12C = 12.00000 are:")
	fmt.Println("     28Si: Abundance 92.21%, Isotopic Mass: 27.97693")
	fmt.Println("     29Si: Abundance  4.70%, Isotopic Mass: 28.97649")
	fmt.Println("     30Si: Abundance  3.09%, Isotopic Mass: 29.97376")
	fmt.Println("     The atomic mass of naturally occurring silicon is:", (chemh.isotopeIsotopicMassBySym("Si", 28)*chemh.isotopeAbundanceBySym("Si", 28))+(chemh.isotopeIsotopicMassBySym("Si", 29)*chemh.isotopeAbundanceBySym("Si", 29))+(chemh.isotopeIsotopicMassBySym("Si", 30)*chemh.isotopeAbundanceBySym("Si", 30)))

	fmt.Println("P20. The natural abundances and isotopic masses of neon (Ne) are:")
	fmt.Println("     20Ne: Abundance 90.00%, Isotopic Mass: 19.99212")
	fmt.Println("     21Ne: Abundance  0.27%, Isotopic Mass: 20.99316")
	fmt.Println("     22Ne: Abundance  9.73%, Isotopic Mass: 21.99132")
	fmt.Println("     The atomic mass of naturally occurring neon is:", (chemh.isotopeIsotopicMassBySym("Ne", 20)*chemh.isotopeAbundanceBySym("Ne", 20))+(chemh.isotopeIsotopicMassBySym("Ne", 21)*chemh.isotopeAbundanceBySym("Ne", 21))+(chemh.isotopeIsotopicMassBySym("Ne", 22)*chemh.isotopeAbundanceBySym("Ne", 22)))

	fmt.Println("P21. Only 2 isotopes of boron (B) occur in nature. One, 10B has abundance 19.61% and atomic mass 10.013. The other, 11B, has abundance of 80.39%. The tabulated relative atomic mass of natural boron is 10.811. What is the atomic mass of 11B?")
	// (10BA * 10BM) + (11BA * 11BM) = TAB
	//                 (11BA * 11BM) = TAB - (10BA * 10BM)
	//                          11BM = (TAB - (10BA * 10BM)) / 11BA
	fmt.Println("     Is: ", (chemh.atomicMassBySym("B")-(chemh.isotopeIsotopicMassBySym("B", 10)*chemh.isotopeAbundanceBySym("B", 10)))/chemh.isotopeAbundanceBySym("B", 11))
	fmt.Println("P22. More than half of zirconium (Zr) atoms naturall occurring are 90Zr. All the other isotopes are:")
	fmt.Println("     91Zr: Abundance: 11.27%, Atomic Mass: 90.9056")
	fmt.Println("     92Zr: Abundance: 17.17%, Atomic Mass: 91.9050")
	fmt.Println("     94Zr: Abundance: 17.33%, Atomic Mass: 92.9063")
	fmt.Println("     96Zr: Abundance:  2.78%, Atomic Mass: 93.9083")
	abundance90Zr := 1.00 - chemh.isotopeAbundanceBySym("Zr", 91) - chemh.isotopeAbundanceBySym("Zr", 92) - chemh.isotopeAbundanceBySym("Zr", 94) - chemh.isotopeAbundanceBySym("Zr", 96)
	totalOtherZr := (chemh.isotopeIsotopicMassBySym("Zr", 91) * chemh.isotopeAbundanceBySym("Zr", 91)) + (chemh.isotopeIsotopicMassBySym("Zr", 92) * chemh.isotopeAbundanceBySym("Zr", 92)) + (chemh.isotopeIsotopicMassBySym("Zr", 94) * chemh.isotopeAbundanceBySym("Zr", 94)) + (chemh.isotopeIsotopicMassBySym("Zr", 96) * chemh.isotopeAbundanceBySym("Zr", 96))
	fmt.Println("     Relative atomic mass of 90Zr is:", (chemh.atomicMassBySym("Zr")-totalOtherZr)/abundance90Zr)

	fmt.Println("P23. If the relative atomic mass of iodine is 126.90447, the mass of one iodine atom is:", chemh.atomicMassBySym("I")/avogadrosNumberPerMole, "g")
	fmt.Println("P24. The mass of 100 million flourine atoms if the relative atomic mass of flourine is 18.998403 on a scale on which exactly 12 is the relative atomic mass of 12C is:", (1e8*chemh.atomicMassBySym("F"))/avogadrosNumberPerMole, "g")

	fmt.Println("P25. Relative atomic masses on 12C scale of the following compounds:")
	fmt.Println("     a. P4O10", chemh.determineMolecularMass(chemh.countFormulaAtoms(chemh.parseFormula("P4O10"))))
	fmt.Println("     b. BrCl:", chemh.determineMolecularMass(chemh.countFormulaAtoms(chemh.parseFormula("BrCl"))))
	fmt.Println("     c. Ca(NO3)2", chemh.determineMolecularMass(chemh.countFormulaAtoms(chemh.parseFormula("Ca(NO3)2"))))
	fmt.Println("     d. KMnO4:", chemh.determineMolecularMass(chemh.countFormulaAtoms(chemh.parseFormula("KMnO4"))))
	fmt.Println("     e. (NH4)2SO4", chemh.determineMolecularMass(chemh.countFormulaAtoms(chemh.parseFormula("(NH4)2SO4"))))

	fmt.Println("P26. Relative atomic masses on 12C scale of the following compounds:")
	fmt.Println("     a. (Ag(NH3)2)Cl", chemh.determineMolecularMassForFormula("(Ag(NH3)2)Cl"))
	fmt.Println("     b. Ca3(Co(CO3)3)2:", chemh.determineMolecularMassForFormula("Ca3(Co(CO3)3)2"))
	fmt.Println("     c. OsO4", chemh.determineMolecularMassForFormula("OsO4"))
	fmt.Println("     d. H2SO4:", chemh.determineMolecularMassForFormula("H2SO4"))
	fmt.Println("     e. Ca3Al2(SiO4)3", chemh.determineMolecularMassForFormula("Ca3Al2(SiO4)3"))

	fmt.Println("P27. Suppose a person counts out gold atoms at the rate of one each second for the entire span of an 80-year life. (Even when he/she sleeps?)")
	fmt.Println("     That's: ", ((secondsPerYear*80)/avogadrosNumberPerMole)*chemh.atomicMassBySym("Au"), "g")
	fmt.Println("     which is", (((secondsPerYear*80)/avogadrosNumberPerMole)*chemh.atomicMassBySym("Au"))/femto, "femtograms")

	fmt.Println("P28. A gold atom has diameter 2.88 x 10^-10 m. Suppose the atoms in 1.00 mol of gold atoms are arrarged just touching their neighbors in a single-file line.")
	fmt.Println("     The length of the line is", 2.88e-10*(chemh.atomicMassBySym("Au")*avogadrosNumberPerMole), "m")
	fmt.Println("     which is", (2.88e-10*(chemh.atomicMassBySym("Au")*avogadrosNumberPerMole))/peta, "petameters")
	fmt.Println("     which is", (2.88e-10*(chemh.atomicMassBySym("Au")*avogadrosNumberPerMole))/metersPerLightYear, "light years.")

	fmt.Println("P29. The Vitamin A molecule has the formula C20H30O, and A2 has the formula C20H28O. Determine how many moles of vitamin A2 contain the same number of atoms as 1.000 mol of vitamin A.")
	vitA := chemh.countFormulaAtoms(chemh.parseFormula("C20H30O"))
	vitA2 := chemh.countFormulaAtoms(chemh.parseFormula("C20H28O"))
	fmt.Println("answer", float64(vitA[1]+vitA[6]+vitA[8])/float64(vitA2[1]+vitA2[6]+vitA2[8]))

	fmt.Println("P30. Arrange the following in order of increasing mass: 1.06 mol SF4; 117g CH4; 8.7 x 10^23 molecules of Cl2O7; 417 x 10^23 atoms of argon (AR).")
}
