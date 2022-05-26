/*
* Identifies possible AGE crosslink candidates from the .dcd trajectory file
* mahditavakol90@gmail.com
*/

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include "Atoms.h"
#include "XST.h"
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
static const std::string vmd = "start vmd";
#else
static const std::string vmd = "vmd";
#endif
using namespace std;
string tempFolder("temp/");

void inputFromUser(string& pdbFileName, string& dcdFileName, string& xstFileName,
	string& dcdFreqstr, string& distanceStr, string& numFrameStr, string& outputFileName);
void vmdFrameExtractor(const string& pdbFileName, const string& dcdFileName,
	const string& tempFolder, bool runVMD);
void readXST(const string& xstFileName, vector<XST>& myXST, bool runMode);
XST extrapolateXST(const int& step, const vector<XST>& myXST, bool runMode);
void addPBCImageX(const int& numberOfFrames, const string& xstFileName,
	const vector<XST>& myXST, const int& dcdFreq, const int& distance,
	bool runMode);
void addPBCImageY(const int& numberOfFrames, const string& xstFileName,
	const vector<XST>& myXST, const int& dcdFreq,
	const int& distance, bool runMode);
int checkPDBLine(const string& line);
int checkGeometery(const float* coordinates, const float* a, const float* b, const float& distance);
void outputAGEInfo(const string& modefiedPDBFileNames, const string& outputFileName,
	const int& numberOfFrames, const float& distance, bool runVMD);
void TestOutputAGEInfo();
void TestvmdFrameExtractor();

int main(int argc, char** argv)
{
	string pdbFileName, dcdFileName, xstFileName, outputFileName, dcdFreqStr, distanceStr, numFrameStr;
	string modifiedPDBFileNames;
	int dcdFreq, numOfFrames;
	float distance;
	vector<XST> xst;
	if (argc == 15)
	{
		for (int i = 0; i < argc / 2; i++)
		{
			if      (string(argv[2 * i + 1]) == "-pdb")     pdbFileName    = string(argv[2 * i + 2]);
			else if (string(argv[2 * i + 1]) == "-dcd")     dcdFileName    = string(argv[2 * i + 2]);
			else if (string(argv[2 * i + 1]) == "-xst")     xstFileName    = string(argv[2 * i + 2]);
			else if (string(argv[2 * i + 1]) == "-dcdfrq")  dcdFreqStr     = string(argv[2 * i + 2]);
			else if (string(argv[2 * i + 1]) == "-agecut")  distanceStr    = string(argv[2 * i + 2]);
			else if (string(argv[2 * i + 1]) == "-numfrm")  numFrameStr    = string(argv[2 * i + 2]);
			else if (string(argv[2 * i + 1]) == "-output")  outputFileName = string(argv[2 * i + 2]);
			else
			{
				inputFromUser(pdbFileName, dcdFileName, xstFileName, dcdFreqStr,distanceStr,numFrameStr, outputFileName);
			}
		}
	}
	else
	{
		inputFromUser(pdbFileName, dcdFileName, xstFileName, dcdFreqStr,distanceStr,numFrameStr, outputFileName);
	}
	stringstream iss(dcdFreqStr);
	stringstream isb(distanceStr);
	stringstream isa(numFrameStr);
	iss >> dcdFreq;
	isb >> distance;
	isa >> numOfFrames;
	modifiedPDBFileNames = "-PBCX-PBCY";

	vmdFrameExtractor(pdbFileName, dcdFileName,tempFolder, true);
	readXST(xstFileName, xst, true);
	string cdScript("cd ");
	cdScript += tempFolder;
	//system(cdScript.c_str());
	addPBCImageX(numOfFrames, xstFileName, xst, dcdFreq, distance, true);
	addPBCImageY(numOfFrames, xstFileName, xst, dcdFreq, distance, true);
	outputAGEInfo(modifiedPDBFileNames, outputFileName,numOfFrames,distance, true);
	//TestOutputAGEInfo();
	//TestvmdFrameExtractor();
	//readXST(xstFileName, xst, true);
	//addPBCImageX(numOfFrames, xstFileName, xst, dcdFreq, distance, true);
	//addPBCImageY(numOfFrames, xstFileName, xst, dcdFreq, distance, true);
	//int step = 75000;
	//extrapolateXST(step, xst, false);

	return 0;
}

void inputFromUser(string& pdbFileName, string& dcdFileName, string& xstFileName,
	string& dcdFreqStr, string& distanceStr, string& numFrameStr, string& outputFileName)
{
	cout << "Enter the pdb file name: " << endl;
	cin >> pdbFileName;
	cout << "Enter the dcd file name: " << endl;
	cin >> dcdFileName;
	cout << "Enter the xst file name: " << endl;
	cin >> xstFileName;
	cout << "Enter the dcd output frequency: " << endl;
	cin >> dcdFreqStr;
	cout << "Enter number of Frames " << endl;
	cin >> numFrameStr;
	cout << "Enter AGE cutoff: " << endl;
	cin >> distanceStr;
	cout << "Enter Output File Name: " << endl;
	cin >> outputFileName;
}

void vmdFrameExtractor(const string& pdbFileName, const string& dcdFileName, const string& tempFolder, bool runVmd)
{
	string tmpFileName("frameExtractor.tcl");
	ofstream tclFrameExtractor(tmpFileName);
	string flderScrpt("mkdir ");
	flderScrpt += tempFolder;
	system(flderScrpt.c_str());
	tclFrameExtractor << "set mol [mol new " << pdbFileName << ".pdb type pdb waitfor all]" << endl;
	tclFrameExtractor << "mol addfile " << dcdFileName << ".dcd type dcd waitfor all molid $mol" << endl;
	tclFrameExtractor << "cd " << tempFolder << endl;
	tclFrameExtractor << "set numFrames [molinfo top get numframes]" << endl;
	tclFrameExtractor << "set firstFrame 0" << endl;
	tclFrameExtractor << "pbc wrap -first $firstFrame -last [expr $numFrames - 1]" << endl;
	tclFrameExtractor << "for {set frame $firstFrame} {$frame < $numFrames } {incr frame} {" << endl;
	tclFrameExtractor << "	set sel [atomselect top all frame $frame]" << endl;
	tclFrameExtractor << "  pbc wrap -first $frame -last $frame" << endl;
	tclFrameExtractor << "  $sel writepdb $frame.pdb" << endl;
	tclFrameExtractor << "}" << endl;
	tclFrameExtractor << "exit" << endl;
	string runScrpt("vmd -dispdev text -eofexit -e ");
	string numFrame("ls -1 temp | wc -l");
	string delScrpt("rm -rf temp");
	string dleScrpt("rm ");
	runScrpt = runScrpt + tmpFileName;
	delScrpt = delScrpt;
	dleScrpt = dleScrpt + tmpFileName;
	if (runVmd)
	{
		system(runScrpt.c_str());
		//system(numFrame.c_str());
		//system(delScrpt.c_str());
		//system(dleScrpt.c_str());
	}
}


void readXST(const string& xstFileName, vector<XST>& myXST, bool runMode)
{
	string fileName = xstFileName + ".xst";
	ifstream file(fileName);
	string line;
	int step;
	float a[3], b[3], c[3], o[3];
	while (getline(file, line))
	{
		if (line.substr(0, 1) == "#") continue; // Skipping the comments section
		XST newXST;
		stringstream iss(line);
		iss >> step;
		iss >> a[0];
		iss >> a[1];
		iss >> a[2];
		iss >> b[0];
		iss >> b[1];
		iss >> b[2];
		iss >> c[0];
		iss >> c[1];
		iss >> c[2];
		iss >> o[0];
		iss >> o[1];
		iss >> o[2];
		newXST.setStep(step);
		newXST.setA(a);
		newXST.setB(b);
		newXST.setC(c);
		newXST.setO(o);
		myXST.push_back(newXST);
	}
	if (runMode == false)
	{
		ofstream file("TestReadXST.txt");
		for (int i = 0; i < myXST.size(); i++)
		{
			int step;
			float a[3], b[3], c[3], o[3];
			myXST[i].getStep(step);
			myXST[i].getA(a);
			myXST[i].getB(b);
			myXST[i].getC(c);
			myXST[i].getO(o);
			file << step << " ";
			file << a[0] << " ";
			file << a[1] << " ";
			file << a[2] << " ";
			file << b[0] << " ";
			file << b[1] << " ";
			file << b[2] << " ";
			file << c[0] << " ";
			file << c[1] << " ";
			file << c[2] << " ";
			file << o[0] << " ";
			file << o[1] << " ";
			file << o[2] << " ";
			file << endl;
		}
	}
}

XST extrapolateXST(const int& step, const vector<XST>& myXST, bool runMode)
{
	XST newXST;
	int n = myXST.size();
	int i = 0;
	float a[3], b[3], c[3], o[3];
	for ( ; i < n; i++)
	{
		int myStep;
		XST xstI = myXST[i];
		xstI.getStep(myStep);
		if (step == myStep)
		{
			if (runMode == false)
			{
				xstI.getA(a);
				xstI.getB(b);
				xstI.getC(c);
				xstI.getO(o);
				cout << "step is:" << endl;
				cout << step << endl;
				cout << "a is: " << endl;
				cout << a[0] << "," << a[1] << "," << a[2] << endl;
				cout << "b is: " << endl;
				cout << b[0] << "," << b[1] << "," << b[2] << endl;
				cout << "c is: " << endl;
				cout << c[0] << "," << c[1] << "," << c[2] << endl;
				cout << "o is: " << endl;
				cout << o[0] << "," << o[1] << "," << o[2] << endl;
			}
			return xstI;
		}
		else if (step > myStep) continue;
		else if (step < myStep) break;
	}
	if (i == n)
	{
		cout << "Error in extrapolateXST !!!" << endl;
		return myXST[i - 1];
	}
	int stepI, stepII;
	float aI[3], bI[3], cI[3], oI[3];
	float aII[3], bII[3], cII[3], oII[3];
	XST xstI  = myXST[i - 1];
	XST xstII = myXST[i];
	xstI.getStep(stepI);
	xstII.getStep(stepII);
	xstI.getA(aI);
	xstII.getA(aII);
	xstI.getB(bI);
	xstII.getB(bII);
	xstI.getC(cI);
	xstII.getC(cII);
	xstI.getO(oI);
	xstII.getO(oII);
	float diff = (float)(step - stepI) / (float)(stepII - stepI);
	for (int j = 0; j < 3; j++)
	{
		a[j] = aI[j] + diff * (aII[j] - aI[j]);
		b[j] = bI[j] + diff * (bII[j] - bI[j]);
		c[j] = cI[j] + diff * (cII[j] - cI[j]);
		o[j] = oI[j] + diff * (oII[j] - oI[j]);
	}
	newXST.setStep(step);
	newXST.setA(a);
	newXST.setB(b);
	newXST.setC(c);
	newXST.setO(o);
	if (runMode) return newXST;
	else
	{
		newXST.getA(aI);
		newXST.getB(bI);
		newXST.getC(cI);
		newXST.getO(oI);
		cout << "step is:" << endl;
		cout << step << endl;
		cout << "a is: " << endl;
		cout << aI[0] << "," << aI[1] << "," << aI[2] << endl;
		cout << "b is: " << endl;
		cout << bI[0] << "," << bI[1] << "," << bI[2] << endl;
		cout << "c is: " << endl;
		cout << cI[0] << "," << cI[1] << "," << cI[2] << endl;
		cout << "o is: " << endl;
		cout << oI[0] << "," << oI[1] << "," << oI[2] << endl;
	}
	return newXST;
}

void addPBCImageX(const int& numberOfFrames, const string& xstFileName,
	const vector<XST>& myXST,  const int& dcdFreq, const int& distance,
	bool runMode)
{
	string line;
	for (int i = 0; i < numberOfFrames; i++)
	{
		vector<Atom> atoms;
		int resID = 0;
		int atomID = 0;
		string commentLine("AGEIdentifier written by Mahdi Tavakol, mahditavakol90@gmail.com");
		stringstream ss;
		ss << i;
		string fileName = tempFolder + ss.str() + ".pdb";
		string outputFileName = tempFolder + ss.str() + "-PBCX.pdb";
		ofstream outputFile(outputFileName);
		Atom dummAtom;
		dummAtom.PrintPDBHeader(outputFile, commentLine);
		ifstream pdbFrame(fileName);
		int step = i*dcdFreq;
		while (getline(pdbFrame, line))
		{
			int success = checkPDBLine(line);
			if (success == 3) break;
			else if (success == 2) continue;
			else if (success == 0) break; // This should never happen.
			else if (success == 1)
			{
				atomID++;
				float a[3], b[3], c[3], o[3];
				float coordinates[3];
				XST newXST;
				Atom newAtom;
				newAtom.ReadPDBLine(line);
				//newAtom.SetResidueID(resID);
				newXST = extrapolateXST(step, myXST,true);
				newXST.getA(a);
				newXST.getB(b);
				newXST.getC(c);
				newXST.getO(o);
				newAtom.GetCoordinates(coordinates);
				//coordinates[0] -= o[0];
				//coordinates[1] -= o[1];
				//coordinates[2] -= o[2];
				newAtom.SetAtomID(atomID);
				newAtom.SetBeta(0.0f);
				newAtom.PrintAtomsPDB(outputFile);
				int testGeo = checkGeometery(coordinates, a, b, distance);
				if (testGeo == 2)
				{
					atomID++;
					Atom imageAtom(newAtom);
					imageAtom.SetAtomID(atomID);
					coordinates[0] += a[0];
					imageAtom.SetCoordinates(coordinates);
					imageAtom.SetBeta(1.0f);
					imageAtom.PrintAtomsPDB(outputFile);
				}

			}
		}
		outputFile << "END" << endl;
	}
}

void addPBCImageY(const int& numberOfFrames, const string& xstFileName,
				  const vector<XST>& myXST, const int& dcdFreq,
	              const int& distance, bool runMode)
{
	string line;
	for (int i = 0; i < numberOfFrames; i++)
	{
		vector<Atom> atoms;
		int atomID = 0;
		string commentLine("AGEIdentifier written by Mahdi Tavakol, mahditavakol90@gmail.com");
		stringstream ss;
		ss << i;
		string fileName = tempFolder + ss.str() + "-PBCX.pdb";
		string outputFileName =  tempFolder + ss.str() + "-PBCX-PBCY.pdb";
		ofstream outputFile(outputFileName);
		Atom dummAtom;
		dummAtom.PrintPDBHeader(outputFile, commentLine);
		ifstream pdbFrame(fileName);
		int step = i*dcdFreq;
		while (getline(pdbFrame, line))
		{
			int success = checkPDBLine(line);
			if (success == 3) break;
			else if (success == 2) continue;
			else if (success == 0) break; // This should never happen.
			else if (success == 1)
			{
				atomID++;
				float a[3], b[3], c[3], o[3];
				float coordinates[3];
				XST newXST;
				Atom newAtom;
				newAtom.ReadPDBLine(line);
				newXST = extrapolateXST(step, myXST, true);
				newXST.getA(a);
				newXST.getB(b);
				newXST.getC(c);
				newXST.getO(o);
				newAtom.GetCoordinates(coordinates);
				//coordinates[0] -= o[0];
				//coordinates[1] -= o[1];
				//coordinates[2] -= o[2];
				newAtom.SetAtomID(atomID);
				newAtom.SetBeta(0.0f);
				newAtom.PrintAtomsPDB(outputFile);
				int testGeo = checkGeometery(coordinates, a, b, distance);
				if (testGeo == 1)
				{
					atomID++;
					Atom imageAtom(newAtom);
					imageAtom.SetAtomID(atomID);
					coordinates[0] += b[0];
					coordinates[1] += b[1];
					imageAtom.SetCoordinates(coordinates);
					imageAtom.SetBeta(1.0f);
					imageAtom.PrintAtomsPDB(outputFile);
				}

			}
		}
		outputFile << "END" << endl;
	}
}

int checkPDBLine(const string& line)
{
	if (line.substr(0, 3) == "END")  return 3;
	else if (line.substr(0, 4) != "ATOM") return 2;
	else if (line.substr(0, 4) == "ATOM") return 1;
	return 0;
}

int checkGeometery(const float* coordinates, const float* a, const float* b, const float& distance)
{
	float x     = coordinates[0];
	float y     = coordinates[1];
	float z     = coordinates[2];
	float a1    = a[0];
	float a2    = a[1];
	float b1    = b[0];
	float b2    = b[1];
	float left  = (b2 / b1)*x - (distance * sqrt(b1*b1 + b2*b2) / b1);
	float right = (b2 / b1)*x;

	if (y >= 0 && y <= distance)
	{
		return 1;
	}

	if (y >= left && y <= right)
	{
		return 2;
	}
	return 0;
}


void outputAGEInfo(const string& modifiedPDBFileNames, const string& outputFileName, const int& numberOfFrames, const float& distance, bool shouldIrunVMD)
{
	string fileName("AGEIden.tcl");
	string scrptfileName("AGEIdentifier.tcl");
	ofstream file(fileName);
	ofstream file2(scrptfileName);
	file << "proc IdentifyAGE {distance pdbFileName outputFileName} {" << endl;
	file << "	set tempFile1 [open temp1.tmp a+]" << endl;
	file << "	set tempFile2 [open temp2.tmp a+]" << endl;
	file << "	set mol [mol new $pdbFileName.pdb type pdb waitfor all]" << endl;
	file << "	set AGEArg [atomselect top \"(resname ARG and (type NH1 or type NH2)) and within $distance of(resname LYS and type NZ)\"]" << endl;
	file << "	foreach ArgID[$AGEArg get resid] ArgChain[$AGEArg get chain]{" << endl;
	file << "		set AGELys[atomselect top \"(resname LYS and type NZ) and within $distance of (resid $ArgID)\"]" << endl;
	file << "		foreach LysID[$AGELys get resid] LysChain[$AGELys get chain]{" << endl;
	file << "			puts - nonewline $tempFile1 \"ARG$ArgID$ArgChain - LYS$LysID$LysChain,\"" << endl;
	file << "			puts - nonewline $tempFile2 \"ARG$ArgID$ArgChain-LYS$LysID$LysChain,\"" << endl;
	file << "		}" << endl;
	file << "	}" << endl;
	file << "	puts $tempFile2 \" \"" << endl;
	file << "}" << endl;


	file2 << "set completeFileName $outputFileName.csv" << endl;
	file2 << "set csvfile   [open $completeFileName w]" << endl;
	file2 << "source " << fileName << endl;
	for (int i = 0; i < numberOfFrames; i++)
	{
		file2 << "IdentifyAGE " << distance << " " << tempFolder << i << modifiedPDBFileNames << " " << outputFileName << endl;
	}
	file2 << endl << endl << endl;
	file2 << "set tempFile1 [open \"temp1.tmp\" r]" << endl <<
		"set tempFile2 [open \"temp2.tmp\" r]" << endl <<
		"set tempLine [read $tempFile1]" << endl <<
		"set keys     [split $tempLine \",\"]" << endl <<
		"set uniqueKeys [lsort - unique $keys]" << endl;

	file2 << "puts -nonewline $file \"frame:, \"" << endl <<
		"foreach key $uniqueKeys{" << endl <<
		"	if {$key != {}} {" << endl <<
		"		puts - nonewline $file \"$key,\"" << endl <<
		"	}" << endl <<
		"}" << endl <<
		"puts $file \" \"" << endl << endl <<
		"while {[gets $tempFile2 line] >= 0} {" << endl <<
		"	set lineKeys [split $line \",\"]" << endl <<
		"	set frame [lindex $lineKeys 0]" << endl <<
		"	puts -nonewline $file \"$frame,\"" << endl <<
		"	foreach key $uniqueKeys{" << endl <<
		"		if {[lsearch -exact $lineKeys $key] >= 0 } {" << endl <<
		"			#puts -nonewline $file \"$key,\"" << endl <<
		"			puts -nonewline $file \"1,\"" << endl <<
		"		} elseif{ $key != {} } {" << endl <<
		"		puts -nonewline $file \"0,\"" << endl <<
		"	}" << endl <<
		"	puts $file \" \"" << endl <<
		"}" << endl <<
		"close $file" << endl <<
		"exit " << endl;
	if (shouldIrunVMD)
	{
		string runScrpt("vmd -dispdev text -eofexit -e ");
		runScrpt += scrptfileName;
		system(runScrpt.c_str());
	}
}

void TestOutputAGEInfo()
{
	outputAGEInfo("-PBCX-PBCY", "AGE", 100, 5, false);
}

void TestvmdFrameExtractor()
{
	vmdFrameExtractor("khar", "kekeh", "temp", false);
}
