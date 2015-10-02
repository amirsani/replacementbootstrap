#include <iostream>
#include <cmath>
#include "functions.h"

using namespace std;
using namespace chrono;

void writeAllSequencesToFile (
		int n,
		int S,
		int** x,
		string file_name)
{
	ogzstream outfile(file_name.c_str(), ios::out | ios::binary);
	writeAllSequencesToStream(n, S, x, outfile);
	outfile.close();
}

void writeAllSequencesToStream (
		int n,
		int S,
		int** x,
		ogzstream& ostream)
{
	for (int i=0; i<S; i++)
		for (int t=0; t<n; t++)
			ostream.write((char*)&x[i][t], sizeof(x[i][t]));
}

void readOneSequenceFromFile (
		int n,
		int seq_idx,
		int* x,
		string file_name)
{
	igzstream infile;
	infile.open(file_name.c_str(), ios::in | ios::binary);
	if (!infile.good())
	{
		cerr << "ERROR: File " << file_name << " doesn't exist..." << endl;
		exit(1);
	}
	readOneSequenceFromStream(n, seq_idx, x, infile);
	infile.close();
}

void readOneSequenceFromStream (
		int n,
		int seq_idx,
		int* x,
		igzstream& istream)
{
	for (int i=0; i<=seq_idx; i++)
	{
		for (int t=0; t<n; t++)
		{
			if (istream.eof( ))
			{
				cerr << "ERROR: Cannot access element [" << i << "; " << t << "], the file is over..." << endl;
				exit(1);
			}
			istream.read((char*)&x[t], sizeof(x[t]));
		}
	}
}

void readOneSequenceFromStreamAtCurPos (
		int n,
		int* x,
		igzstream& istream)
{
	for (int t=0; t<n; t++)
	{
		if (istream.eof( ))
		{
			cerr << "ERROR: Cannot access element [" << t << "], the file is over..." << endl;
			exit(1);
		}
		istream.read((char*)&x[t], sizeof(x[t]));
	}
}

void readAllSequencesFromFile (
		int n,
		int S,
		int** x,
		string file_name)
{
	igzstream infile;
	infile.open(file_name.c_str(), ios::in | ios::binary);
	if (!infile.good())
	{
		cerr << "ERROR: File " << file_name << " doesn't exist..." << endl;
		exit(1);
	}

	readAllSequencesFromStream(n, S, x, infile);
	infile.close();
}

void readAllSequencesFromStream (
		int n,
		int S,
		int** x,
		igzstream& istream)
{
	for (int i=0; i<S; i++)
		for (int t=0; t<n; t++)
		{
			if (istream.eof( ))
			{
				cerr << "ERROR: Cannot access element [" << i << "; " << t << "], the file is over..." << endl;
				exit(1);
			}
			istream.read((char*)&x[i][t], sizeof(x[i][t]));
		}
}
