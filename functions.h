#ifndef FUNCTIONS_H
#define FUNCTIONS_H 1

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>
#include <chrono>
#include <sys/stat.h>
#include "gzstream.h"

using namespace std;
using namespace chrono;

void writeAllSequencesToFile(int n, int S, int** x, string file_name);

void writeAllSequencesToStream(int n, int S, int** x, ogzstream& ostream);

void readOneSequenceFromFile(int n, int seq_idx, int* x, string file_name);

void readOneSequenceFromStream(int n, int seq_idx, int* x, igzstream& istream);

void readOneSequenceFromStreamAtCurPos(int n, int* x, igzstream& istream);

void readAllSequencesFromFile(int n, int S, int** x, string file_name);

void readAllSequencesFromStream(int n, int S, int** x, igzstream& istream);

#endif // FUNCTIONS_H
