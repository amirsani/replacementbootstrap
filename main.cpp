#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "functions.h"
#include "bootstrap_methods.h"
#include <math.h>

using namespace std;

int main(int argc,char** argv)
{
    string	input_file = argv[1];

    int bootstraps = atoi(argv[2]),
        n = atoi(argv[3]),
        replacements = atoi(argv[4]),
        change_option = 1,
        x[n];

    const int A = atoi(argv[5]);

    readOneSequenceFromFile(n, 0, x, input_file);

    RBoot(n, A, x, replacements, bootstraps, floor(1.5*log(n)), change_option);

    return 0;
}
