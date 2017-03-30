//wham.cpp

/*

The modified version of weighted histogram analysis method (WHAM)

Calculating PMF for each window from umbrella sampling data,
and calculate the weight for each frame.
(see reference (1) eq.(8), and )

Free energy (PMF) and weight are related by following equation:
(free energy) = -1 * ln(weight) / k
where k is Boltzman constant

Reference:
(1) Moradi, M., Enkavi, G., & Tajkhorshid, E. (2015). 
Atomic-level characterization of transport cycle thermodynamics in the glycerol-3-phosphate:phosphate antiporter. 
Nature Communications, 6, 8393. http://doi.org/10.1038/ncomms9393
eq.(7), eq(8) and eq.(9)
w^t corresponds to the weight

(2) Bartels, C. (2000). 
Analyzing biased Monte Carlo and molecular dynamics simulations. 
Chemical Physics Letters, 331(5–6), 446–454. http://doi.org/10.1016/S0009-2614(00)01215-X
rho_{i, t} is corresponding to the weight in this code.


usage: input parameters as stdin

>make
>./weight < param.txt

Parameters:
TOLERANCE
N_REPLICA
N_DATA ... number of time points in the input file
INPUT_PATTERN ... input file name pattern
WEIGHT_OUT ... filename for output weights
TEMPERATURE
MAX_ITER
CENTERS ... list of umbrella centers
SPRINGS ... list of spring constans of umbrella potential

input file should have following format:
[n_step] [value]\n
...

*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <cmath>
using std::cout;
using std::endl;
using std::cin;
using std::string;
using std::ifstream;
using std::ofstream;


double bias(double current, double center, double spring);
//double calcRMSD(double old[], double weight[]);


int main(int argc, char const *argv[])
{
    const double boltzman = 0.001987191683; // kcal/mol

    // parameters 
    int N; // num of data points
    int n_replica = 0;
    int max_iter = 100000;

    double temp = 310; // temperature in K

    // io file params
    string outputfile = "weight.dat";
    string input_pattern;
    double tolerance = 1e-12;

    // arrays for data sets
    double *df;
    double *weight;
    double *old_weight;
    double *pmf;
    double *center;
    double *spring;
    int *n_samples;

    // wham params
    double rmsd = 10.0;
    double denom;
    double pot;
    double coor;
    double tmp;
    double sum;
    double delta;
    double first_val;
    int n = 0;

    //////
    cout << "Load params from stdin ..." << endl;
    string key;
    int value;
    double dvalue;
    while(!cin.eof()) {
        cin >> key;
        if (key == "TOLERANCE") { 
            cin >> tolerance;
        } else if (key == "N_REPLICA") {
            cin >> n_replica;
        } else if (key == "N_DATA") {
            cin >> N;
        } else if (key == "INPUT_PATTERN") {
            cin >> input_pattern;
        } else if (key == "WEIGHT_OUT") {
            cin >> outputfile;
        } else if (key == "TEMPERATURE") {
            cin >> temp;
        } else if (key == "MAX_ITER") {
            cin >> max_iter;
        } else if (key == "CENTERS") {
            if (! n_replica) {
                cout << "Define N_REPLICA first!" << endl;
                exit(1);
            }
            center = new double[n_replica];
            for (int i = 0; i < n_replica; ++i) {
                cin >> center[i];
            }
        } else if (key == "SPRINGS") {
            if (! n_replica) {
                cout << "Define N_REPLICA first!" << endl;
                exit(1);
            }
            spring = new double[n_replica];
            for (int i = 0; i < n_replica; ++i) {
                cin >> spring[i];
            }
        }
    }

    // read param from ARG file
    cout << "Load data ..." << endl;

    int step;
    double val;
    char filename[256];
    string buffer;
    df = new double[N*n_replica];
    n_samples = new int[N*n_replica];
    const char* INPUT_PATTERN = input_pattern.c_str();
    for (int i = 0; i < n_replica; ++i)
    {
        sprintf(filename, INPUT_PATTERN, i);

        cout << "Loading file " << filename << "..." << endl;
        ifstream ifs(filename);
        int n = 0;
        while(getline(ifs, buffer))
        {
            step = 0;
            val = 0.0;
            sscanf(buffer.data(), "%d %lf", &step, &val);
            df[N*i + n] = val;
            n++;
        }
        n_samples[i] = n;
        ifs.close();
    }

    // setup
    weight = new double[N*n_replica];
    old_weight = new double[N*n_replica];
    pmf = new double[n_replica];

    // calc
    const double beta = 1 / (temp * boltzman);
    int n_run = 0;
    while (rmsd > tolerance) {
        printf("WHAM run #%d\n", n_run);

        // backup old result
        memcpy(old_weight, weight, N*n_replica*(sizeof(double)));

        // update weight
        for (int i = 0; i < n_replica; ++i)
        {
            for (int t = 0; t < N; ++t)
            {
                denom = 0.0;
                for (int j = 0; j < n_replica; ++j)
                {
                    pot = pmf[j] - bias(df[N*i + t], center[j], spring[j]);
                    denom += n_samples[j] * exp(beta * pot);
                }
                weight[N*i + t] = 1.0 / denom;
            }
        }

        // update pmf
        for (int i = 0; i < n_replica; ++i)
        {
            tmp = 0.0;
            for (int t = 0; t < N; ++t)
            {
                for (int j = 0; j < n_replica; ++j)
                {
                    pot = bias(df[N*j + t], center[i], spring[i]);
                    tmp += weight[N*j + t] * exp(-1*beta*pot);
                }
            }
            pmf[i] = -1 * log(tmp) / beta;
        }

        // calc change of weights from last iteration
        sum = 0.0;
        n = 0;
        for (int i = 0; i < n_replica; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                delta = weight[N*i + j] - old_weight[N*i + j];
                sum += delta * delta;
                n++;
            }
        }
        rmsd = sqrt(sum / n);
        printf("Last change is %e\n", rmsd);

        // reset origin of pmf
        first_val = pmf[0];
        for (int i = 0; i < n_replica; ++i)
        {
            pmf[i] -= first_val;
        }

        n_run++;
        if (n_run > max_iter)
        {
            cout << "Not converged! Break." << endl;
            break;
        }
    }

    // wirte out the weight file
    ofstream w_out;
    w_out.open(outputfile);

    cout << "Writing weight file ..." << endl;
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < n_replica; ++i)
        {
            w_out << weight[N*i + j] << " ";
        }
        w_out << endl;
    }

    cout << "PMF per window:" << endl;
    for (int i = 0; i < n_replica; ++i)
    {
        //cout << "Window: " << i << ", PMF: " << pmf[i] << endl;
        printf("Window: %2d,PMF: %f\n", i, pmf[i]);
    }

    w_out.close();

    delete[] weight;
    delete[] old_weight;
    delete[] df;
    delete[] pmf;
    delete[] n_samples;
    delete[] center;
    delete[] spring;

    return 0;
}


double bias(double current, double center, double spring) {
    double bias = 0.0;
    double delta = 0.0;
    delta = current - center;
    bias = 0.5 * spring * delta * delta;
    return bias;
}


