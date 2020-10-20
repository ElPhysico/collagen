#include <chrono>           /* Time */
#include <iostream>         /* cout, ... */
#include <fstream>          /* ifstream, ... */
#include <math.h>           /* ceil, ... */
#include <vector>           /* vector, ... */

using namespace std;


/* Variables */
int N = 0;                  /* number of atoms in one molecule */

double distance_atoms = 0.255;  /* distance between neighbouring atoms */
double diameter_atom = 1.12;    /* diameter of atom */
double L;                   /* length of collagen molecule */
double cd_cutoff = 2.0;     /* cutoff of cd potential */
double lj_cutoff = 2.0;     /* cutoff of lj potential */
double max_cutoff = max(cd_cutoff, lj_cutoff);
double lat_gap;             /* lateral gap = top <-> bottom distance */
double rad_gap;             /* distance end focus molecue to end of box */
double offset;              /* distance top molecule to start of box */
double box;                 /* length of periodic box */

vector<double> charge, charge_ind;  /* vector for the charges of each atom */
vector<vector<vector<double>>> latmin, radmin, offmin, eCDmin, eLJmin, eTotmin;
vector<vector <double>> CDLatmin, CDRadmin, CDOffmin, CDmin;

/* Functions */
void readAtoms(string &file);
void headerBoth(string &file);
void headerCD(string &file);
double factorCD(double q1, double q2, double d);
double factorLJ(double d);
double distance(double pos, double first, int n);
double LJ_per_mol(double pos, double dx, double ref);
double totalLJfactor();
double CD_per_mol(double pos, double q1, double dx, double ref);
double totalCDfactor();




int main(int argc, char const *argv[])
{
    auto start = chrono::high_resolution_clock::now();
    /************************************************************************/
    string file;
    // file = "./dis/charge_distribution_36";            /* Anne's collagen */
    file = "./dis/charge_distribution_1054";       /* Native collagen */
    readAtoms(file);
    L = (N - 1) * distance_atoms;
    // for (int i = 0; i < (int) charge_ind.size(); i++) {
    //   cout << "\n " << charge_ind[i];
    // }
    /************************************************************************/
    file = "./data/output.dat";    /* where the simulation results go */
    // file = "outputCD.dat";  /* only CD potential turned on */
    headerBoth(file);
    // headerCD(file);
    /************************************************************************/
    cout << "\n\n:: running computation...";
    /* Parameterspace exploration */
    // int lat_max = (int) ceil(max_cutoff / (0.1 * diameter_atom));
    int lat_max = 0;
    // int rad_max = (int) ceil(10.0 / 0.1);
    int rad_max = (int) ceil((L + 2 * max_cutoff) / distance_atoms);
    // int rad_max = (int) ceil((L + 2 * max_cutoff) / 0.1);
    int off_max;
    int maxRuns = (lat_max + 1) * (rad_max + 1);
    int counter = 0;
    double lj, lje;
    double cd, cde;

    /* Assign needed data vectors */
    // latmin.assign(10, vector<vector<double>> (50, vector<double> (20, 0)));
    radmin.assign(10, vector<vector<double>> (50, vector<double> (20, 0)));
    offmin.assign(10, vector<vector<double>> (50, vector<double> (20, 0)));
    eLJmin.assign(10, vector<vector<double>> (50, vector<double> (20, 1e6)));
    eCDmin.assign(10, vector<vector<double>> (50, vector<double> (20, 1e6)));
    eTotmin.assign(10, vector<vector<double>> (50, vector<double> (20, 1e6)));

    /* 10 to account for 10 minima */
    // CDLatmin.assign(10, vector<double> (20, 0));
    // CDRadmin.assign(10, vector<double> (20, 0));
    // CDOffmin.assign(10, vector<double> (20, 0));
    // CDmin.assign(10, vector<double> (20, 1e6));

    for (int lat = 0; lat <= lat_max; lat++) {
        // lat_gap = lat * 0.1 * diameter_atom;
        lat_gap = diameter_atom;

        for (int rad = 0; rad <= rad_max; rad++) {
            cout << "\n --> " << (100. * counter) / (1. * maxRuns);
            cout << "\% done ...";
            rad_gap = rad * distance_atoms;
            // rad_gap = rad * 0.1;
            box = L + rad_gap;
            // off_max = (int) ceil((box - 2.0) / distance_atoms);
            // off_max = (int) ceil(box / 0.1);
            off_max = (int) ceil((L - rad_gap) / distance_atoms);
            // cout << "\n\noffmax = " << off_max;
            for (int off = 0; off <= off_max; off++) {
                offset = rad_gap + off * distance_atoms;
                // offset = off * 0.1;
                /* separating the offset range into 1/10ths */
                // int tmp = (int) floor(10.0 * off / off_max);
                // tmp = min(tmp, 9);
                int tmp = 0;

                /* Energy calculation for given parameters above */
                lj = totalLJfactor();
                cd = totalCDfactor();

                /* Anne's original range:
                 * LJepsilon from 0.01 to 0.5 in 0.01 steps
                 * CDepsilon from 10 to 200 in 10 steps
                 */
                /* Both CD and LJ potential */
                for (int i = 0; i < 50; i++) {
                    lje = (i + 1) * 0.01 * lj;
                    for (int j = 0; j < 20; j++) {
                        cde = cd / ((j + 1) * 10);

                        if (lje + cde < eTotmin[tmp][i][j]) {
                            eLJmin[tmp][i][j] = lje;
                            eCDmin[tmp][i][j] = cde;
                            eTotmin[tmp][i][j] = lje + cde;
                            // latmin[tmp][i][j] = lat_gap;
                            radmin[tmp][i][j] = rad_gap;
                            offmin[tmp][i][j] = offset;
                        }
                    }
                }

                /* Only CD potential */
                // for (int i = 0; i < 20; i++) {
                //     cde = cd / ((i + 1) * 10);
                //
                //     if (cde < CDmin[tmp][i]) {
                //         CDmin[tmp][i] = cde;
                //         CDLatmin[tmp][i] = lat_gap;
                //         CDRadmin[tmp][i] = rad_gap;
                //         CDOffmin[tmp][i] = offset;
                //     }
                // }

            }
            counter++;
        }
    }
    /************************************************************************/
    /* Testparameter */
    // lat_gap = 10 * 0.1 * diameter_atom;
    // // rad_gap = 0.25 * 12 * L + 0.0 * max_cutoff + 0;
    // // offset = 1.0 * 0 * L + 268.0;
    // // lat_gap = 1.12;
    // rad_gap = 1.1; //L + 2 * max_cutoff;
    // offset = 0 * distance_atoms + 192.2;
    // box = L + rad_gap;
    //
    // cout << "\n\nlat_gap = " << lat_gap;
    // cout << "\nrad_gap = " << rad_gap;
    // cout << "\noffset = " << offset;
    // cout << "\nbox length = " << box;
    //
    // file = "outputTest.dat";
    // FILE *test;
    // test = fopen(file.c_str(), "w");
    // int max_k = (int) ceil(box / 0.1);
    // double factor = 10;
    // fprintf(test, "#atoms per molecule N = %i", N);
    // fprintf(test, "\n#distance_atoms = %.3f", distance_atoms);
    // fprintf(test, "\n#molecule_length = %.3f", L);
    // fprintf(test, "\n#lat_gap = %.3f", lat_gap);
    // fprintf(test, "\n#rad_gap = %.3f", rad_gap);
    // fprintf(test, "\n\n#offset");
    // fprintf(test, "\tE_LJ");
    // fprintf(test, "\tE_CD");
    // fprintf(test, "\tE_total");
    // for (int k = 0; k <= max_k; k++) {
    //     offset = k * 0.1;
    //     fprintf(test, "\n%.3f", offset);
    //     double elj = totalLJfactor() * 0.01;
    //     double ecd = totalCDfactor() / factor;
    //     fprintf(test, "\t%.3f", elj);
    //     fprintf(test, "\t%.3f", ecd);
    //     fprintf(test, "\t%.3f", elj + ecd);
    // }
    // fclose(test);

    // cout << "\n\nTotal LJ energy = " << 0.01 * totalLJfactor();
    // cout << "\nTotal CD energy = " << totalCDfactor() / 50.0;
    /************************************************************************/
    FILE *outf;
    outf = fopen(file.c_str(), "a");
    /* Both LJ and CD potential */
    for (int i = 0; i < 50; i++) {
        for (int j = 0; j < 20; j++) {
            fprintf(outf, "\n");
            fprintf(outf, "%.3f", (i + 1) * 0.01);
            fprintf(outf, "\t%.3f", (j + 1) * 10.);
            // fprintf(outf, "\t%.3f", latmin[0][i][j]);
            fprintf(outf, "\t%.3f", lat_gap);
            for (int k = 0; k < 1; k++) {
                fprintf(outf, "\t%.3f", radmin[k][i][j]);
                fprintf(outf, "\t%.3f", offmin[k][i][j]);
                fprintf(outf, "\t%.3f", L + radmin[k][i][j] - offmin[k][i][j]);
                fprintf(outf, "\t%.3f", eCDmin[k][i][j]);
                fprintf(outf, "\t%.3f", eLJmin[k][i][j]);
                fprintf(outf, "\t%.3f", eTotmin[k][i][j]);
            }
        }
        fprintf(outf, "\n");
    }
    /* Only CD potential */
    // for (int j = 0; j < 20; j++) {
    //     fprintf(outf, "\n");
    //     fprintf(outf, "%.3f", (j + 1) * 10.);
    //     for (int i = 0; i < 10; i++) {
    //         fprintf(outf, "\t\t\t%.3f", CDLatmin[i][j]);
    //         fprintf(outf, "\t\t\t%.3f", CDRadmin[i][j]);
    //         fprintf(outf, "\t\t%.3f", CDOffmin[i][j]);
    //         fprintf(outf, "\t\t\t%.3f", L + CDRadmin[i][j] - CDOffmin[i][j]);
    //         fprintf(outf, "\t\t\t%.3f", CDmin[i][j]);
    //     }
    // }
    // fclose(outf);
    /************************************************************************/
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start);
    cout << "\n\n\n:: ... finished in " << duration.count() << " s.\n\n\n";
    /************************************************************************/

    return 0;
}




/* Functions */
void readAtoms(string &file)
{
    cout << "\n\n:: reading atom configuration...";
    string line;
    ifstream myfile;
    myfile.open(file.c_str(), ios::in);
    while (myfile.peek() != EOF) {
        getline(myfile, line);
        N++;
        charge.push_back(atof(line.c_str()));
        if (charge[N - 1] != 0) {
          charge_ind.push_back(N - 1);
        }
        // cout << "\natom " << N << " has charge " << line;
    }
    myfile.close();
    cout << "\n -> " << N << " atoms read.";
    L = (N - 1) * distance_atoms;
    cout << " Molecule length L = " << L << ".";
}

void headerBoth(string &file)
{
    FILE *outf;
    outf = fopen(file.c_str(), "w");
    fprintf(outf, "#atoms per molecule N = %i", N);
    fprintf(outf, "\n#distance_atoms = %.3f", distance_atoms);
    fprintf(outf, "\n#molecule_length = %.3f", L);
    fprintf(outf, "\n\n#LJ_epsilon");
    fprintf(outf, "\tCD_epsilon");
    fprintf(outf, "\tlateral_gap_min");
    for (int i = 0; i < 10; i++) {
        fprintf(outf, "\tradial_gap_min_%i", i);
        fprintf(outf, "\toffset_min_%i", i);
        fprintf(outf, "\tD-periodicity_min_%i", i);
        fprintf(outf, "\tE_CD_min_%i", i);
        fprintf(outf, "\tE_LJ_min_%i", i);
        fprintf(outf, "\tE_total_min_%i", i);
    }
    fclose(outf);
}

void headerCD(string &file)
{
    FILE *outf;
    outf = fopen(file.c_str(), "w");
    fprintf(outf, "#atoms per molecule N = %i", N);
    fprintf(outf, "\n#distance_atoms = %.3f", distance_atoms);
    fprintf(outf, "\n#molecule_length = %.3f", L);
    fprintf(outf, "\n\n#CD_epsilon");
    for (int i = 0; i < 10; i++) {
        fprintf(outf, "\tlateral_gap_min_%i", i);
        fprintf(outf, "\tradial_gap_min_%i", i);
        fprintf(outf, "\toffset_min_%i", i);
        fprintf(outf, "\tD-periodicity_min_%i", i);
        fprintf(outf, "\tE_CD_min_%i", i);
    }
    fclose(outf);
}

double factorCD(double q1, double q2, double d)
{
    double tmp = 0;
    double CDconst = 501.7619;
    double CDdamping = 1;

    tmp = (CDconst * q1 * q2) / d;
    tmp *= exp(-1.0 * CDdamping * d);

    return tmp;
}

double factorLJ(double d)
{
    return 4 * (pow(1 / d, 12) - pow(1 / d, 6));
}

double distance(double pos, double first, int n)
{
    double d;
    d = first + n * distance_atoms - pos;
    d *= d;
    d += lat_gap * lat_gap;
    d = sqrt(d);
    return d;
}

double LJ_per_mol(double pos, double dx, double ref)
{
    double d, sum;
    int left, right;
    bool first = true;
    sum = 0;

    left = ceil(((pos - dx) - ref) / distance_atoms);
    right = floor(((pos + dx) - ref) / distance_atoms);

    if (left < N && right >= 0) {
        /* left end */
        if (left < 0 && right < N) {
            /* left is replaced by 0 here */
            if (first) {
                sum = 0;
                for (int i = 0; i <= right; i++) {
                    d = distance(pos, ref, i);
                    sum += factorLJ(d);
                }
                first = false;
            } else {
                d = distance(pos, ref, 0);
                sum += factorLJ(d);
            }
            if (left == -1) first = true;
        }

        /* bulk */
        if (left >= 0 && right < N) {
            if (first) {
                sum = 0;
                for (int i = 0; i <= right - left; i++) {
                    d = distance(pos, ref, left + i);
                    sum += factorLJ(d);
                }
                first = false;
            }
            if (right == N - 1) first = true;
        }

        /* right end */
        if (left >= 0 && right >= N) {
            /* right is replaced by N - 1 */
            if (first) {
                sum = 0;
                for (int i = 0; i <= N - 1 - left; i++) {
                    d = distance(pos, ref, left + i);
                    sum += factorLJ(d);
                }
                first = false;
            } else {
                d = distance(pos, ref, N);
                sum -= factorLJ(d);
            }
            if (left == N - 1) first = true;
        }
    }
    return sum;
}

double totalLJfactor()
{
    double pos, sum;
    int left, right;
    double dx = sqrt(lj_cutoff * lj_cutoff - lat_gap * lat_gap);
    sum = 0;
    for (int atom = 1; atom <= N; atom++) {
        pos = (atom - 1) * distance_atoms;

        /* Atom atom feels from top left molecule: */
        sum += LJ_per_mol(pos, dx, offset - box);

        /* Atom atom feels from top molecule: */
        sum += LJ_per_mol(pos, dx, offset);

        /* Atom atom feels from top right molecule: */
        sum += LJ_per_mol(pos, dx, offset + box);

        /* Atom atom feels from left molecule: */
        left = ceil(((pos - lj_cutoff) + box) / distance_atoms);
        right = floor(((pos + lj_cutoff) + box) / distance_atoms);
        if (left < N && right >= 0) {
            left = max(left, 0);
            right = min(right, N - 1);
            for (int i = 0; i <= right - left; i++) {
                double d = abs(pos + box - (left + i) * distance_atoms);
                sum += factorLJ(d);
            }
        }

        /* Atom atom feels from right molecule: */
        left = ceil(((pos - lj_cutoff) - box) / distance_atoms);
        right = floor(((pos + lj_cutoff) - box) / distance_atoms);
        if (left < N && right >= 0) {
            left = max(left, 0);
            right = min(right, N - 1);
            for (int i = 0; i <= right - left; i++) {
                double d = abs(pos - box - (left + i) * distance_atoms);
                sum += factorLJ(d);
            }
        }
    }
    return sum;
}

double CD_per_mol(double pos, double q1, double dx, double ref)
{
    double d, sum, q2;
    int left, right;
    sum = 0;

    left = ceil(((pos - dx) - ref) / distance_atoms);
    right = floor(((pos + dx) - ref) / distance_atoms);

    if (left < N && right >= 0) {
        left = max(left, 0);
        right = min(right, N - 1);
        for (int i = 0; i <= right - left; i++) {
            q2 = charge[left + i];
            if (q2 != 0) {
                d = distance(pos, ref, left + i);
                // cout << "\nhit with " << factorCD(q1, q2, d) / 150;
                // cout << " to change sum = " << sum / 150;
                sum += factorCD(q1, q2, d);
                // cout << " to sum = " << sum / 150;
            }
        }
    }
    return sum;
}

double totalCDfactor()
{
    double pos, sum, q1, q2, d;
    int left, right;
    double dx = sqrt(cd_cutoff * cd_cutoff - lat_gap * lat_gap);
    sum = 0;
    for (int atom = 0; atom < (int) charge_ind.size(); atom++) {
        q1 = charge[charge_ind[atom]];
        pos = charge_ind[atom] * distance_atoms;

        /* Atom atom feels from top left molecule: */
        sum += CD_per_mol(pos, q1, dx, offset - box);

        /* Atom atom feels from top molecule: */
        sum += CD_per_mol(pos, q1, dx, offset);

        /* Atom atom feels from top right molecule: */
        sum += CD_per_mol(pos, q1, dx, offset + box);

        /* Atom atom feels from left molecule: */
        left = ceil(((pos - cd_cutoff) + box) / distance_atoms);
        right = floor(((pos + cd_cutoff) + box) / distance_atoms);
        if (left < N && right >= 0) {
            left = max(left, 0);
            right = min(right, N - 1);
            for (int i = 0; i <= right - left; i++) {
                q2 = charge[left + i];
                if (q2 != 0) {
                    d = abs(pos + box - (left + i) * distance_atoms);
                    sum += factorCD(q1, q2, d);
                }
            }
        }

        /* Atom atom feels from right molecule: */
        left = ceil(((pos - cd_cutoff) - box) / distance_atoms);
        right = floor(((pos + cd_cutoff) - box) / distance_atoms);
        if (left < N && right >= 0) {
            left = max(left, 0);
            right = min(right, N - 1);
            for (int i = 0; i <= right - left; i++) {
                q2 = charge[left + i];
                if (q2 != 0) {
                    d = abs(pos - box - (left + i) * distance_atoms);
                    sum += factorCD(q1, q2, d);
                }
            }
        }

        // cout << "\nAfter atom " << atom << " the CD energy is " << sum / 150.0;
        // cout << "\t";
    }
    return sum;
}
