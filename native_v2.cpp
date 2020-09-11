#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <chrono>

using namespace std;


/* Variables */
int N;  /* atoms per molecue */
double L;   /* length of molecule */

double A_init = 5;
double A_end = 50.0;
double CD_epsilon = 1.;
double CD_const = 1.;
double CD_damping = 1.;
double CD_cutoff = 2.;
int CD_count = 21;
double CD_inc = (A_end - A_init) / (1. * CD_count);

double LJ_init = 0.01;
double LJ_end = 0.5;
double LJ_epsilon = 1.;
double LJ_cutoff = 2.;
int LJ_count = 51;
double LJ_inc = (LJ_end - LJ_init) / (1. * LJ_count);

double cutoff_max = max(CD_cutoff, LJ_cutoff);

double distance_atoms = 0.255;
double diameter_atom = 1.12;
double lat_gap, rad_gap, offset, box_length;
double rad_inc = 0.1;
double offset_inc = 1.0;

// double LJ_energy_min, CD_energy_min, energy_min;

vector<double> charge;
vector<vector<double>> eTotal, eCD, eLJ, eTotalmin, eCDmin, eLJmin;
vector<vector<int>> latmin, radmin, offmin;



/* Functions */
void collectAtoms();
void header(string &file);
double e_CD(double charge1, double charge2, double d);
double e_LJ(double d);
vector<vector<double>> calc_LJ_top(double x);
vector<vector<double>> calc_CD_top(double x);
vector<vector<double>> calc_LJ_side();
vector<vector<double>> calc_CD_side();
vector<vector<double>> add_E(vector<vector<double>> v1, vector<vector<double>> v2);


int main(int argc, char const *argv[])
{
    cout << ":: Starting...\n";
    auto start = chrono::high_resolution_clock::now();

    /* Initialize tracking vectors */
    eTotal.assign(CD_count, vector<double> (LJ_count, 0));
    eCD.assign(CD_count, vector<double> (LJ_count, 0));
    eLJ.assign(CD_count, vector<double> (LJ_count, 0));
    eTotalmin.assign(CD_count, vector<double> (LJ_count, 1e6));
    eCDmin.assign(CD_count, vector<double> (LJ_count, 1e6));
    eLJmin.assign(CD_count, vector<double> (LJ_count, 1e6));
    latmin.assign(CD_count, vector<int> (LJ_count, 0));
    radmin.assign(CD_count, vector<int> (LJ_count, 0));
    offmin.assign(CD_count, vector<int> (LJ_count, 0));

    /* Collecting atoms, charge distribution from external file */
    collectAtoms();
    L = (N - 1) * distance_atoms;

    /* Writing header */
    string file = "./energy_min_v2.dat";
    header(file);

    /* Calculation of energy landscape */
    cout << "\n\n:: calculation of energy landscape...";
    int counter = 0;
    int lat_gap_max = (int) ceil(cutoff_max / (0.1 * diameter_atom));
    int rad_gap_max = (int) ceil(10 / rad_inc);
    // int rad_gap_max = (int) ceil((L + 2 * cutoff_max) / rad_inc);
    int max_runs = (lat_gap_max + 1) * (rad_gap_max + 1);
    for (int lat = 0; lat <= lat_gap_max; lat++) {
        lat_gap = lat * 0.1 * diameter_atom;
        for (int rad = 0; rad <= rad_gap_max; rad++) {
            cout << "\n -> " << (100. * counter) / (1. * max_runs) << "\% done...";
            rad_gap = rad * rad_inc;
            box_length = L + rad_gap;
            int offset_max = (int) ceil(box_length / offset_inc);
            // cout << "\nlatgapmax = " << lat_gap_max << " " << lat_gap_max * 0.1 * diameter_atom;
            // cout << "\nradgapmax = " << rad_gap_max << " " << rad_gap_max * 0.1;
            // cout << "\noffsetmax = " << offset_max << " " << offset_max * 0.1;
            for (int off = 0; off <= offset_max; off++) {
                offset = off * offset_inc;
                /* reset tracking vectors */
                // lat_gap = diameter_atom;
                // rad_gap = diameter_atom;
                // offset = 0;
                eCD.assign(CD_count, vector<double> (LJ_count, 0));
                eLJ.assign(CD_count, vector<double> (LJ_count, 0));

                /* top middle molecule */
                eLJ = add_E(eLJ, calc_LJ_top(offset));
                eCD = add_E(eCD, calc_CD_top(offset));
                /* top left molecule */
                eLJ = add_E(eLJ, calc_LJ_top(offset - box_length));
                eCD = add_E(eCD, calc_CD_top(offset - box_length));
                /* top right molecule */
                eLJ = add_E(eLJ, calc_LJ_top(offset + box_length));
                eCD = add_E(eCD, calc_CD_top(offset + box_length));
                /* side molecules */
                eLJ = add_E(eLJ, calc_LJ_side());
                eCD = add_E(eCD, calc_CD_side());

                for (int cd = 0; cd < CD_count; cd++) {
                    for (int lj = 0; lj < LJ_count; lj++) {
                        double cLJ = eLJ[cd][lj];
                        double cCD = eCD[cd][lj];
                        // if (cLJ < 0) {
                        //     cout << "\ncLJ = " << cLJ;
                        //     cout << "\ncCD = " << cCD;
                        // }
                        if (cLJ + cCD < eTotalmin[cd][lj]) {
                            eCDmin[cd][lj] = cCD;
                            eLJmin[cd][lj] = cLJ;
                            eTotalmin[cd][lj] = cCD + cLJ;
                            latmin[cd][lj] = lat;
                            radmin[cd][lj] = rad;
                            offmin[cd][lj] = off;
                        }
                    }
                }
            }
            counter++;
        }
    }

    FILE *outf;
    outf = fopen(file.c_str(), "a");
    for (int lj = 0; lj < LJ_count; lj++) {
        for (int cd = 0; cd < CD_count; cd++) {
            fprintf(outf, "\n");
            fprintf(outf, "%.3f", LJ_epsilon = LJ_init + lj * LJ_inc); // (lj + 1) * 0.01;);
            fprintf(outf, "\t\t%.3f", A_init + cd * CD_inc); //(cd + 1) * 10;);
            fprintf(outf, "\t\t\t%.3f", latmin[cd][lj] * 0.1 * diameter_atom);
            fprintf(outf, "\t\t\t%.3f", radmin[cd][lj] * rad_inc);
            fprintf(outf, "\t\t%.3f", offmin[cd][lj] * offset_inc);
            fprintf(outf, "\t\t\t%.3f", L + radmin[cd][lj] * rad_inc - offmin[cd][lj] * offset_inc);
            fprintf(outf, "\t\t\t%.3f", eCDmin[cd][lj]);
            fprintf(outf, "\t\t%.3f", eLJmin[cd][lj]);
            fprintf(outf, "\t\t%.3f", eTotalmin[cd][lj]);
        }
        fprintf(outf, "\n");
    }
    fclose(outf);



    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start);
    cout << "\n\n\n:: ... finished in " << duration.count() << " s.\n\n\n";
    return 0;
}


/* Functions */
void collectAtoms()
{
    cout << "\n\n:: collecting atoms, charge distribution...";
    charge.clear(); /* removes all elements of vector -> making its size 0 */
    N = 0;
    string input;
    ifstream myfile;
    myfile.open("charge_distribution", ios::in);
    while (myfile.peek() != EOF) {
        getline(myfile, input);
        N++;
        charge.push_back(atof(input.c_str()));
    }
    myfile.close();
    cout << "\n -> " << N << " atoms collected.";
}
void header(string &file)
{
    FILE *outf;
    outf = fopen(file.c_str(), "w");
    fprintf(outf, "#atoms per molecule N = %i", N);
    fprintf(outf, "\n#distance_atoms = %.3f", distance_atoms);
    fprintf(outf, "\n#molecule_length = %.3f", L);
    fprintf(outf, "\n\n#LJ_epsilon");
    fprintf(outf, "\tCD_A");
    fprintf(outf, "\tlateral_gap_min");
    fprintf(outf, "\tradial_gap_min");
    fprintf(outf, "\toffset_min");
    fprintf(outf, "\tD-periodicity_min");
    fprintf(outf, "\tE_CD_min");
    fprintf(outf, "\tE_LJ_min");
    fprintf(outf, "\tE_total_min");
    fclose(outf);
}
double e_CD(double charge1, double charge2, double d)
{
    // double tmp = (CD_const * charge1 * charge2) / (CD_epsilon * d);
    double tmp = (CD_const * CD_epsilon * charge1 * charge2) / d;
    return tmp * exp(-1. * CD_damping * d);
}
double e_LJ(double d)
{
    return 4. * LJ_epsilon * (pow(1. / d, 12.) - pow(1. / d, 6.));
}
vector<vector<double>> calc_LJ_top(double x)
{
    double atom_offset; /* offset from the imaginary "distance_atom" grid */
    int first_atom; /* reference atom with respect to grid */
    int left_atoms, right_atoms; /* how many interact to left and right */
    int current_atom, start_bulk, end_bulk;
    vector<vector<double>> current_atom_energy, LJ_energy;
    current_atom_energy.assign(CD_count, vector<double> (LJ_count, 0));
    LJ_energy.assign(CD_count, vector<double> (LJ_count, 0));
    double d, tmp;
    double cutoff_x = sqrt(LJ_cutoff * LJ_cutoff - lat_gap * lat_gap);

    first_atom = (int) floor(x / distance_atoms);
    atom_offset = x - first_atom * distance_atoms;
    tmp = cutoff_x + atom_offset - distance_atoms;
    left_atoms = (int) ceil(tmp / distance_atoms);
    tmp = cutoff_x - atom_offset;
    right_atoms = (int) ceil(tmp / distance_atoms);
    /* left border */
    for (int i = 0; i < left_atoms + right_atoms; i++) {
        current_atom = first_atom - right_atoms + 1 + i;
        if (i < right_atoms) {
            d = (right_atoms - 1 - i) * distance_atoms + atom_offset;
        } else {
            d = (i - right_atoms + 1) * distance_atoms - atom_offset;
        }
        d *= d;
        d += lat_gap * lat_gap;
        d = sqrt(d);
        for (int cd = 0; cd < CD_count; cd++) {
            CD_epsilon = A_init + cd * CD_inc; //(cd + 1) * 10;
            for (int lj = 0; lj < LJ_count; lj++) {
                LJ_epsilon = LJ_init + lj * LJ_inc; // (lj + 1) * 0.01;
                current_atom_energy[cd][lj] += e_LJ(d);
                if (current_atom >= 0 && current_atom <= N - 1) {
                    LJ_energy[cd][lj] += current_atom_energy[cd][lj];
                }
            }
        }
    }
    start_bulk = current_atom;
    /* right border */
    current_atom_energy.assign(CD_count, vector<double> (LJ_count, 0));
    for (int i = 0; i < left_atoms + right_atoms; i++) {
        current_atom = first_atom + N - 1 + left_atoms - i;
        if (i < left_atoms) {
            d = (left_atoms - i) * distance_atoms - atom_offset;
        } else {
            d = (i - left_atoms) * distance_atoms + atom_offset;
        }
        d *= d;
        d += lat_gap * lat_gap;
        d = sqrt(d);
        for (int cd = 0; cd < CD_count; cd++) {
            CD_epsilon = A_init + cd * CD_inc; //(cd + 1) * 10;
            for (int lj = 0; lj < LJ_count; lj++) {
                LJ_epsilon = LJ_init + lj * LJ_inc; // (lj + 1) * 0.01;
                current_atom_energy[cd][lj] += e_LJ(d);
                if (current_atom >= 0 && current_atom <= N - 1) {
                    LJ_energy[cd][lj] += current_atom_energy[cd][lj];
                }
            }
        }
    }
    end_bulk = current_atom;
    /* bulk atoms -> max energy contribution */
    if (start_bulk < N && end_bulk >= 0) {
        if(start_bulk < 0) start_bulk = 0;
        if (end_bulk > N - 1) end_bulk = N - 1;
        for (int cd = 0; cd < CD_count; cd++) {
            for (int lj = 0; lj < LJ_count; lj++) {
                int help = end_bulk - start_bulk - 1;
                LJ_energy[cd][lj] += help * current_atom_energy[cd][lj];
            }
        }
    }

    return LJ_energy;
}
vector<vector<double>> calc_CD_top(double x)
{
    double atom_offset; /* offset from the imaginary "distance_atom" grid */
    int first_atom; /* reference atom with respect to grid */
    int left_atoms, right_atoms; /* how many interact to left and right */
    vector<vector<double>> current_atom_energy, CD_energy;
    current_atom_energy.assign(CD_count, vector<double> (LJ_count, 0));
    CD_energy.assign(CD_count, vector<double> (LJ_count, 0));
    vector<double> d_left, d_right;
    double tmp;
    double cutoff_x = sqrt(CD_cutoff * CD_cutoff - lat_gap * lat_gap);

    first_atom = (int) floor(x / distance_atoms);
    atom_offset = x - first_atom * distance_atoms;
    tmp = cutoff_x + atom_offset - distance_atoms;
    left_atoms = (int) ceil(tmp / distance_atoms);
    tmp = cutoff_x - atom_offset;
    right_atoms = (int) ceil(tmp / distance_atoms);

    /* calculate distances of atoms to the left and right of focus atom */
    for (int i = 1; i <= left_atoms; i++) {
        double d;
        d = i * distance_atoms - atom_offset;
        d *= d;
        d += lat_gap * lat_gap;
        d = sqrt(d);
        d_left.push_back(d);
    }
    for (int i = 0; i < right_atoms; i++) {
        double d;
        d = i * distance_atoms + atom_offset;
        d *= d;
        d += lat_gap * lat_gap;
        d = sqrt(d);
        d_right.push_back(d);
    }
    /* cycle through all focus molecule atoms */
    for (int i = 0; i < N; i++) {
        current_atom_energy.assign(CD_count, vector<double> (LJ_count, 0));
        if (charge[i] != 0) {
            /* relevant atoms to the left */
            for (int j = 1; j <= left_atoms; j++) {
                int atom = i - j;
                atom -= first_atom;
                if (atom < N && atom >= 0 && charge[atom] != 0) {
                    for (int cd = 0; cd < CD_count; cd++) {
                        CD_epsilon = A_init + cd * CD_inc; //(cd + 1) * 10;
                        for (int lj = 0; lj < LJ_count; lj++) {
                            LJ_epsilon = LJ_init + lj * LJ_inc; // (lj + 1) * 0.01;
                            current_atom_energy[cd][lj] += e_CD(charge[i], charge[atom], d_left[j - 1]);
                        }
                    }
                }
            }
            /* relevant atoms to the right */
            for (int j = 0; j < right_atoms; j++) {
                int atom = i + j;
                atom -= first_atom;
                if (atom < N && atom >= 0 && charge[atom] != 0) {
                    for (int cd = 0; cd < CD_count; cd++) {
                        CD_epsilon = A_init + cd * CD_inc; //(cd + 1) * 10;
                        for (int lj = 0; lj < LJ_count; lj++) {
                            LJ_epsilon = LJ_init + lj * LJ_inc; // (lj + 1) * 0.01;
                            current_atom_energy[cd][lj] += e_CD(charge[i], charge[atom], d_right[j]);
                        }
                    }
                }
            }
            for (int cd = 0; cd < CD_count; cd++) {
                for (int lj = 0; lj < LJ_count; lj++) {
                    CD_energy[cd][lj] += current_atom_energy[cd][lj];
                }
            }
        }
    }

    return CD_energy;
}
vector<vector<double>> calc_LJ_side()
{
    vector<vector<double>> current_atom_energy, LJ_energy;
    current_atom_energy.assign(CD_count, vector<double> (LJ_count, 0));
    LJ_energy.assign(CD_count, vector<double> (LJ_count, 0));
    if (rad_gap <= LJ_cutoff) {
        int max_atoms = ceil(LJ_cutoff / distance_atoms);
        for (int i = 0; i <= max_atoms; i++) {
            int active_atom = max_atoms - i;
            double d = active_atom * distance_atoms + rad_gap;
            if (d <= LJ_cutoff) {
                for (int cd = 0; cd < CD_count; cd++) {
                    CD_epsilon = A_init + cd * CD_inc; //(cd + 1) * 10;
                    for (int lj = 0; lj < LJ_count; lj++) {
                        LJ_epsilon = LJ_init + lj * LJ_inc; // (lj + 1) * 0.01;
                        current_atom_energy[cd][lj] += e_LJ(d);
                    }
                }
            }
            for (int cd = 0; cd < CD_count; cd++) {
                for (int lj = 0; lj < LJ_count; lj++) {
                    LJ_energy[cd][lj] += current_atom_energy[cd][lj];
                }
            }
        }
    }
    for (int cd = 0; cd < CD_count; cd++) {
        for (int lj = 0; lj < LJ_count; lj++) {
            LJ_energy[cd][lj] *= 2;
        }
    }

    return LJ_energy;
}
vector<vector<double>> calc_CD_side()
{
    vector<vector<double>> current_atom_energy, CD_energy;
    current_atom_energy.assign(CD_count, vector<double> (LJ_count, 0));
    CD_energy.assign(CD_count, vector<double> (LJ_count, 0));
    if (rad_gap <= CD_cutoff) {
        int max_atoms = ceil(CD_cutoff / distance_atoms);
        for (int i = 0; i < max_atoms; i++) {
            current_atom_energy.assign(CD_count, vector<double> (LJ_count, 0));
            for (int j = 0; j < max_atoms - i; j++) {
                double d = (i + j) * distance_atoms + rad_gap;
                if (d <= CD_cutoff) {
                    int c1 = charge[i];
                    int c2 = charge[N - 1 - j];
                    for (int cd = 0; cd < CD_count; cd++) {
                        CD_epsilon = A_init + cd * CD_inc; //(cd + 1) * 10;
                        for (int lj = 0; lj < LJ_count; lj++) {
                            LJ_epsilon = LJ_init + lj * LJ_inc; // (lj + 1) * 0.01;
                            current_atom_energy[cd][lj] += e_CD(c1, c2, d);
                        }
                    }
                }
            }
            for (int cd = 0; cd < CD_count; cd++) {
                for (int lj = 0; lj < LJ_count; lj++) {
                    CD_energy[cd][lj] += current_atom_energy[cd][lj];
                }
            }
        }
    }
    for (int cd = 0; cd < CD_count; cd++) {
        for (int lj = 0; lj < LJ_count; lj++) {
            CD_energy[cd][lj] *= 2;
        }
    }

    return CD_energy;
}
vector<vector<double>> add_E(vector<vector<double>> v1, vector<vector<double>> v2)
{
    for (int cd = 0; cd < CD_count; cd++) {
        for (int lj = 0; lj < LJ_count; lj++) {
            v1[cd][lj] += v2[cd][lj];
        }
    }
    return v1;
}
