#include <chrono>           /* Time */
#include <iostream>         /* cout, ... */
#include <fstream>          /* ifstream, ... */
#include <math.h>           /* ceil, ... */
#include <vector>           /* vector, ... */

using namespace std;


/* Variables */
string file;
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

/* Functions */
void readAtoms(string &file);
void headerSingle(string &file);
double factorCD(double q1, double q2, double d);
double factorLJ(double d);
double distance(double pos, double first, int n);
double LJ_per_mol(double pos, double dx, double ref);
double totalLJfactor();
double CD_per_mol(double pos, double q1, double dx, double ref);
double totalCDfactor();
void singleEmin();
void multipleEmin(int number);




int main(int argc, char const *argv[])
{
    auto start = chrono::high_resolution_clock::now();
    /************************************************************************/
    // file = "./dis/charge_distribution_36";            /* Anne's collagen */
    file = "./dis/charge_distribution_1054";       /* Native collagen */
    readAtoms(file);
    L = (N - 1) * distance_atoms;
    // for (int i = 0; i < (int) charge_ind.size(); i++) {
    //   cout << "\n " << charge_ind[i];
    // }
    /************************************************************************/
    cout << "\n\n:: running computation...";
    /* Parameterspace exploration */
    lat_gap = diameter_atom;
    // singleEmin();
    multipleEmin(200);

    // // double Z = 0.;
    // // for (int i = 0; i < (int) totalEnergies.size(); i++) {
    // //   Z += totalEnergies[i];
    // //   cout << "\ntotalE[" << i << "] = " << totalEnergies[i];
    // // }
    // // double avg_rad = 0.;
    // // double avg_off = 0.;
    // // double prob = 0.;
    // // int index = 0;
    // // for (int rad = 0; rad <= rad_max; rad++) {
    // //   rad_gap = rad * distance_atoms;
    // //   off_max = (int) ceil((L - rad_gap - 10.) / distance_atoms);
    // //   for (int off = 0; off <= off_max; off++) {
    // //     offset = rad_gap + 5. + off * distance_atoms;
    // //     prob = exp(-totalEnergies[index]) / Z;
    // //     avg_rad += prob * rad_gap;
    // //     avg_off += prob * off;
    // //     index++;
    // //   }
    // //   cout << "\navg_rad = " << avg_rad;
    // //   cout << "\navg_off = " << avg_off;
    // // }
    // //
    // // cout << "\n\nThe weighted average D-periodicity is:";
    // // cout << " " << L + avg_rad - avg_off;
    // // cout << "\n\nZ = " << Z;
    // // cout << "\n\nTotal energies.size() = " << totalEnergies.size();

    /************************************************************************/
    /* Testparameter */
    // lat_gap = 10 * 0.1 * diameter_atom;
    // // rad_gap = 0.25 * 12 * L + 0.0 * max_cutoff + 0;
    // // offset = 1.0 * 0 * L + 268.0;
    // // lat_gap = 1.12;
    // rad_gap = 1.02; //L + 2 * max_cutoff;
    // // offset = 0 * distance_atoms + 267.75;
    // box = L + rad_gap;
    //
    // cout << "\n\nlat_gap = " << lat_gap;
    // cout << "\nrad_gap = " << rad_gap;
    // cout << "\noffset = " << offset;
    // cout << "\nbox length = " << box;
    //
    // file = "./data/outputTest.dat";
    // FILE *test;
    // test = fopen(file.c_str(), "w");
    // int max_k = (int) ceil((L - rad_gap - 0.) / distance_atoms);
    // double factor = 100;
    // double factor2 = 0.05;
    // fprintf(test, "#atoms per molecule N = %i", N);
    // fprintf(test, "\n#distance_atoms = %.3f", distance_atoms);
    // fprintf(test, "\n#molecule_length = %.3f", L);
    // fprintf(test, "\n#CD divisor = %.3f", factor);
    // fprintf(test, "\n#LJ factor = %.3f", factor2);
    // fprintf(test, "\n#lat_gap = %.3f", lat_gap);
    // fprintf(test, "\n#rad_gap = %.3f", rad_gap);
    // fprintf(test, "\n\n#offset");
    // fprintf(test, "\tE_LJ");
    // fprintf(test, "\tE_CD");
    // fprintf(test, "\tE_total");
    // for (int k = 0; k <= max_k; k++) {
    //     offset = rad_gap + 0. + k * distance_atoms;
    //     fprintf(test, "\n%.3f", offset);
    //     double elj = totalLJfactor() * factor2;
    //     double ecd = totalCDfactor() / factor;
    //     fprintf(test, "\t%.3f", elj);
    //     fprintf(test, "\t%.3f", ecd);
    //     fprintf(test, "\t%.3f", elj + ecd);
    // }
    // fclose(test);

    // cout << "\n\nTotal LJ energy = " << 0.01 * totalLJfactor();
    // cout << "\nTotal CD energy = " << totalCDfactor() / 50.0;
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

void headerSingle(string &file)
{
  FILE *outf;
  outf = fopen(file.c_str(), "w");
  fprintf(outf, "#atoms per molecule N = %i", N);
  fprintf(outf, "\n#distance_atoms = %.3f", distance_atoms);
  fprintf(outf, "\n#molecule_length = %.3f", L);
  fprintf(outf, "\n\n#LJ_epsilon");
  fprintf(outf, "\tCD_epsilon");
  fprintf(outf, "\tlateral_gap_min");
  fprintf(outf, "\tradial_gap_min");
  fprintf(outf, "\toffset_min");
  fprintf(outf, "\tD-periodicity_min");
  fprintf(outf, "\tE_CD_min");
  fprintf(outf, "\tE_LJ_min");
  fprintf(outf, "\tE_total_min");
  fclose(outf);
}

// void headerMultiple(string &file)
// {
//   FILE *outf;
//   outf = fopen(file.c_str(), "w");
//   fprintf(outf, "#atoms per molecule N = %i", N);
//   fprintf(outf, "\n#distance_atoms = %.3f", distance_atoms);
//   fprintf(outf, "\n#molecule_length = %.3f", L);
//   fprintf(outf, "\n\n#LJ_epsilon");
//   fprintf(outf, "\tCD_epsilon");
//   fprintf(outf, "\tlateral_gap_min");
//   for (int i = 0; i < 10; i++) {
//       fprintf(outf, "\tradial_gap_min_%i", i);
//       fprintf(outf, "\toffset_min_%i", i);
//       fprintf(outf, "\tD-periodicity_min_%i", i);
//       fprintf(outf, "\tE_CD_min_%i", i);
//       fprintf(outf, "\tE_LJ_min_%i", i);
//       fprintf(outf, "\tE_total_min_%i", i);
//   }
//   fclose(outf);
// }

double factorCD(double q1, double q2, double d)
{
    double tmp = 0;
    double CDconst = 22.4 * 22.4;
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

void singleEmin()
{
  int rad_max = (int) ceil((L + 2 * max_cutoff) / distance_atoms);
  int off_max;
  double lj, lje;
  double cd, cde;
  vector<vector<double>> latmin, radmin, offmin, eCDmin, eLJmin, eTotmin;

  file = "./data/singleEmin.dat";

  /* Write header for outputfile */
  headerSingle(file);

  /* Assign needed data vectors */
  // latmin.assign(50, vector<double> (20, 0));
  radmin.assign(50, vector<double> (20, 0));
  offmin.assign(50, vector<double> (20, 0));
  eLJmin.assign(50, vector<double> (20, 1e6));
  eCDmin.assign(50, vector<double> (20, 1e6));
  eTotmin.assign(50, vector<double> (20, 1e6));

  for (int rad = 0; rad <= rad_max; rad++) {
    cout << "\n --> " << (100. * rad) / (1. * rad_max);
    cout << "\% done ...";
    rad_gap = rad * distance_atoms;
    box = L + rad_gap;
    off_max = (int) ceil((L - rad_gap - 0.) / distance_atoms);
    for (int off = 0; off <= off_max; off++) {
      offset = rad_gap + 0. + off * distance_atoms;
      /* Energy factors calculation for given parameters above */
      lj = totalLJfactor();
      cd = totalCDfactor();
      /* Both CD and LJ potential */
      for (int i = 0; i < 50; i++) {
          lje = (i + 1) * 0.01 * lj;
          for (int j = 0; j < 20; j++) {
              cde = cd / ((j + 1) * 10);

              if (lje + cde < eTotmin[i][j]) {
                  eLJmin[i][j] = lje;
                  eCDmin[i][j] = cde;
                  eTotmin[i][j] = lje + cde;
                  // latmin[i][j] = lat_gap;
                  radmin[i][j] = rad_gap;
                  offmin[i][j] = offset;
              }
          }
      }
    }
  }
  /************************************************************************/
  FILE *outf;
  outf = fopen(file.c_str(), "a");
  /* Both LJ and CD potential */
  for (int i = 0; i < 50; i++) {
      for (int j = 0; j < 20; j++) {
          fprintf(outf, "\n");
          fprintf(outf, "%.3f", (i + 1) * 0.01);
          fprintf(outf, "\t%.3f", (j + 1) * 10.);
          fprintf(outf, "\t%.3f", lat_gap);
          fprintf(outf, "\t%.3f", radmin[i][j]);
          fprintf(outf, "\t%.3f", offmin[i][j]);
          fprintf(outf, "\t%.3f", L + radmin[i][j] - offmin[i][j]);
          fprintf(outf, "\t%.3f", eCDmin[i][j]);
          fprintf(outf, "\t%.3f", eLJmin[i][j]);
          fprintf(outf, "\t%.3f", eTotmin[i][j]);
      }
      fprintf(outf, "\n");
  }
}

void multipleEmin(int number)
{
  int rad_max = (int) ceil((L + 2 * max_cutoff) / distance_atoms);
  int off_max;
  double lj, lje;
  double cd, cde;
  vector<vector<vector<double>>> latmin, radmin, offmin;
  vector<vector<vector<double>>> eCDmin, eLJmin, eTotmin;

  /* Write header for outputfile */
  for (int i = 1; i <= number; i++) {
    file = "./data/multipleEmins/Emin_" + to_string(i) + ".dat";
    headerSingle(file);
  }

  /* Assign needed data vectors */
  // latmin.assign(10, vector<vector<double>> (50, vector<double> (20, 0)));
  radmin.assign(number, vector<vector<double>> (50, vector<double> (20, 0)));
  offmin.assign(number, vector<vector<double>> (50, vector<double> (20, 0)));
  eLJmin.assign(number, vector<vector<double>> (50, vector<double> (20, 1e6)));
  eCDmin.assign(number, vector<vector<double>> (50, vector<double> (20, 1e6)));
  eTotmin.assign(number, vector<vector<double>> (50, vector<double> (20, 1e6)));

  for (int rad = 0; rad <= rad_max; rad++) {
    cout << "\n --> " << (100. * rad) / (1. * rad_max);
    cout << "\% done ...";
    rad_gap = rad * distance_atoms;
    box = L + rad_gap;
    off_max = (int) ceil((L - rad_gap - 0.) / distance_atoms);
    for (int off = 0; off <= off_max; off++) {
      offset = rad_gap + 0. + off * distance_atoms;
      /* Energy factors calculation for given parameters above */
      lj = totalLJfactor();
      cd = totalCDfactor();
      /* Both CD and LJ potential */
      for (int i = 0; i < 50; i++) {
          lje = (i + 1) * 0.01 * lj;
          for (int j = 0; j < 20; j++) {
              cde = cd / ((j + 1) * 10);

              for (int k = 0; k < number; k++) {
                if (lje + cde < eTotmin[k][i][j]) {
                    eLJmin[k][i][j] = lje;
                    eCDmin[k][i][j] = cde;
                    eTotmin[k][i][j] = lje + cde;
                    // latmin[k][i][j] = lat_gap;
                    radmin[k][i][j] = rad_gap;
                    offmin[k][i][j] = offset;
                    break;
                }
              }
          }
      }
    }
  }

  FILE *outf;
  for (int l = 1; l <= number; l++) {
    int k = l - 1;
    file = "./data/multipleEmins/Emin_" + to_string(l) + ".dat";
    outf = fopen(file.c_str(), "a");
    /* Both LJ and CD potential */
    for (int i = 0; i < 50; i++) {
        for (int j = 0; j < 20; j++) {
            fprintf(outf, "\n");
            fprintf(outf, "%.3f", (i + 1) * 0.01);
            fprintf(outf, "\t%.3f", (j + 1) * 10.);
            fprintf(outf, "\t%.3f", lat_gap);
            fprintf(outf, "\t%.3f", radmin[k][i][j]);
            fprintf(outf, "\t%.3f", offmin[k][i][j]);
            fprintf(outf, "\t%.3f", L + radmin[k][i][j] - offmin[k][i][j]);
            fprintf(outf, "\t%.3f", eCDmin[k][i][j]);
            fprintf(outf, "\t%.3f", eLJmin[k][i][j]);
            fprintf(outf, "\t%.3f", eTotmin[k][i][j]);
        }
        fprintf(outf, "\n");
    }
  }
}
