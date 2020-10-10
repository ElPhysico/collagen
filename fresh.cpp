#include <chrono>           /* Time */
#include <iostream>         /* cout, ... */
#include <fstream>          /* ifstream, ... */
#include <math.h>           /* ceil, ... */

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
double dx;                  /* max delta x to be in cutoff range */

/* Functions */
void readAtoms(string &file);
void header(string &file);


int main(int argc, char const *argv[])
{
    auto start = chrono::high_resolution_clock::now();
    /************************************************************************/
    string file;
    file = "charge_distribution_36";            /* Anne's collagen */
    // file = "charge_distribution_1054";       /* Native collagen */
    readAtoms(file);
    L = (N - 1) * distance_atoms;
    /************************************************************************/
    file = "output.dat";    /* where the simulation results go */
    header(file);
    /************************************************************************/
    cout << "\n\n:: running computation...";
    // int lat_max = (int) ceil(max_cutoff / (0.1 * diameter_atom));
    // for (int lat = 0; lat < lat_max; lat++) {
    //     lat_gap = lat * 0.1 * diameter_atom;
    //     for (int j = 0; j < 10; j++) {
    //         for (int k = 0; k < 10; k++) {
    //
    //         }
    //     }
    // }
    lat_gap = 10 * 0.1 * diameter_atom;
    dx = sqrt(max_cutoff * max_cutoff - lat_gap * lat_gap);
    rad_gap = 0.25 * 0 * L + 0.5 * max_cutoff;
    offset = 0.25 * 1 * L;
    box = L + rad_gap;
    cout << "\n\nlat_gap = " << lat_gap;
    cout << "\nrad_gap = " << rad_gap;
    cout << "\noffset = " << offset;
    cout << "\nbox length = " << box;

    cout << "\n\ndx = " << dx;
    cout << "\ntop first atom at: " << offset;
    cout << "\ntop left first atom at: " << offset - box;
    cout << "\ntop right first atom at: " << offset + box;

    // int focusAtom = 3;
    // double pos = (focusAtom - 1) * distance_atoms;
    // cout << "\n\nTest Atom Postion: " << pos;
    // cout << "\nlinks top von " << ceil(((pos - dx) - (offset - box)) / distance_atoms);
    // cout << " bis " << floor(((pos + dx) - (offset - box)) / distance_atoms);
    // cout << "\ntop von " << ceil(((pos - dx) - offset) / distance_atoms);
    // cout << " bis " << floor(((pos + dx) - offset) / distance_atoms);
    // cout << "\nrechts top von " << ((pos - dx) - (offset + box)) / distance_atoms;
    // cout << " bis " << ((pos + dx) - (offset + box)) / distance_atoms;
    // cout << "\nlinks von " << ((pos - max_cutoff) - (0 - box)) / distance_atoms;
    // cout << " bis " << ((pos + max_cutoff) - (0 - box)) / distance_atoms;
    // cout << "\nrechts von " << ((pos - max_cutoff) - (0 + box)) / distance_atoms;
    // cout << " bis " << ((pos + max_cutoff) - (0 + box)) / distance_atoms;

    double pos;
    int left, right;
    for (int atom = 1; atom <= N; atom++) {
        pos = (atom - 1) * distance_atoms;
        /* top left molecule */
        left = ceil(((pos - dx) - (offset - box)) / distance_atoms);
        right = floor(((pos + dx) - (offset - box)) / distance_atoms);
        cout << "\nAtom " << atom << " feels atoms " << left << " to ";
        cout << right << " from the top left molecule,";
        if (left < N && right >= 0) {
            left = max(left, 0);
            right = min(right, N - 1);
            cout << " --> " << left << " to " << right;
        }
        /* top molecule */
        left = ceil(((pos - dx) - (offset)) / distance_atoms);
        right = floor(((pos + dx) - (offset)) / distance_atoms);
        cout << " atoms " << left << " to ";
        cout << right << " from the top molecule,";
        if (left < N && right >= 0) {
            left = max(left, 0);
            right = min(right, N - 1);
            cout << " --> " << left << " to " << right;
        }
        /* top right molecule */
        left = ceil(((pos - dx) - (offset + box)) / distance_atoms);
        right = floor(((pos + dx) - (offset + box)) / distance_atoms);
        cout << " atoms " << left << " to ";
        cout << right << " from the top right molecule,";
        if (left < N && right >= 0) {
            left = max(left, 0);
            right = min(right, N - 1);
            cout << " --> " << left << " to " << right;
        }
        /* left molecule */
        left = ceil(((pos - max_cutoff) - (-box)) / distance_atoms);
        right = floor(((pos + max_cutoff) - (-box)) / distance_atoms);
        cout << " atoms " << left << " to ";
        cout << right << " from the left molecule,";
        if (left < N && right >= 0) {
            left = max(left, 0);
            right = min(right, N - 1);
            cout << " --> " << left << " to " << right;
        }
        /* right molecule */
        left = ceil(((pos - max_cutoff) - (box)) / distance_atoms);
        right = floor(((pos + max_cutoff) - (box)) / distance_atoms);
        cout << " atoms " << left << " to ";
        cout << right << " from the right molecule.";
        if (left < N && right >= 0) {
            left = max(left, 0);
            right = min(right, N - 1);
            cout << " --> " << left << " to " << right;
        }
    }

    /************************************************************************/
    /************************************************************************/
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start);
    cout << "\n\n\n:: ... finished in " << duration.count() << " s.\n\n\n";
    /************************************************************************/
    return 0;
}


/* Functions */
void readAtoms(string &file) {
    cout << "\n\n:: reading atom configuration...";
    string line;
    ifstream myfile;
    myfile.open(file.c_str(), ios::in);
    while (myfile.peek() != EOF) {
        getline(myfile, line);
        N++;
        // cout << "\natom " << N << " has charge " << line;
    }
    myfile.close();
    cout << "\n -> " << N << " atoms read.";
    L = (N - 1) * distance_atoms;
    cout << " Molecule length L = " << L << ".";
}

void header(string &file) {
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
