#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>

#include <iomanip>
#include <sstream>

using namespace std;
using namespace arma;

void get_atom(char line[], char atom[]){
    int atom_index = 0;
    for(int i=12; i<16; i++){
        if(isalnum(line[i])){
            atom[atom_index++] = line[i];
        }
    }
    atom[atom_index] = '\0';
}

bool correct_atom(char atom[]){
    return strcmp(atom,"C")==0 || strcmp(atom,"CA")==0 || strcmp(atom,"N")==0 || strcmp(atom,"O")==0;
}

double get_double(string line, int start, int length){
    double x = atof(line.substr(start, length).c_str());
    return x;
}

int atom_index(char atom[]){
    if(strcmp(atom,"N")==0){
        return 0;
    } else if(strcmp(atom,"CA")==0){     
        return 1;
    } else if(strcmp(atom,"C")==0){
        return 2;
    } else if(strcmp(atom,"O")==0){
        return 3;
    }

    cout << "Unrecognized atom " << atom << endl;
}

int get_index(char* atom, char atoms[][5], int n_atoms){
    for(int i=0; i<n_atoms; i++){
        if(strcmp(atom, atoms[i]) == 0){
            return i;
        }
    }
    return -1;
}

void read_number(char line[], int from, int length, double *result){
    char num[length+1];
    memcpy(num, &line[from], length);
    num[length] = '\0';
    *result = atof(num);
}

void read_file(const char *file, mat &points, rowvec &g_atom_found, char atoms[][5], int n_atoms){
    ifstream input(file);
    char buf[81];
    char c_aa[6], l_aa[6], atom[5];    
    int aa_index = -1;
    
    c_aa[5] = '\0';
    l_aa[5] = '\0';
    atom[4] = '\0';
    mat aa_atoms = zeros<mat>(3, n_atoms);
    //mat empty = zeros<mat>(n_atoms,3);
    rowvec l_atom_found = zeros<rowvec>(n_atoms);
    //rowvec empty_vec = zeros<rowvec>(n_atoms);
    
    //cout << "Reading file " << file << endl;
    
    while(input.good()){
        input.getline(buf, 81);
                
        // Change Amino Acid
        memcpy(l_aa, &buf[22], 5);
        //cout << c_aa << endl;
        
        if(strcmp(c_aa, l_aa) != 0){
            //cout << "Change from" << c_aa << "to" << l_aa << " " << strcmp(c_aa, l_aa) << endl;
            memcpy(c_aa, &l_aa[0], 5);            
            aa_index++;
            //cout << aa_index << '\n';
            if(aa_index != 0){
                points = join_horiz(points, aa_atoms);
                
                //cout << points << '\n';
                
                g_atom_found = join_horiz(g_atom_found, l_atom_found);
            }
			
            l_atom_found = zeros<rowvec>(n_atoms);
            aa_atoms = zeros<mat>(3, n_atoms);
            
        }
        
        get_atom(buf, atom);
        int atom_index = get_index(atom, atoms, n_atoms);
        
        if(atom_index != -1){
            double x, y, z;
            read_number(buf,30, 8,&x);
            read_number(buf,38, 8,&y);
            read_number(buf,46, 8,&z);
        
            colvec point = zeros<colvec>(3);
            point(0) = x;
            point(1) = y;
            point(2) = z;
        
        
            aa_atoms.col(atom_index) = point;
            
            l_atom_found(atom_index) = 1.0;
        }
    }
    
    points = join_horiz(points, aa_atoms);
    
    g_atom_found = join_horiz(g_atom_found, l_atom_found);
}



bool read_files(const char* file1, const char* file2, mat &points1, mat &points2, char atoms[][5], int n_atoms){
    mat pos1, pos2;
    rowvec atom_found1, atom_found2;
    rowvec total_atoms_found;
    uvec logic_vec;
	
    thread thread1(read_file, file1, ref(pos1), ref(atom_found1), atoms, n_atoms);
    thread thread2(read_file, file2, ref(pos2), ref(atom_found2), atoms, n_atoms);
    
    thread1.join();
    thread2.join();
    //~ 
    //~ read_file( file1, pos1, atom_found1, atoms, n_atoms);
    //~ read_file( file2, pos2, atom_found2, atoms, n_atoms);
    
    
    // for(int i=0; i< pos1.size(); i++){
    //     cout << "Amino acid " << i << endl;
    //     for(int j=0; j<pos1[i].size(); j++){
    //         cout << "\t";
    //         for(int k=0; k<3; k++){
    //             cout << pos1[i][j][k] << " ";
    //         }
    //         cout << endl;
    //     }
    // }
    
    total_atoms_found = atom_found1 % atom_found2;
    
    logic_vec = find(total_atoms_found==1);
	
	points1 = pos1.cols(logic_vec);
    
    points2 = pos2.cols(logic_vec);
    
    
    
    if(atom_found1.n_elem != atom_found2.n_elem){
        cout << "The PDB files have different number of Amino acids" << endl;
        return false;
    } 
    //~ 
    //~ int index1=0;
    //~ int index2=0;
    //~ for(int i=0; i<atom_found1.size(); i++){
        //~ for(int j=0; j<n_atoms; j++){
			//~ 
            //~ if(atom_found1[i][j] && atom_found2[i][j]){
                //~ points1 = join_horiz(points1, pos1.col(index1));
                //~ points2 = join_horiz(points2, pos2.col(index2));
                //~ //cout << i*n_atoms+j << ' ';
            //~ }
            //~ else{
				//~ cout << i << ' '<< j << ' ' << atom_found1[i][j] <<'\n';
				//~ cout << i << ' '<< j << ' ' << atom_found2[i][j] <<'\n';
			//~ }
			
			//~ 
			//~ index1+=atom_found1[i][j];
			//~ 
			//~ index2+=atom_found2[i][j];
        //~ }
    //~ }
    
    
    return true;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


void CalculateRmsd(mat P1, mat P2, double* rmsd)
{
	int n_points = P1.n_cols;
	
	mat Diff = P1-P2;
	
	*rmsd = sqrt(sum(sum((Diff % Diff)))/n_points);
}

// P and Q must contain the points in the columns and coordinates in the rows i.e. number of dimensions = number of rows, number of points = number of columns
mat Kabsch(mat P, mat Q, double* rmsd)
{
	mat Middle_mat;
	
	mat P_translated;
	mat Q_translated;
	mat P_rotated_translated;
	mat P_finish;
	
	mat U;
	vec s;
	mat V;
	
	mat R;
	
	colvec P_centroid;
	colvec Q_centroid;
	rowvec v_ones;
	
	mat Covariance_matrix;
	
	int d;
	int n_points;
	
	// Calculate the number of datapoints
	n_points = P.n_cols;
	
	//Calculate the coordinates of the centres of mass
	P_centroid=sum(P,1)/n_points;
	
	Q_centroid=sum(Q,1)/n_points;
	
	v_ones = ones<rowvec>(n_points);
	
	//Translate the points so that the cetre of mass coincides with the origin
	P_translated = P - P_centroid*v_ones;
	
	Q_translated = Q - Q_centroid*v_ones;
	
	//Calculate the covariance matrix	
	Covariance_matrix=P_translated*Q_translated.t()/n_points;
	
	//Calculate svd
	svd(U,s,V,Covariance_matrix);
	
	//Ensure a right-handed coordinate system
	d = sgn(det( V*U.t() ));
	
	Middle_mat = eye<mat>(3,3);
	
	Middle_mat(2,2) =  d;
	
	//Calculate the rotation matrix
	R=V*Middle_mat*U.t();
	
	//Rotate P
	P_rotated_translated = R*P_translated;
	
	//Translate P to the centre of mass of Q
	P_finish = P_rotated_translated + Q_centroid*v_ones;
	
	//Calculate rmsd
	CalculateRmsd(P_rotated_translated, Q_translated, rmsd);
	
	
	//Done!
	return P_finish;
}

string format_number(double Number_to_format)
{
	int wspc_missing;
	
	string Result;
	
	string complete_string;
	
	ostringstream Convert;
	
	Convert << fixed << setprecision(3) << Number_to_format;
	
	Result = Convert.str();
	
	wspc_missing = 7 - Result.length();
	
	complete_string = string(wspc_missing, ' ') + Result;
	
	return complete_string;

}

void write_modified_file(const char *original_file, mat P_mat, const char *result_file)
{
	
	ifstream infile(original_file);
	
	ofstream pdb_modified;
	
	double x,y,z;
	
	int i;
	
	pdb_modified.open(result_file);
	
	string line;
	
	string x_string;
	
	string y_string;
	
	string z_string;
	
	i=0;
	
	while(getline(infile, line))
	{
		
		if (line.compare(13, 3, "CA ")==0 or line.compare(13, 3, "N  ")==0 or line.compare(13, 3, "C  ")==0 or line.compare(13, 3, "O  ")==0)
		
		{
			x = P_mat(0,i);
			
			y = P_mat(1,i);
			
			z = P_mat(2,i);
			
			x_string = format_number(x);
			
			y_string = format_number(y);
			
			z_string = format_number(z);
			
			line.replace(31,7,x_string);
			
			line.replace(39,7,y_string);
			
			line.replace(47,7,z_string);
			
			//~ cout << line << line.length() << '\n';
			
			pdb_modified << line << '\n';
			
			i++;
		}
	}
	
	pdb_modified.close();
}

int main(int argc, char** argv)
{
	colvec a;
	mat P1_mat;
	colvec b;
	mat P2_mat;
	
	mat P1_kabsch;
	
	mat Diff;
	
	int d;
	int len;
	
	double rmsd;
	
	char atoms[4][5];
	
	
	//strcpy(atoms[0],"CA\0");
	strcpy(atoms[0],"N");
	strcpy(atoms[1],"CA");
	strcpy(atoms[2],"C");
	strcpy(atoms[3],"O");
	
	const char *file1 = "p1.pdb";
	const char *file2 = "p2.pdb";
	const char *result_file1 = "p1_modif.pdb";
	const char *result_file2 = "p2_modif.pdb";
	
	
	//~ for(std::vector<vector<double> >::iterator it = p1.begin(); it != p1.end(); ++it) {
			//~ a=conv_to<colvec>::from((*it));
			//~ P1_mat = join_horiz(P1_mat, a);
	//~ }
	//~ 
	//~ for(std::vector<vector<double> >::iterator it = p2.begin(); it != p2.end(); ++it) {
			//~ b=conv_to<colvec>::from((*it));
			//~ P2_mat = join_horiz(P2_mat, b);
	//~ }
	

	read_files(file1, file2, P1_mat, P2_mat, atoms, 4);
	P1_kabsch=Kabsch(P1_mat, P2_mat , &rmsd);
	
	write_modified_file(file1, P1_kabsch, result_file1);
	write_modified_file(file2, P2_mat, result_file2);
	
	cout << rmsd << '\n';
	
	return 0;
}

