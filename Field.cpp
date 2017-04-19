#include"Field.h"

Field::Field(int width, int height, double h_u, int pml)
:WIDTH(width),
HEIGHT(height),
H_U(h_u),
H_S(1),
N_X(width/h_u),
N_Y(height/h_u),
N_PML(pml), 
N_PX(width/h_u + 2*pml),
N_PY(height/h_u + 2*pml), 
N_CELL( (width/h_u + 2*pml) * (height/h_u + 2*pml))
{
/*
	cout << "Nx = "<<N_X << endl;
	cout << "Ny = "<<N_Y << endl;
	cout << "N_PML = "<<N_PML << endl;
	cout << "N_CELL = "<<N_CELL << endl;
	*/
};


//B-PML‚É‚¨‚¯‚é“±“d—¦ƒÐ
double Field::sigmaX(const int &i, const int &j){
	//ƒ¢x = h = 1
	if(i<N_PML) 
		return 1.0*(N_PML - i)/N_PML;

	else if(N_PML<=i && i<(N_X+N_PML)) 
		return 0.0;

	else
		return 1.0*(i - (N_PX - N_PML - 1))/N_PML;
	//‚±‚ê‚¾‚Æˆês‚Å‘‚¯‚é‚¯‚Ç, ‰Â“Ç«‚ª‚æ‚­‚È‚¢
	//	int k = i%(N_X+N_PML);
	//return (N_PML*H_S - k)/(N_PML*H_S);
}

double Field::sigmaY(const int &i, const int &j){
	//ƒ¢y = h = 1
	if(j<N_PML) 
		return 1.0*(N_PML - j)/N_PML;

	else if(j>=N_PML && j<(N_Y+N_PML))
		return 0.0;

	else
		return 1.0*(j - (N_PY - N_PML - 1)*H_S)/(N_PML*H_S);
}