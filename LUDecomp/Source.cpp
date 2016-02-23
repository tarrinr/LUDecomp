#include "Twin.h"
#include <sstream>

std::vector<std::vector<double>> getA(Twin&);
std::vector<double> getb(Twin&, const int&);
void dmat(Twin&, const std::vector<std::vector<double>>&);
void dvec(Twin&, const std::vector<double>&);
std::vector<double> ludecomp(const std::vector<std::vector<double>>&);

int main() {

	Twin t("LU Decomposition");

	while (true) {

		std::vector<std::vector<double>> A = getA(t);
		std::vector<double> b = getb(t, A.size());

		std::vector<std::vector<double>> aug(A);

		for (int i = 0; i < int(aug.size()); i++)
			aug[i].push_back(b[i]);

		std::vector<double> x = ludecomp(aug);
		t.println("Vector x:");
		t.println();

		dvec(t, x);
		t.println("Continue? [y/n]");
		t.display();

		std::string input;
		std::getline(std::cin, input);

		if (input != "y" && input != "Y") break;
	}

	return EXIT_SUCCESS;
}

std::vector<std::vector<double>> getA(Twin& t) {

	std::vector<std::vector<double>> A;

	bool enter = true;
	while (enter) {

		while (true) {
			t.println("Enter the first row of matrix A, in Ax = b.");
			t.println("A must be a square, full ranking matrix.");
			t.println("Example: 1 2 3 4");
			t.display();

			double in;
			std::string input;
			std::getline(std::cin, input);

			if (input == "") {
				t.println("Try again");
				t.println();
			}
			else {
				std::stringstream iss(input);
				A.push_back(std::vector<double>());
				while (iss >> in) A[0].push_back(in);
				break;
			}
		}

		int m = 1;

		while (m < int(A[0].size())) {

			dmat(t, A);

			t.println("Enter the next row of the matrix.");
			t.display();

			double in;
			std::string input;
			std::getline(std::cin, input);
			std::stringstream iss(input);

			A.push_back(std::vector<double>());
			while (iss >> in) A[m].push_back(in);

			if (A[m].size() == A[0].size())
				m++;
			else {
				A.pop_back();
				t.println("Invalid row.");
				t.println();
			}
		}

		t.println("Matrix A:");
		t.println();
		dmat(t, A);
		t.println("Save and continue? [y/n]");
		t.display();

		std::string input;
		std::getline(std::cin, input);

		if (input != "y" && input != "Y") {
			A.clear();
		}
		else enter = false;
	}
	return A;
}

std::vector<double> getb(Twin& t, const int& m) {

	std::vector<double> b;

	while (true) {

		t.println("Enter the vector b, in Ax = b.");
		t.println("Example: 1 2 3 4");

		t.display();

		double in;
		std::string input;
		std::getline(std::cin, input);

		std::stringstream iss(input);
		while (iss >> in) b.push_back(in);

		if (int(b.size()) == m) {
			t.println("Vector b:");
			t.println();
			dvec(t, b);
			t.println("Save and continue? [y/n]");
			t.display();

			std::getline(std::cin, input);
			if (input == "y" || input == "Y") {
				break;
			}
			else b.clear();
		}
		else {
			b.clear();
			t.println("Invalid vector.");
			t.println();
		}
	}

	return b;
}

void dmat(Twin& t, const std::vector<std::vector<double>>& mat) {

	t.println();

	for (int i = 0; i < int(mat.size()); i++) {

		t.print("[ ");

		for (int j = 0; j < int(mat[0].size()); j++) {
			t.print(mat[i][j]);
			if (j < int(mat[0].size()) - 1) t.print(",");
			t.print(" ");
		}

		t.print("]");
		t.println();

	}
}

void dvec(Twin& t, const std::vector<double>& vec) {

	for (auto&& i : vec) {
		t.println("[ ");
		t.print(i);
		t.print(" ]");
	}

	t.println();
}

std::vector<double> ludecomp(const std::vector<std::vector<double>>& aug) {

	int n = aug.size();

	std::vector<std::vector<double>> l, u(aug);
	std::vector<double> x(n), y(n), b;

	for (int i = 0; i < n; i++) {
		b.push_back(u[i].back());
		u[i].pop_back();
		l.push_back(std::vector<double>(n, 0));
		l[i][i] = 1;
	}

	for (int i = 0; i < n; i++) {

		double maxEl = abs(u[i][i]);
		int maxRow = i;
		for (int k = i + 1; k < n; k++) {
			if (abs(u[k][i]) > maxEl) {
				maxEl = abs(u[k][i]);
				maxRow = k;
			}
		}

		for (int k = i; k < n; k++) {
			double tmp = u[maxRow][k];
			u[maxRow][k] = u[i][k];
			u[i][k] = tmp;
		}

		double tmp = b[maxRow];
		b[maxRow] = b[i];
		b[i] = tmp;

		for (int k = i + 1; k < n; k++) {
			l[k][i] = u[k][i] / u[i][i];

			for (int j = i; j < n; j++) u[k][j] += -l[k][i] * u[i][j];

		}
	}

	for (int i = 0; i < n; i++) {
		y[i] = b[i] / l[i][i];
		for (int k = 0; k < n; k++) {
			b[k] -= l[k][i] * y[i];
		}
	}

	for (int i = n - 1; i >= 0; i--) {
		x[i] = y[i] / u[i][i];
		for (int k = i - 1; k >= 0; k--) {
			y[k] -= u[k][i] * x[i];
		}
	}

	return x;
	
}