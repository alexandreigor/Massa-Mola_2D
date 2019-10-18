#include "vetor.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

inline vetor newVetor(double x, double y) {
	vetor _vetor;
	_vetor.x = x;
	_vetor.y = y;

	return _vetor;
}

inline vetor newVetorByDirection(double modulo, vetor direcao) {
	double fat = modulo / length(direcao);

	return multEscalar(fat, direcao);
}

inline vetor add(const vetor vetor1, const vetor vetor2) {
	vetor vetorResp;
	vetorResp.x = vetor1.x + vetor2.x;
	vetorResp.y = vetor1.y + vetor2.y;

	return vetorResp;
}
inline vetor sub(const vetor vetor1, const vetor vetor2) {
	vetor vetorResp;
	vetorResp.x = vetor1.x - vetor2.x;
	vetorResp.y = vetor1.y - vetor2.y;

	return vetorResp;
}
inline double sqrLength(const vetor _vetor) {
	return prodEscalar(_vetor, _vetor);
}

inline double length(const vetor _vetor) {
	return sqrt(prodEscalar(_vetor, _vetor));
}

inline vetor multEscalar(const double k, const vetor _vetor) {
	vetor vetorResp;
	vetorResp.x = k * _vetor.x;
	vetorResp.y = k * _vetor.y;
	return vetorResp;
}
inline vetor divEscalar(const double k, const vetor _vetor) {
	vetor vetorResp;
	vetorResp.x = _vetor.x / k;
	vetorResp.y = _vetor.y / k;
	return vetorResp;
}

inline double prodEscalar(const vetor vetor1, const vetor vetor2) {
	return vetor1.x * vetor2.x + vetor1.y * vetor2.y;
}

inline void setZero(vetor *vetor1) {
	(*vetor1).x = 0.0;
	(*vetor1).y = 0.0;
}

inline void vprint(char *nome, vetor _vetor) {
	printf("%s=( %f , %f )\n", nome, _vetor.x, _vetor.y);
}

inline char* toString(vetor _vetor) {
	char *str;
	str = (char*) malloc(sizeof(char[50]));
	sprintf(str, "(%03.5f, %03.5f)", _vetor.x, _vetor.y);
	return str;
}

