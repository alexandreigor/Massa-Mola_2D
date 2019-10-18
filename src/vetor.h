#ifndef VETOR_H_
#define VETOR_H_

typedef struct _vetor vetor;

vetor newVetor(double x, double y);
vetor newVetorByDirection(double modulo, vetor direcao);
vetor add(const vetor vetor1, const vetor vetor2);
vetor sub(const vetor vetor1, const vetor vetor2);
vetor multEscalar(const double k, const vetor _vetor);
vetor divEscalar(const double k, const vetor _vetor);
double prodEscalar(const vetor vetor1, const vetor vetor2);
double sqrLength(const vetor vetor);
double length(const vetor vetor);
void setZero(vetor *_vetor);
void vprint(char *nome, vetor _vetor);
char* toString(vetor _vetor);

struct _vetor {
	double x;
	double y;
};

#endif /* VETOR_H_ */
