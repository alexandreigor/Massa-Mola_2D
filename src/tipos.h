/*
 * tipos.h
 *
 *  Created on: 17/03/2011
 *      Author: igor
 */
#include "vetor.h"

#ifndef TIPOS_H_
#define TIPOS_H_

typedef struct _molaS molaSimples;
typedef struct _massaS massaSimples;
typedef enum _estado estado;
typedef enum _printable printable;

enum _estado {
	TERMINAR = (char) 0, RODANDO = (char) 1, PAUSA = (char) 2
};

enum _printable {
	FORCAS = (char) 0,
	POSICOES = (char) 1,
	VELOCIDADES = (char) 2,
	ACELERACOES = (char) 3,
};

struct _molaS {
	double k;
	double a;
};

struct _massaS {
	double m; // Massa
	vetor r; // Posição
	vetor v; // Velocidade
	vetor a; // Aceleração
};

char* nameOf(printable p) {
	switch (p) {
	case FORCAS:
		return "FORCAS";
	case POSICOES:
		return "POSICOES";
	case VELOCIDADES:
		return "VELOCIDADES";
	case ACELERACOES:
		return "ACELERACOES";
	default:
		printf("Argumento 'printable p = %i' inválido em nameOf!\n", p);
		return "null";
	}
}

#endif /* TIPOS_H_ */
