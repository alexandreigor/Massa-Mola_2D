#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
//#include <GL/glut.h>
#include <GL/freeglut.h>
#include "tipos.h"
#include "vetor.h"
#define NUM_M_X 14 // Número de massas
#define NUM_M_Y 14
#define DEBUG 0

int viewport_x = 30, viewport_y = 30;

const float L_X = 100, L_Y = 100; // Largura do "array" de massas-molas
massaSimples massas[NUM_M_Y][NUM_M_X];
molaSimples molasH[NUM_M_Y][NUM_M_X];
molaSimples molasV[NUM_M_Y][NUM_M_X];
float tempo = 400;
vetor forcas[NUM_M_Y][NUM_M_X];
FILE *arquivo;
char *nome_arq = "saida.txt";
unsigned long num_passos;
unsigned long passo;
double tol_min, tol_min_3, dt, dt2, dth, dt2h, fator_queda, f_max;
const double coef_arrasto = 0E-2;
estado estadoExecucao;
double mMax, mMin, kMax, kMin;

void init();
void singleStep();
void evolve();
void desenha();
void desenhaMassa(double x, double y, float m, float lMassa);
void desenhaMola(float l, float r, float a, float k, float y, float lMassa);
void keyboard(unsigned char key, int x, int y);
void openAndInitFile();
void writeStepInFile(double U, double K);
void closeAndPlotFile();
void finalize();
void finalize2();
inline double fMod(double x, double y);
inline int fDiv(double x, double y);
inline double menorDistancia(double num1, double num2, double l);
inline vetor menorDistanciaV(vetor v1, vetor v2, double lx, double ly);

int main(int argc, char **argv) {

	printf("Tamanho de char: %li\n\n", sizeof(char));

	init();
	atexit(finalize);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Sistema massa-mola");

	glClearColor(1.0, 1.0, 1.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	float borders = 0.0;
	gluOrtho2D(-borders * viewport_x, (1 + borders) * viewport_x, -borders
			* viewport_y, (1 + borders) * viewport_y);

	glutDisplayFunc(desenha);
	glutKeyboardFunc(keyboard);
	glutCloseFunc(finalize2);
	//glutIdleFunc(singleStep);
	glutIdleFunc(NULL);
	estadoExecucao = PAUSA;

	glutMainLoop();
	//*/

	return 0;
}

void openAndInitFile() {
	arquivo = fopen(nome_arq, "w");
	fprintf(arquivo, "t\tU\tK\tE\t\n");
}

void reverterTempo() {
	//	register int i, j;
	//	for (i = 0; i < NUM_M_Y; i++) {
	//			for (j = 0; j < NUM_M_X; j++) {
	//				massas[i][j].v = multEscalar(-1,massas[i][j].v);
	//				//massas[i][j].a = multEscalar(-1,massas[i][j].a);
	//			}
	//	}
	dt = -dt;
	dth = -dth;
}

// Inicializa as posições, velocidades, acelerações das massas e
// disposições das molas
void init() {
	num_passos = 200 * tempo;
	tol_min = L_X / (NUM_M_X * 50);
	tol_min_3 = 3 * tol_min;
	dt = tempo / num_passos;
	dt2 = dt * dt;
	dth = dt / 2;
	dt2h = dt2 / 2;
	fator_queda = 5 / tol_min;
	f_max = 7.0;

	const double m = 1.0, k = 0.01, dx = L_X / NUM_M_X, a = 5.0, dy = L_Y
			/ NUM_M_Y;
	vetor p_total, v_cm;
	double m_total = 0.0;
	register int i, j;

	setZero(&p_total);

	for (i = 0; i < NUM_M_Y; i++) {
		for (j = 0; j < NUM_M_X; j++) {
			massas[i][j].m = m;
			massas[i][j].r = newVetor(j * dx, i * dy);
			setZero(&(massas[i][j].v));
			setZero(&(massas[i][j].a));
			molasH[i][j].k = k;
			molasH[i][j].a = a;
			molasV[i][j].k = k;
			molasV[i][j].a = a;
		}
	}

	massas[0][0].r = newVetor(0.5 * dx, 0);
	massas[0][1].r = newVetor(0.4 * dx, 0.5 * dy);

	mMax = mMin = massas[0][0].m;
	kMax = kMin = molasH[0][0].k;

	for (i = 0; i < NUM_M_Y; i++) {
		for (j = 0; j < NUM_M_X; j++) {
			p_total = add(p_total, multEscalar(massas[i][j].m, massas[i][j].v));
			m_total += massas[i][j].m;

			if (massas[i][j].m > mMax) {
				mMax = massas[i][j].m;
			}
			if (massas[i][j].m < mMin) {
				mMin = massas[i][j].m;
			}

		}
	}

	v_cm = divEscalar(m_total, p_total);
	vprint("p_total", p_total);
	vprint("v_cm", v_cm);

	setZero(&p_total);
	vprint("p_total", p_total);
	/*
	 for (i = 0; i < NUM_M_Y; i++) {
	 for (j = 0; j < NUM_M_X; j++) {
	 massas[i][j].v.x -= v_cm.x;
	 massas[i][j].v.y -= v_cm.y;
	 p_total = add(p_total, multEscalar(massas[i][j].m, massas[i][j].v));
	 }
	 }
	 vprint("p_total", p_total);
	 */
	openAndInitFile();

}

inline void writeStepInFile(double U, double K) {
	fprintf(arquivo, "%f\t%f\t%f\t%f\t\n", passo * fabs(dt), U, K, (U + K));
}

void atualizarPosicoes() {
	register int i, j;
	for (i = 0; i < NUM_M_Y; i++) {
		for (j = 0; j < NUM_M_X; j++) {
			massas[i][j].r = add(massas[i][j].r, add(multEscalar(dt,
					massas[i][j].v), multEscalar(dt2h, massas[i][j].a)));

			//massas[i][j].r.x = fMod(massas[i][j].r.x, L_X);
			//massas[i][j].r.y = fMod(massas[i][j].r.y, L_Y);
		}
	}
}

inline void printDados(printable p) {
	register int i, j;
	char *alvo = nameOf(p);
	printf("\nImprimindo %s:\n", alvo);
	for (i = 0; i < NUM_M_Y; i++) {
		for (j = 0; j < NUM_M_X; j++) {
			vetor alvoImpressao;
			switch (p) {
			case FORCAS:
				alvoImpressao = forcas[i][j];
				break;
			case POSICOES:
				alvoImpressao = massas[i][j].r;
				break;
			case VELOCIDADES:
				alvoImpressao = massas[i][j].v;
				break;
			case ACELERACOES:
				alvoImpressao = massas[i][j].a;
				break;
			default:
				printf(
						"Argumento 'printable p = %i' inválido em printDados!\n",
						p);
				return;
			}
			printf("%s\t", toString(alvoImpressao));
		}
		printf("\n");
	}
}

double atualizarVelocidadeEAceleracaoECinetica() {
	double K = 0;
	register int i, j;
	for (i = 0; i < NUM_M_Y; i++) {
		for (j = 0; j < NUM_M_X; j++) {
			massas[i][j].v = add(massas[i][j].v, multEscalar(dth, add(
					forcas[i][j], massas[i][j].a)));
			massas[i][j].a = forcas[i][j];
			K += (massas[i][j].m) * sqrLength(massas[i][j].v) * (0.5);
		}
	}
	if (DEBUG)
		printDados(FORCAS);
	return K;
}

double calcularForcasEPotencial() {
	double U = 0.0;
	register int i, j;
	vetor dF;
	double mod_dF;

	for (i = 0; i < NUM_M_Y; i++) {
		for (j = 0; j < NUM_M_X; j++) {
			setZero(&forcas[i][j]);
		}
	}
	for (i = 0; i < NUM_M_Y; i++) {

		int nextY = (i + 1) % NUM_M_Y;

		for (j = 0; j < NUM_M_X; j++) {

			int nextX = (j + 1) % NUM_M_X;

			vetor dr_X = menorDistanciaV(massas[i][nextX].r, massas[i][j].r,
					L_X, L_Y);
			vetor dr_Y = menorDistanciaV(massas[nextY][j].r, massas[i][j].r,
					L_X, L_Y);

			double mod_dr_X = length(dr_X);
			double mod_dr_Y = length(dr_Y);
			double dxH = mod_dr_X - molasH[i][j].a;
			double dxV = mod_dr_Y - molasV[i][j].a;

			/* Forças entre massas lado a lado */
			mod_dF = molasH[i][j].k * (dxH);
			dF = newVetorByDirection(mod_dF, dr_X);
			forcas[i][j] = add(forcas[i][j], divEscalar(massas[i][j].m, dF));
			forcas[i][nextX] = sub(forcas[i][nextX], divEscalar(
					massas[i][nextX].m, dF));

			/* Forças entre massas acima/abaixo */
			mod_dF = molasV[i][j].k * (dxV);
			dF = newVetorByDirection(mod_dF, dr_Y);
			forcas[i][j] = add(forcas[i][j], divEscalar(massas[i][j].m, dF));
			forcas[nextY][j] = sub(forcas[nextY][j], divEscalar(
					massas[nextY][j].m, dF));

			/* Efeito de arrasto (dissipação de enrgia pelo ar) */
			forcas[i][j] = sub(forcas[i][j], multEscalar(coef_arrasto,
					massas[i][j].v));

			U += 0.5 * (molasH[i][j].k * pow(dxH, 2) + molasV[i][j].k * pow(
					dxV, 2.0));

			// TODO: Tratar os choques em duas dimensões
			/*
			 if (mod_dr_X < tol_min_3) {
			 float expon = exp(-fator_queda * fabs(mod_dr_X)) * f_max;
			 forcas[0][j].x -= fator_queda * expon;
			 forcas[0][nextX].x += fator_queda * expon;

			 U += expon;
			 }
			 */
		}
	}
	return U;
}

void singleStep() {
	if (estadoExecucao == RODANDO) {
		double U = 0.0, K = 0.0;

		if (DEBUG && !(passo % 15)) {
			estadoExecucao = PAUSA;
			glutIdleFunc(NULL);
		}

		atualizarPosicoes();
		U = calcularForcasEPotencial();
		K = atualizarVelocidadeEAceleracaoECinetica();

		writeStepInFile(U, K);

		if (!(passo % 50)) {
			desenha();
		}
		printf("passo=%li\tU=%f\tK=%f\tE=%f\n", passo, U, K, (U + K));
		passo++;
	} else {
		closeAndPlotFile();
	}

}

void closeAndPlotFile() {
	//*
	char *nome_arq_plot1 = "E_molas.gnu";
	char comando[150] = "";

	fclose(arquivo);
	arquivo = fopen(nome_arq_plot1, "w");
	fprintf(arquivo, "set encoding iso_8859_1\n"
		"plot 'saida.txt' using 1:2 title 'U' with lines, \\\n"
		"     'saida.txt' using 1:3 title 'K' with lines, \\\n"
		"     'saida.txt' using 1:4 title 'E' with lines\n"
		"#set terminal png\n"
		"#set output 'energias.png'\n"
		"#replot\n");
	fclose(arquivo);

	sprintf(comando, "gnome-terminal -e \"gnuplot '%s' - \" & ", nome_arq_plot1);
	printf(comando);
	system(comando);

}

void evolve() {
	for (passo = 0; passo < num_passos; passo++) {
		singleStep();
	}
}
void finalize() {
	printf("exit(0);");
	closeAndPlotFile();
	exit(0);
}

void finalize2() {
	printf("exit(2);");
	estadoExecucao = TERMINAR;

}

void desenha(void) {
	int i, j;
	float lMassa = L_X / (NUM_M_X * 15); // Proporcional à distância de interação do "potencial de contato"
	glClear(GL_COLOR_BUFFER_BIT);

	glColor3f(0.5, 0.5, 1.0);

	glPushMatrix();
	glScalef(viewport_x / L_X, viewport_y / L_Y, 1.0);
	glTranslated(0.0, L_Y, 0.0);
	glRotated(180, 1.0, 0.0, 0.0);
	//glTranslated(0.0, -L_Y/2 , 0.0);

	/* Primeira massa de cor diferente
	 glColor4f(0.3, 0.3, 0.3, 0.5);
	 desenhaMassa(massas[0][0].r.x, y, massas[0][0].m, lMassa);
	 // */

	glBegin(GL_LINE_STRIP);
	{
		glVertex2f(0.0, 0.0);
		glVertex2f(L_X, 0.0);
		glVertex2f(L_X, L_Y);
		glVertex2f(0.0, L_Y);
		glVertex2f(0.0, 0.0);
	}
	glEnd();

	glColor3f(0.5, 0.5, 1.0);
	for (i = 0; i < NUM_M_Y; i++) {
		int nextY = (i + 1) % NUM_M_Y;

		for (j = 0; j < NUM_M_X; j++) {
			int nextX = (j + 1) % NUM_M_X;
			double x0, y0;
			double xl, yl;
			double xd, yd;

			x0 = fMod(massas[i][j].r.x, L_X);
			y0 = fMod(massas[i][j].r.y, L_Y);

			xl = fMod(massas[i][nextX].r.x, L_X);
			yl = fMod(massas[i][nextX].r.y, L_Y);

			xd = fMod(massas[nextY][j].r.x, L_X);
			yd = fMod(massas[nextY][j].r.y, L_Y);

			if (!(i + j)) {
				glColor3f(1.0, 0.2, 0.2);
			} else if (!i) {
				glColor3f(0.8, 0.4, 0.0);
			} else if (!j) {
				glColor3f(0.4, 0.8, 0.0);
			} else {
				glColor3f(0.2, 0.2, 0.0);

			}

			desenhaMassa(x0, y0, massas[i][j].m, lMassa);

			vetor dr_X = menorDistanciaV(massas[i][nextX].r, massas[i][j].r,
					L_X, L_Y);
			vetor dr_Y = menorDistanciaV(massas[nextY][j].r, massas[i][j].r,
					L_X, L_Y);

			glLineWidth(0.1);
			glBegin(GL_LINES);
			{
				glColor3f(0.5, 1.0, 1.0);
				glVertex2d(x0, y0);
				glVertex2d(x0 + dr_X.x, y0 + dr_X.y);

				glVertex2d(xl, yl);
				glVertex2d(xl - dr_X.x, yl - dr_X.y);

				glColor3f(1.0, 0.0, 1.0);
				glVertex2d(x0, y0);
				glVertex2d(x0 + dr_Y.x, y0 + dr_Y.y);

				glVertex2d(xd, yd);
				glVertex2d(xd - dr_Y.x, yd - dr_Y.y);
			}
			glEnd();

		}
	}

	glPopMatrix();

	glFlush();

}

void desenhaCirculo(double xc, double yc, double radius) {
	register char i;
	char num = 16;
	double da = 2.0 * M_PI / (num - 1);
	glPushMatrix();
	glTranslatef(xc, yc, 0.0);
	glBegin(GL_POLYGON);
	for (i = 0; i < num; i++) {
		float ang = i * da;
		glVertex2d(radius * cos(ang), radius * sin(ang));
	}
	glEnd();
	glPopMatrix();
}

void desenhaMassa(double xc, double yc, float m, float lMassa) {
	double rMax = 0.7;
	double rMin = rMax / 2;
	double radius = rMin;//rMin + (rMax - rMin) * (m - mMin) / (mMax - mMin);

	//if (mMax > mMin) radius = rMin + (rMax - rMin) * (m - mMin) / (mMax - mMin);

	desenhaCirculo(xc, yc, radius);
	desenhaCirculo(xc + L_X, yc + L_Y, radius);
	desenhaCirculo(xc + L_X, yc, radius);
	desenhaCirculo(xc, yc + L_X, radius);
	desenhaCirculo(xc - L_X, yc - L_Y, radius);
	desenhaCirculo(xc - L_X, yc, radius);
	desenhaCirculo(xc, yc - L_X, radius);
}

void desenhaMola(float l, float r, float a, float k, float y, float lMassa) {
	float dxMola = 0.5 * lMassa;
	float w, ang, rad;
	float xr = l + lMassa + dxMola, xl = r - lMassa - dxMola;
	float Dx = xl - xr;

	ang = (2 * a) * M_PI;
	rad = 0.5;

	glBegin(GL_LINE_STRIP);

	glVertex2f(l + lMassa, y + 0.5);
	glVertex2f(l + lMassa + dxMola, y + 0.5);

	for (w = 0; w <= ang; w += 0.1) {
		float x = xr + Dx * w / ang;
		glVertex2f(x, y + 0.5 + rad * sin(w));
	}

	glVertex2f(r - lMassa - dxMola, y + 0.5);
	glVertex2f(r - lMassa, y + 0.5);

	glEnd();

}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		printf("<ESC>\n");
		finalize();
		return;
	case 'r': {
		char old = estadoExecucao;
		estadoExecucao = PAUSA;
		reverterTempo();
		estadoExecucao = old;
		break;
	}
	case ' ':
		if (estadoExecucao == RODANDO) {
			glutIdleFunc(NULL);
			estadoExecucao = PAUSA;
		} else {
			glutIdleFunc(singleStep);
			estadoExecucao = RODANDO;
		}
		break;
	case 'p':
		printDados(POSICOES);
		return;
	case 'f':
		printDados(FORCAS);
		return;
	case 'v':
		printDados(VELOCIDADES);
		return;
	case 'a':
		printDados(ACELERACOES);
		return;
	default:
		printf("<default>\n");
		return;
	};
	glutPostRedisplay();
}

inline double fMod(double x, double y) {
	double ret = fmod(x, y);
	return (ret < 0) ? (ret + y) : ret;
}

inline int fDiv(double x, double y) {
	return floor(x / y);
}

inline double menorDistancia(double num1, double num2, double l) {
	double diff, diffAbs, aux, auxAbs;

	diff = num1 - num2 + l;
	diffAbs = fabs(diff);

	aux = num1 - num2;
	auxAbs = fabs(aux);

	if (auxAbs < diffAbs) {
		diff = aux;
		diffAbs = auxAbs;
	}

	aux = num1 - num2 - l;
	auxAbs = fabs(aux);

	if (auxAbs < diffAbs) {
		diff = aux;
		diffAbs = auxAbs;
	}

	return diff;

}

inline vetor menorDistanciaV(vetor v1, vetor v2, double lx, double ly) {
	vetor menorDist = newVetor(0, 0);

	menorDist.x = menorDistancia(v1.x, v2.x, lx);
	menorDist.y = menorDistancia(v1.y, v2.y, ly);

	return menorDist;
}
