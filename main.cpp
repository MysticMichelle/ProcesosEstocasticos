#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_ESTADOS    20
#define MAX_DECISIONES 20
#define MAX_POLITICAS  10000
#define COSTO_INVALIDO 999999999.0   // Costo "blindado" para decisiones no aplicables

// ============================================================
//  Prototipos
// ============================================================
void   Portada();
void   menuPrincipal();
void   algoritmo1();
void   algoritmo2();

double leerValor();                  // Lee decimal o fraccion a/b

void leerCostos(int m, int k,
                double costos[MAX_ESTADOS][MAX_DECISIONES],
                int   aplicable[MAX_DECISIONES][MAX_ESTADOS]);

void leerMatrizTransicion(int decision, int m,
                          double P[MAX_ESTADOS][MAX_ESTADOS],
                          int    aplicable[MAX_ESTADOS]);

void construirMatrizPolitica(int m,
                             int    politica[MAX_ESTADOS],
                             double Pdec[MAX_DECISIONES][MAX_ESTADOS][MAX_ESTADOS],
                             double Ppol[MAX_ESTADOS][MAX_ESTADOS]);

int  resolverSistemaEstacionario(int m,
                                 double P[MAX_ESTADOS][MAX_ESTADOS],
                                 double pi[MAX_ESTADOS]);

double calcularEsperanza(int m,
                         int    politica[MAX_ESTADOS],
                         double pi[MAX_ESTADOS],
                         double costos[MAX_ESTADOS][MAX_DECISIONES]);

void imprimirMatriz(const char* titulo, int n,
                    double M[MAX_ESTADOS][MAX_ESTADOS]);

void imprimirSistemaEstacionario(int m,
                                 double P[MAX_ESTADOS][MAX_ESTADOS],
                                 double pi[MAX_ESTADOS]);

// ============================================================
//  Portada
// ============================================================
void Portada() {
    printf("*-*-*-*-  Procesos Estocasticos  -*-*-*-*\n");
    printf("          Proyecto Final I\n\n");
    printf("Profesora : Ada Ruth Cuellar Aguayo\n\n");
    printf("Integrantes:\n");
    printf("  Deciga Torres Melanie\n");
    printf("  De Yahve Heredia Oliver\n");
    printf("  Lopez Martinez Valeria\n");
    printf("  Moctezuma Isidro Michelle\n");
    printf("\nPresione Enter para continuar...\n");
    getchar();
    system("cls");
}

// ============================================================
//  Menu Principal
// ============================================================
void menuPrincipal() {
    int x, repetir = 1;
    do {
        system("cls");
        printf("\n========== MENU PRINCIPAL ==========\n");
        printf("1. Enumeracion Exhaustiva de Politicas\n");
        printf("2. Metodo de Programacion Lineal (pendiente)\n");
        printf("3. Salir\n");
        printf("=====================================\n");
        printf("Selecciona una opcion: ");
        scanf("%d", &x);
        getchar();

        switch (x) {
            case 1: algoritmo1(); break;
            case 2: algoritmo2(); break;
            case 3:
                printf("Adios, vuelva pronto.\n");
                repetir = 0;
                break;
            default:
                printf("Opcion invalida, intente de nuevo.\n");
        }
        if (repetir) {
            printf("\nPresione Enter para volver al menu...\n");
            getchar();
        }
    } while (repetir);
}

// ============================================================
//  leerValor: acepta "0.375" o "3/8" indistintamente
// ============================================================
double leerValor() {
    char buffer[64];
    scanf("%s", buffer);
    char *slash = strchr(buffer, '/');
    if (slash != NULL) {
        *slash = '\0';
        double num = atof(buffer);
        double den = atof(slash + 1);
        if (fabs(den) < 1e-15) {
            printf("  [!] Denominador cero, se usa 0.\n");
            return 0.0;
        }
        return num / den;
    }
    return atof(buffer);
}

// ============================================================
//  Leer costos.
//  Si la decision no aplica al estado => COSTO_INVALIDO
//  automaticamente (ya marcado en aplicable[][]).
// ============================================================
void leerCostos(int m, int k,
                double costos[MAX_ESTADOS][MAX_DECISIONES],
                int    aplicable[MAX_DECISIONES][MAX_ESTADOS]) {

    printf("\n--- Ingrese los costos C[estado][decision] ---\n");
    printf("    (Los pares no aplicables quedan en %.0f automaticamente)\n",
           COSTO_INVALIDO);

    for (int i = 0; i <= m; i++) {
        for (int d = 1; d <= k; d++) {
            if (aplicable[d][i]) {
                printf("  C[%d][%d] = ", i, d);
                costos[i][d] = leerValor();
            } else {
                costos[i][d] = COSTO_INVALIDO;
                printf("  C[%d][%d] = %.0f  [decision no aplicable, asignado automaticamente]\n",
                       i, d, COSTO_INVALIDO);
            }
        }
    }
    getchar();
}

// ============================================================
//  Leer matriz de transicion para una decision.
//  Pregunta estado por estado si aplica.
//  Si NO aplica => renglon fantasma P[i][i]=1, resto=0
// ============================================================
void leerMatrizTransicion(int decision, int m,
                          double P[MAX_ESTADOS][MAX_ESTADOS],
                          int    aplicable[MAX_ESTADOS]) {
    int filas = m + 1;
    printf("\n  === Matriz de transicion para decision k=%d ===\n", decision);

    for (int i = 0; i < filas; i++) {
        int resp;
        printf("  La decision k=%d aplica al estado %d? (1=Si / 0=No): ",
               decision, i);
        scanf("%d", &resp);
        getchar();

        if (resp == 1) {
            aplicable[i] = 1;
            printf("  Ingrese renglon %d (puede usar fracciones a/b):\n", i);
            for (int j = 0; j < filas; j++) {
                printf("    P[%d][%d] = ", i, j);
                P[i][j] = leerValor();
            }
            getchar();
        } else {
            aplicable[i] = 0;
            for (int j = 0; j < filas; j++)
                P[i][j] = (i == j) ? 1.0 : 0.0;
            printf("  [AUTO] Renglon %d: P[%d][%d]=1, resto=0  (renglon fantasma)\n",
                   i, i, i);
        }
    }
}

// ============================================================
//  Construir la matriz de transicion de una politica
// ============================================================
void construirMatrizPolitica(int m,
                             int    politica[MAX_ESTADOS],
                             double Pdec[MAX_DECISIONES][MAX_ESTADOS][MAX_ESTADOS],
                             double Ppol[MAX_ESTADOS][MAX_ESTADOS]) {
    int filas = m + 1;
    for (int i = 0; i < filas; i++) {
        int d = politica[i];
        for (int j = 0; j < filas; j++)
            Ppol[i][j] = Pdec[d][i][j];
    }
}

// ============================================================
//  Resolver pi*P = pi, sum(pi)=1  (Gauss-Jordan)
//  Retorna 1=ok, 0=singular
// ============================================================
int resolverSistemaEstacionario(int m,
                                double P[MAX_ESTADOS][MAX_ESTADOS],
                                double pi[MAX_ESTADOS]) {
    int n = m + 1;
    double A[MAX_ESTADOS][MAX_ESTADOS + 1];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            A[i][j] = P[j][i] - (i == j ? 1.0 : 0.0);
        A[i][n] = 0.0;
    }
    for (int j = 0; j < n; j++) A[n-1][j] = 1.0;
    A[n-1][n] = 1.0;

    for (int col = 0; col < n; col++) {
        int    pivot = -1;
        double best  = 0.0;
        for (int row = col; row < n; row++) {
            if (fabs(A[row][col]) > best) {
                best  = fabs(A[row][col]);
                pivot = row;
            }
        }
        if (pivot == -1 || best < 1e-12) return 0;

        if (pivot != col) {
            for (int j = 0; j <= n; j++) {
                double tmp  = A[col][j];
                A[col][j]   = A[pivot][j];
                A[pivot][j] = tmp;
            }
        }
        for (int row = 0; row < n; row++) {
            if (row == col) continue;
            double factor = A[row][col] / A[col][col];
            for (int j = col; j <= n; j++)
                A[row][j] -= factor * A[col][j];
        }
    }
    for (int i = 0; i < n; i++)
        pi[i] = A[i][n] / A[i][i];

    return 1;
}

// ============================================================
//  Imprimir una matriz n x n con etiqueta
// ============================================================
void imprimirMatriz(const char* titulo, int n,
                    double M[MAX_ESTADOS][MAX_ESTADOS]) {
    printf("  %s\n", titulo);
    printf("       ");
    for (int j = 0; j < n; j++) printf("   j=%-4d", j);
    printf("\n");
    for (int i = 0; i < n; i++) {
        printf("  i=%d [ ", i);
        for (int j = 0; j < n; j++)
            printf("%8.4f ", M[i][j]);
        printf("]\n");
    }
    printf("\n");
}

// ============================================================
//  Imprimir el sistema pi*P = pi y su solucion
//  Muestra la matriz aumentada (P^T - I | 0) con la
//  restriccion sum=1 en el ultimo renglon, y luego el
//  vector pi resultante en formato columna.
// ============================================================
void imprimirSistemaEstacionario(int m,
                                 double P[MAX_ESTADOS][MAX_ESTADOS],
                                 double pi[MAX_ESTADOS]) {
    int n = m + 1;

    printf("  Sistema de ecuaciones estacionarias (pi * P = pi):\n");
    printf("  Ecuaciones (P^T - I) * pi = 0  +  sum(pi) = 1\n");
    printf("  %-30s | RHS\n", "Coeficientes");
    printf("  ");
    for (int j = 0; j < n; j++) printf("  pi_%d    ", j);
    printf("  | RHS\n");
    printf("  ");
    for (int j = 0; j <= n; j++) printf("---------");
    printf("\n");

    // Renglones del sistema (P^T - I)
    for (int i = 0; i < n - 1; i++) {
        printf("  ");
        for (int j = 0; j < n; j++) {
            double coef = P[j][i] - (i == j ? 1.0 : 0.0);
            printf("%8.4f ", coef);
        }
        printf("  |  0\n");
    }
    // Ultimo renglon: sum = 1
    printf("  ");
    for (int j = 0; j < n; j++) printf("%8.4f ", 1.0);
    printf("  |  1\n\n");

    // Solucion
    printf("  Solucion — vector estacionario pi:\n");
    for (int i = 0; i < n; i++)
        printf("    pi_%d = %10.6f\n", i, pi[i]);
    printf("\n");
}

// ============================================================
//  E(C) = sum_i  pi[i] * C[i][ decision_i ]
// ============================================================
double calcularEsperanza(int m,
                         int    politica[MAX_ESTADOS],
                         double pi[MAX_ESTADOS],
                         double costos[MAX_ESTADOS][MAX_DECISIONES]) {
    double ec = 0.0;
    for (int i = 0; i <= m; i++)
        ec += pi[i] * costos[i][ politica[i] ];
    return ec;
}

// ============================================================
//  ALGORITMO 1: Enumeracion Exhaustiva de Politicas
// ============================================================
void algoritmo1() {
    system("cls");
    printf("===== ENUMERACION EXHAUSTIVA DE POLITICAS =====\n");

    // 1. Estados
    int m;
    printf("\nCuantos estados (i) tiene tu problema? (estados van de 0 a m): ");
    scanf("%d", &m);
    getchar();
    int numEstados = m + 1;

    // 2. Decisiones
    int k;
    printf("Cuantas decisiones (k) tiene tu problema?: ");
    scanf("%d", &k);
    getchar();

    printf("\nNumero maximo de politicas (si todas fueran viables): %d^%d = %.0f\n",
           k, numEstados, pow((double)k, (double)numEstados));

    // 3. Matrices de transicion por decision
    double Pdec[MAX_DECISIONES][MAX_ESTADOS][MAX_ESTADOS];
    memset(Pdec, 0, sizeof(Pdec));

    int aplicable[MAX_DECISIONES][MAX_ESTADOS];
    memset(aplicable, 0, sizeof(aplicable));

    printf("\n--- Ingrese las matrices de transicion ---\n");
    printf("    (puede usar fracciones: ej. 3/8)\n");
    for (int d = 1; d <= k; d++)
        leerMatrizTransicion(d, m, Pdec[d], aplicable[d]);

    // 4. Costos
    double costos[MAX_ESTADOS][MAX_DECISIONES];
    memset(costos, 0, sizeof(costos));
    leerCostos(m, k, costos, aplicable);

    // 5. Politicas
    int numPoliticas;
    printf("\nCuantas politicas desea evaluar?: ");
    scanf("%d", &numPoliticas);
    getchar();

    if (numPoliticas > MAX_POLITICAS) {
        printf("Se limita a %d politicas.\n", MAX_POLITICAS);
        numPoliticas = MAX_POLITICAS;
    }

    int politicas[MAX_POLITICAS][MAX_ESTADOS];

    printf("\nIngrese cada politica como vector de %d decisiones (estado 0..%d):\n",
           numEstados, m);

    for (int n = 0; n < numPoliticas; n++) {
        printf("\n  Politica R%d:\n", n + 1);
        for (int i = 0; i < numEstados; i++) {
            printf("    Estado %d -> decision (1..%d): ", i, k);
            scanf("%d", &politicas[n][i]);
        }
        getchar();
    }

    // 6. Minimizar o maximizar
    int objetivo;
    printf("\nDesea minimizar (1) o maximizar (2) el costo esperado?: ");
    scanf("%d", &objetivo);
    getchar();

    // 7. Evaluar cada politica
    double ec[MAX_POLITICAS];
    double Ppol[MAX_ESTADOS][MAX_ESTADOS];
    double pi[MAX_ESTADOS];

    int    politicaOptima = -1;
    double ecOptimo = (objetivo == 1) ? 1e18 : -1e18;

    printf("\n====== RESULTADOS ======\n");

    for (int n = 0; n < numPoliticas; n++) {
        printf("\n-- Politica R%d: (", n + 1);
        for (int i = 0; i < numEstados; i++)
            printf("%d%s", politicas[n][i], i < numEstados - 1 ? "," : "");
        printf(") --\n");

        // Verificar que cada decision aplique al estado correspondiente
        int absurda = 0;
        for (int i = 0; i < numEstados; i++) {
            int d = politicas[n][i];
            if (d < 1 || d > k || !aplicable[d][i]) {
                printf("  [!] Estado %d usa decision %d que NO le aplica => ABSURDA, omitida.\n",
                       i, d);
                absurda = 1;
                break;
            }
        }
        if (absurda) {
            ec[n] = (objetivo == 1) ? 1e18 : -1e18;
            continue;
        }

        // Construir matriz de politica y mostrarla
        construirMatrizPolitica(m, politicas[n], Pdec, Ppol);
        imprimirMatriz("Matriz de transicion de la politica:", numEstados, Ppol);

        // Resolver sistema estacionario
        if (!resolverSistemaEstacionario(m, Ppol, pi)) {
            printf("  Sistema singular => politica omitida.\n");
            ec[n] = (objetivo == 1) ? 1e18 : -1e18;
            continue;
        }

        // Mostrar sistema y vector pi
        imprimirSistemaEstacionario(m, Ppol, pi);

        // Calcular E(C)
        ec[n] = calcularEsperanza(m, politicas[n], pi, costos);
        printf("  E(C) = %.4f\n", ec[n]);

        // Actualizar optimo
        if ((objetivo == 1 && ec[n] < ecOptimo) ||
            (objetivo == 2 && ec[n] > ecOptimo)) {
            ecOptimo       = ec[n];
            politicaOptima = n;
        }
    }

    // 8. Resultado final
    printf("\n====== POLITICA OPTIMA ======\n");
    if (politicaOptima >= 0) {
        printf("R%d: (", politicaOptima + 1);
        for (int i = 0; i < numEstados; i++)
            printf("%d%s", politicas[politicaOptima][i],
                   i < numEstados - 1 ? "," : "");
        printf(")\n");
        printf("E(C) %s = %.4f\n",
               objetivo == 1 ? "minimo" : "maximo", ecOptimo);
    } else {
        printf("No se pudo determinar politica optima (todas absurdas o singulares).\n");
    }
}

// ============================================================
//  ALGORITMO 2: Placeholder Programacion Lineal
// ============================================================
void algoritmo2() {
    system("cls");
    printf("===== PROGRAMACION LINEAL (pendiente) =====\n");
    printf("Esta funcionalidad aun no esta implementada.\n");
}

// ============================================================
//  Main
// ============================================================
int main() {
    Portada();
    menuPrincipal();
    return 0;
}