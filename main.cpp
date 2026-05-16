#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_ESTADOS    20
#define MAX_DECISIONES 20
#define MAX_POLITICAS  10000
#define COSTO_INVALIDO 999999999.0

#define MAX_VARS_PPL   (MAX_ESTADOS * MAX_DECISIONES)
#define MAX_RESTR_PPL  (MAX_ESTADOS + 5)

// ============================================================
//  Variables globales compartidas entre algoritmos
// ============================================================
static int    g_m         = -1;
static int    g_k         = -1;
static int    g_capturado =  0;

static double g_Pdec[MAX_DECISIONES][MAX_ESTADOS][MAX_ESTADOS];
static double g_costos[MAX_ESTADOS][MAX_DECISIONES];
static int    g_aplicable[MAX_DECISIONES][MAX_ESTADOS];

// ============================================================
//  Tableau del Simplex (estatico)
// ============================================================
static double SP_T[MAX_RESTR_PPL][MAX_VARS_PPL * 2 + 1];
static int    SP_base[MAX_RESTR_PPL];

// ============================================================
//  simplex_min: minimiza c^T x  s.a. Ax = b, x >= 0
//  Metodo de dos fases con variables artificiales.
//  Retorna 1=optimo encontrado, 0=infactible/no acotado.
// ============================================================
static int simplex_min(int nv, int nr,
                       double c[],
                       double A[MAX_RESTR_PPL][MAX_VARS_PPL],
                       double b[],
                       double sol[], double *zOpt) {
    int tv = nv + nr;

    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nv; j++)  SP_T[i][j] = A[i][j];
        for (int j = nv; j < tv; j++) SP_T[i][j] = (j - nv == i) ? 1.0 : 0.0;
        SP_T[i][tv] = b[i];
        SP_base[i]  = nv + i;
    }
    for (int j = 0; j < nv; j++) {
        double s = 0;
        for (int i = 0; i < nr; i++) s += A[i][j];
        SP_T[nr][j] = -s;
    }
    for (int j = nv; j < tv; j++) SP_T[nr][j] = 0.0;
    double sb = 0;
    for (int i = 0; i < nr; i++) sb += b[i];
    SP_T[nr][tv] = -sb;

    /* Fase 1 */
    for (int it = 0; it < 30000; it++) {
        int col = -1; double mv = -1e-9;
        for (int j = 0; j < tv; j++) if (SP_T[nr][j] < mv) { mv = SP_T[nr][j]; col = j; }
        if (col < 0) break;
        int row = -1; double mr = 1e18;
        for (int i = 0; i < nr; i++) if (SP_T[i][col] > 1e-9) { double r = SP_T[i][tv] / SP_T[i][col]; if (r < mr) { mr = r; row = i; } }
        if (row < 0) return 0;
        SP_base[row] = col;
        double pv = SP_T[row][col];
        for (int j = 0; j <= tv; j++) SP_T[row][j] /= pv;
        for (int i = 0; i <= nr; i++) { if (i == row) continue; double f = SP_T[i][col]; for (int j = 0; j <= tv; j++) SP_T[i][j] -= f * SP_T[row][j]; }
    }
    if (fabs(SP_T[nr][tv]) > 1e-6) return 0;

    /* Fase 2 */
    for (int j = 0; j <= tv; j++) SP_T[nr][j] = 0.0;
    for (int j = 0; j < nv;  j++) SP_T[nr][j] = c[j];
    SP_T[nr][tv] = 0.0;
    for (int i = 0; i < nr; i++) { int bv = SP_base[i]; if (bv < nv) { double f = SP_T[nr][bv]; for (int j = 0; j <= tv; j++) SP_T[nr][j] -= f * SP_T[i][j]; } }

    for (int it = 0; it < 30000; it++) {
        int col = -1; double mv = -1e-9;
        for (int j = 0; j < nv; j++) if (SP_T[nr][j] < mv) { mv = SP_T[nr][j]; col = j; }
        if (col < 0) break;
        int row = -1; double mr = 1e18;
        for (int i = 0; i < nr; i++) if (SP_T[i][col] > 1e-9) { double r = SP_T[i][tv] / SP_T[i][col]; if (r < mr) { mr = r; row = i; } }
        if (row < 0) return 0;
        SP_base[row] = col;
        double pv = SP_T[row][col];
        for (int j = 0; j <= tv; j++) SP_T[row][j] /= pv;
        for (int i = 0; i <= nr; i++) { if (i == row) continue; double f = SP_T[i][col]; for (int j = 0; j <= tv; j++) SP_T[i][j] -= f * SP_T[row][j]; }
    }

    for (int j = 0; j < nv; j++) sol[j] = 0.0;
    for (int i = 0; i < nr; i++) if (SP_base[i] < nv) sol[SP_base[i]] = SP_T[i][tv];
    *zOpt = SP_T[nr][tv];
    return 1;
}

// ============================================================
//  Prototipos
// ============================================================
void   Portada();
void   menuPrincipal();
void   algoritmo1();
void   algoritmo2();
double leerValor();
void   leerCostos(int m, int k, double costos[MAX_ESTADOS][MAX_DECISIONES], int aplicable[MAX_DECISIONES][MAX_ESTADOS]);
void   leerMatrizTransicion(int decision, int m, double P[MAX_ESTADOS][MAX_ESTADOS], int aplicable[MAX_ESTADOS]);
void   construirMatrizPolitica(int m, int politica[MAX_ESTADOS], double Pdec[MAX_DECISIONES][MAX_ESTADOS][MAX_ESTADOS], double Ppol[MAX_ESTADOS][MAX_ESTADOS]);
int    resolverSistemaEstacionario(int m, double P[MAX_ESTADOS][MAX_ESTADOS], double pi[MAX_ESTADOS]);
double calcularEsperanza(int m, int politica[MAX_ESTADOS], double pi[MAX_ESTADOS], double costos[MAX_ESTADOS][MAX_DECISIONES]);
void   imprimirMatriz(const char* titulo, int n, double M[MAX_ESTADOS][MAX_ESTADOS]);
void   imprimirSistemaEstacionario(int m, double P[MAX_ESTADOS][MAX_ESTADOS], double pi[MAX_ESTADOS]);

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
        if (g_capturado)
            printf("2. Programacion Lineal  [datos cargados: %d estados, %d decisiones]\n",
                   g_m + 1, g_k);
        else
            printf("2. Programacion Lineal  [ejecute opcion 1 primero para cargar datos]\n");
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
//  leerValor: acepta "0.375" o "3/8"
// ============================================================
double leerValor() {
    char buffer[64];
    scanf("%s", buffer);
    char *slash = strchr(buffer, '/');
    if (slash != NULL) {
        *slash = '\0';
        double num = atof(buffer);
        double den = atof(slash + 1);
        if (fabs(den) < 1e-15) { printf("  [!] Denominador cero, se usa 0.\n"); return 0.0; }
        return num / den;
    }
    return atof(buffer);
}

// ============================================================
//  Leer costos
// ============================================================
void leerCostos(int m, int k,
                double costos[MAX_ESTADOS][MAX_DECISIONES],
                int    aplicable[MAX_DECISIONES][MAX_ESTADOS]) {
    printf("\n--- Ingrese los costos C[estado][decision] ---\n");
    printf("    (Los pares no aplicables quedan en %.0f automaticamente)\n", COSTO_INVALIDO);
    for (int i = 0; i <= m; i++) {
        for (int d = 1; d <= k; d++) {
            if (aplicable[d][i]) {
                printf("  C[%d][%d] = ", i, d);
                costos[i][d] = leerValor();
            } else {
                costos[i][d] = COSTO_INVALIDO;
                printf("  C[%d][%d] = %.0f  [no aplicable, asignado automaticamente]\n",
                       i, d, COSTO_INVALIDO);
            }
        }
    }
    getchar();
}

// ============================================================
//  Leer matriz de transicion
// ============================================================
void leerMatrizTransicion(int decision, int m,
                          double P[MAX_ESTADOS][MAX_ESTADOS],
                          int    aplicable[MAX_ESTADOS]) {
    int filas = m + 1;
    printf("\n  === Matriz de transicion para decision k=%d ===\n", decision);
    for (int i = 0; i < filas; i++) {
        int resp;
        printf("  La decision k=%d aplica al estado %d? (1=Si / 0=No): ", decision, i);
        scanf("%d", &resp); getchar();
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
            for (int j = 0; j < filas; j++) P[i][j] = (i == j) ? 1.0 : 0.0;
            printf("  [AUTO] Renglon %d: P[%d][%d]=1, resto=0  (renglon fantasma)\n", i, i, i);
        }
    }
}

// ============================================================
//  Construir la matriz de transicion de una politica
// ============================================================
void construirMatrizPolitica(int m, int politica[MAX_ESTADOS],
                             double Pdec[MAX_DECISIONES][MAX_ESTADOS][MAX_ESTADOS],
                             double Ppol[MAX_ESTADOS][MAX_ESTADOS]) {
    int filas = m + 1;
    for (int i = 0; i < filas; i++) {
        int d = politica[i];
        for (int j = 0; j < filas; j++) Ppol[i][j] = Pdec[d][i][j];
    }
}

// ============================================================
//  Resolver pi*P = pi  (Gauss-Jordan)
// ============================================================
int resolverSistemaEstacionario(int m, double P[MAX_ESTADOS][MAX_ESTADOS], double pi[MAX_ESTADOS]) {
    int n = m + 1;
    double A[MAX_ESTADOS][MAX_ESTADOS + 1];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) A[i][j] = P[j][i] - (i == j ? 1.0 : 0.0);
        A[i][n] = 0.0;
    }
    for (int j = 0; j < n; j++) A[n-1][j] = 1.0;
    A[n-1][n] = 1.0;
    for (int col = 0; col < n; col++) {
        int pivot = -1; double best = 0.0;
        for (int row = col; row < n; row++) if (fabs(A[row][col]) > best) { best = fabs(A[row][col]); pivot = row; }
        if (pivot == -1 || best < 1e-12) return 0;
        if (pivot != col) { for (int j = 0; j <= n; j++) { double tmp = A[col][j]; A[col][j] = A[pivot][j]; A[pivot][j] = tmp; } }
        for (int row = 0; row < n; row++) {
            if (row == col) continue;
            double factor = A[row][col] / A[col][col];
            for (int j = col; j <= n; j++) A[row][j] -= factor * A[col][j];
        }
    }
    for (int i = 0; i < n; i++) pi[i] = A[i][n] / A[i][i];
    return 1;
}

// ============================================================
//  Imprimir matriz
// ============================================================
void imprimirMatriz(const char* titulo, int n, double M[MAX_ESTADOS][MAX_ESTADOS]) {
    printf("  %s\n       ", titulo);
    for (int j = 0; j < n; j++) printf("   j=%-4d", j);
    printf("\n");
    for (int i = 0; i < n; i++) {
        printf("  i=%d [ ", i);
        for (int j = 0; j < n; j++) printf("%8.4f ", M[i][j]);
        printf("]\n");
    }
    printf("\n");
}

// ============================================================
//  Imprimir sistema estacionario
// ============================================================
void imprimirSistemaEstacionario(int m, double P[MAX_ESTADOS][MAX_ESTADOS], double pi[MAX_ESTADOS]) {
    int n = m + 1;
    printf("  Sistema de ecuaciones estacionarias (pi * P = pi):\n");
    printf("  Ecuaciones (P^T - I) * pi = 0  +  sum(pi) = 1\n  ");
    for (int j = 0; j < n; j++) printf("  pi_%d    ", j);
    printf("  | RHS\n  ");
    for (int j = 0; j <= n; j++) printf("---------");
    printf("\n");
    for (int i = 0; i < n - 1; i++) {
        printf("  ");
        for (int j = 0; j < n; j++) printf("%8.4f ", P[j][i] - (i == j ? 1.0 : 0.0));
        printf("  |  0\n");
    }
    printf("  ");
    for (int j = 0; j < n; j++) printf("%8.4f ", 1.0);
    printf("  |  1\n\n  Solucion � vector estacionario pi:\n");
    for (int i = 0; i < n; i++) printf("    pi_%d = %10.6f\n", i, pi[i]);
    printf("\n");
}

// ============================================================
//  E(C)
// ============================================================
double calcularEsperanza(int m, int politica[MAX_ESTADOS],
                         double pi[MAX_ESTADOS], double costos[MAX_ESTADOS][MAX_DECISIONES]) {
    double ec = 0.0;
    for (int i = 0; i <= m; i++) ec += pi[i] * costos[i][ politica[i] ];
    return ec;
}

// ============================================================
//  ALGORITMO 1: Enumeracion Exhaustiva de Politicas
// ============================================================
void algoritmo1() {
    system("cls");
    printf("===== ENUMERACION EXHAUSTIVA DE POLITICAS =====\n");

    printf("\nCuantos estados (i) tiene tu problema? (estados van de 0 a m): ");
    scanf("%d", &g_m); getchar();
    int numEstados = g_m + 1;

    printf("Cuantas decisiones (k) tiene tu problema?: ");
    scanf("%d", &g_k); getchar();

    printf("\nNumero maximo de politicas (si todas fueran viables): %d^%d = %.0f\n",
           g_k, numEstados, pow((double)g_k, (double)numEstados));

    memset(g_Pdec,      0, sizeof(g_Pdec));
    memset(g_aplicable, 0, sizeof(g_aplicable));

    printf("\n--- Ingrese las matrices de transicion ---\n");
    printf("    (puede usar fracciones: ej. 3/8)\n");
    for (int d = 1; d <= g_k; d++)
        leerMatrizTransicion(d, g_m, g_Pdec[d], g_aplicable[d]);

    memset(g_costos, 0, sizeof(g_costos));
    leerCostos(g_m, g_k, g_costos, g_aplicable);

    g_capturado = 1;   /* datos disponibles para algoritmo 2 */

    int numPoliticas;
    printf("\nCuantas politicas desea evaluar?: ");
    scanf("%d", &numPoliticas); getchar();
    if (numPoliticas > MAX_POLITICAS) { printf("Se limita a %d.\n", MAX_POLITICAS); numPoliticas = MAX_POLITICAS; }

    int politicas[MAX_POLITICAS][MAX_ESTADOS];
    printf("\nIngrese cada politica como vector de %d decisiones (estado 0..%d):\n", numEstados, g_m);
    for (int n = 0; n < numPoliticas; n++) {
        printf("\n  Politica R%d:\n", n + 1);
        for (int i = 0; i < numEstados; i++) {
            printf("    Estado %d -> decision (1..%d): ", i, g_k);
            scanf("%d", &politicas[n][i]);
        }
        getchar();
    }

    int objetivo;
    printf("\nDesea minimizar (1) o maximizar (2) el costo esperado?: ");
    scanf("%d", &objetivo); getchar();

    double ec[MAX_POLITICAS];
    double Ppol[MAX_ESTADOS][MAX_ESTADOS], pi[MAX_ESTADOS];
    int    politicaOptima = -1;
    double ecOptimo = (objetivo == 1) ? 1e18 : -1e18;

    printf("\n====== RESULTADOS ======\n");
    for (int n = 0; n < numPoliticas; n++) {
        printf("\n-- Politica R%d: (", n + 1);
        for (int i = 0; i < numEstados; i++) printf("%d%s", politicas[n][i], i < numEstados - 1 ? "," : "");
        printf(") --\n");

        int absurda = 0;
        for (int i = 0; i < numEstados; i++) {
            int d = politicas[n][i];
            if (d < 1 || d > g_k || !g_aplicable[d][i]) {
                printf("  [!] Estado %d usa decision %d que NO le aplica => ABSURDA, omitida.\n", i, d);
                absurda = 1; break;
            }
        }
        if (absurda) { ec[n] = (objetivo == 1) ? 1e18 : -1e18; continue; }

        construirMatrizPolitica(g_m, politicas[n], g_Pdec, Ppol);
        imprimirMatriz("Matriz de transicion de la politica:", numEstados, Ppol);

        if (!resolverSistemaEstacionario(g_m, Ppol, pi)) {
            printf("  Sistema singular => politica omitida.\n");
            ec[n] = (objetivo == 1) ? 1e18 : -1e18; continue;
        }
        imprimirSistemaEstacionario(g_m, Ppol, pi);
        ec[n] = calcularEsperanza(g_m, politicas[n], pi, g_costos);
        printf("  E(C) = %.4f\n", ec[n]);

        if ((objetivo == 1 && ec[n] < ecOptimo) || (objetivo == 2 && ec[n] > ecOptimo)) {
            ecOptimo = ec[n]; politicaOptima = n;
        }
    }

    printf("\n====== POLITICA OPTIMA ======\n");
    if (politicaOptima >= 0) {
        printf("R%d: (", politicaOptima + 1);
        for (int i = 0; i < numEstados; i++) printf("%d%s", politicas[politicaOptima][i], i < numEstados - 1 ? "," : "");
        printf(")\n");
        printf("E(C) %s = %.4f\n", objetivo == 1 ? "minimo" : "maximo", ecOptimo);
    } else {
        printf("No se pudo determinar politica optima (todas absurdas o singulares).\n");
    }
}

// ============================================================
//  ALGORITMO 2: Programacion Lineal
//
//  Reutiliza: g_m, g_k, g_Pdec, g_costos, g_aplicable.
//
//  PPL con variables y_{i,k} solo para pares (i,k) aplicables:
//
//  Min  Z = sum_{i,k} C[i][k] * y_{i,k}
//
//  s.a.
//    (Edo j)  sum_k y_{j,k} - sum_{i,k} P[k][i][j]*y_{i,k} = 0
//             para j = 0 .. m-1
//    (Suma)   sum_{i,k} y_{i,k} = 1
//    y_{i,k} >= 0
//
//  Politica optima: D_{i,k} = y_{i,k} / sum_k y_{i,k}
// ============================================================
void algoritmo2() {
    system("cls");
    printf("===== PROGRAMACION LINEAL =====\n");

    if (!g_capturado) {
        printf("\n[!] No hay datos cargados.\n");
        printf("    Ejecute primero la Enumeracion Exhaustiva (opcion 1)\n");
        printf("    para ingresar las matrices de transicion y costos.\n");
        return;
    }

    int numEstados = g_m + 1;

    /* ---- Enumerar variables aplicables ---- */
    int var_i[MAX_VARS_PPL], var_kv[MAX_VARS_PPL];
    int var_idx[MAX_ESTADOS][MAX_DECISIONES];
    memset(var_idx, -1, sizeof(var_idx));
    int numVars = 0;

    for (int i = 0; i < numEstados; i++)
        for (int d = 1; d <= g_k; d++)
            if (g_aplicable[d][i]) {
                var_i[numVars]  = i;
                var_kv[numVars] = d;
                var_idx[i][d]   = numVars;
                numVars++;
            }

    /* ---- Construir sistema ---- */
    /* numEstados restricciones:
       - filas 0..m-2 : ecuaciones de estado estable para j=0..m-1
       - fila  m-1    : suma total = 1  (reemplaza la ecuacion j=m que es redundante) */
    int numRestr = numEstados;

    static double A_mat[MAX_RESTR_PPL][MAX_VARS_PPL];
    static double b_rhs[MAX_RESTR_PPL];
    static double c_obj[MAX_VARS_PPL];
    memset(A_mat, 0, sizeof(A_mat));
    memset(b_rhs, 0, sizeof(b_rhs));
    memset(c_obj, 0, sizeof(c_obj));

    for (int v = 0; v < numVars; v++)
        c_obj[v] = g_costos[ var_i[v] ][ var_kv[v] ];

    /* Ecuaciones de estado estable j = 0 .. m-1 */
    for (int j = 0; j < numEstados - 1; j++) {
        for (int d = 1; d <= g_k; d++) { int v = var_idx[j][d]; if (v >= 0) A_mat[j][v] += 1.0; }
        for (int v = 0; v < numVars; v++) A_mat[j][v] -= g_Pdec[ var_kv[v] ][ var_i[v] ][j];
        b_rhs[j] = 0.0;
    }
    /* Suma = 1 */
    for (int v = 0; v < numVars; v++) A_mat[numEstados - 1][v] = 1.0;
    b_rhs[numEstados - 1] = 1.0;

    /* ============================================================
       Mostrar el PPL
       ============================================================ */
    printf("\n");
    printf("============================================================\n");
    printf("  FORMULACION DEL PPL\n");
    printf("  y_{i,k} = prob. incondicional estacionaria de estar\n");
    printf("            en estado i y tomar la decision k\n");
    printf("============================================================\n\n");

    printf("  Min Z =");
    int primero = 1;
    for (int v = 0; v < numVars; v++) {
        if (!primero) printf("\n         +");
        printf(" %.4f * y_{%d,%d}", c_obj[v], var_i[v], var_kv[v]);
        primero = 0;
    }
    printf("\n\n  s.a.\n\n");

    for (int j = 0; j < numEstados - 1; j++) {
        printf("  (Edo %d): ", j);
        primero = 1;
        for (int v = 0; v < numVars; v++) {
            double coef = A_mat[j][v];
            if (fabs(coef) < 1e-12) continue;
            if (!primero) {
                if (coef > 0) printf(" + %.4f*y_{%d,%d}", coef, var_i[v], var_kv[v]);
                else          printf(" - %.4f*y_{%d,%d}", -coef, var_i[v], var_kv[v]);
            } else {
                printf("%.4f*y_{%d,%d}", coef, var_i[v], var_kv[v]);
                primero = 0;
            }
        }
        printf(" = 0\n");
    }

    printf("  (Suma)  : ");
    primero = 1;
    for (int v = 0; v < numVars; v++) {
        if (!primero) printf(" + ");
        printf("y_{%d,%d}", var_i[v], var_kv[v]);
        primero = 0;
    }
    printf(" = 1\n\n  y_{i,k} >= 0\n");
    printf("\n============================================================\n");
    printf("\n  Presione Enter para resolver...\n");
    getchar();

    /* ============================================================
       Resolver
       ============================================================ */
    static double sol[MAX_VARS_PPL];
    double zOpt;
    memset(sol, 0, sizeof(sol));

    printf("\n  Resolviendo PPL con Simplex de dos fases...\n");

    int ok = simplex_min(numVars, numRestr, c_obj, A_mat, b_rhs, sol, &zOpt);

    if (!ok) {
        printf("\n  [!] El Simplex no encontro solucion factible.\n");
        printf("      Verifique que las matrices de transicion sean consistentes.\n");
        return;
    }

    /* ============================================================
       Mostrar variables y_{i,k}
       ============================================================ */
    printf("\n  ============================================================\n");
    printf("  SOLUCION OPTIMA  --  Variables y_{i,k}\n");
    printf("  ============================================================\n");
    printf("  %-10s  %-10s  %-16s\n", "Estado i", "Decision k", "y_{i,k}");
    printf("  %-10s  %-10s  %-16s\n", "--------", "----------", "----------------");
    for (int v = 0; v < numVars; v++) {
        double yv = sol[v] < 0.0 ? 0.0 : sol[v];
        printf("  %-10d  %-10d  %-16.8f\n", var_i[v], var_kv[v], yv);
    }
    printf("\n  E(C) optimo = Z = %.6f\n", zOpt);

    /* ============================================================
       Calcular D_{i,k}
       ============================================================ */
    printf("\n  ============================================================\n");
    printf("  POLITICA OPTIMA  --  D_{i,k} = y_{i,k} / sum_k(y_{i,k})\n");
    printf("  ============================================================\n\n");

    printf("  %6s", "Estado");
    for (int d = 1; d <= g_k; d++) printf("    D[i,%d] ", d);
    printf("\n  %6s", "------");
    for (int d = 1; d <= g_k; d++) printf("  --------");
    printf("\n");

    int polOptima[MAX_ESTADOS];
    memset(polOptima, 0, sizeof(polOptima));

    for (int i = 0; i < numEstados; i++) {
        double sumaYi = 0.0;
        for (int d = 1; d <= g_k; d++) { int v = var_idx[i][d]; if (v >= 0 && sol[v] > 0) sumaYi += sol[v]; }

        printf("  %6d", i);
        double maxD = -1.0;
        for (int d = 1; d <= g_k; d++) {
            int v = var_idx[i][d];
            if (v < 0) { printf("  %8s", "N/A"); continue; }
            double yv  = sol[v] < 0 ? 0.0 : sol[v];
            double Dik = (sumaYi > 1e-12) ? yv / sumaYi : 0.0;
            printf("  %8.4f", Dik);
            if (Dik > maxD) { maxD = Dik; polOptima[i] = d; }
        }
        printf("\n");
    }

    printf("\n  Politica determinista resultante R*(");
    for (int i = 0; i < numEstados; i++)
        printf("%d%s", polOptima[i], i < numEstados - 1 ? "," : "");
    printf(")\n");

    /* Verificacion */
    double Ppol[MAX_ESTADOS][MAX_ESTADOS], pi[MAX_ESTADOS];
    construirMatrizPolitica(g_m, polOptima, g_Pdec, Ppol);
    if (resolverSistemaEstacionario(g_m, Ppol, pi)) {
        double ecVerif = calcularEsperanza(g_m, polOptima, pi, g_costos);
        printf("\n  Verificacion con distribucion estacionaria de R*:\n");
        for (int i = 0; i < numEstados; i++) printf("    pi_%d = %.6f\n", i, pi[i]);
        printf("  E(C) verificado = %.6f\n", ecVerif);
    }
    printf("\n  ============================================================\n");
}

// ============================================================
//  Main
// ============================================================
int main() {
    Portada();
    menuPrincipal();
    return 0;
}
