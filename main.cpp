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

static int    g_m         = -1;
static int    g_k         = -1;
static int    g_capturado =  0;

static double g_Pdec[MAX_DECISIONES][MAX_ESTADOS][MAX_ESTADOS];
static double g_costos[MAX_ESTADOS][MAX_DECISIONES];
static int    g_aplicable[MAX_DECISIONES][MAX_ESTADOS];

static double SP_T[MAX_RESTR_PPL][MAX_VARS_PPL * 2 + 1];
static int    SP_base[MAX_RESTR_PPL];

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

void   Portada();
void   menuPrincipal();
void   algoritmo1();
void   algoritmo2();
void   algoritmo3();
void   algoritmo4();
double leerValor();
void   leerCostos(int m, int k, double costos[MAX_ESTADOS][MAX_DECISIONES], int aplicable[MAX_DECISIONES][MAX_ESTADOS]);
void   leerMatrizTransicion(int decision, int m, double P[MAX_ESTADOS][MAX_ESTADOS], int aplicable[MAX_ESTADOS]);
void   construirMatrizPolitica(int m, int politica[MAX_ESTADOS], double Pdec[MAX_DECISIONES][MAX_ESTADOS][MAX_ESTADOS], double Ppol[MAX_ESTADOS][MAX_ESTADOS]);
int    resolverSistemaEstacionario(int m, double P[MAX_ESTADOS][MAX_ESTADOS], double pi[MAX_ESTADOS]);
double calcularEsperanza(int m, int politica[MAX_ESTADOS], double pi[MAX_ESTADOS], double costos[MAX_ESTADOS][MAX_DECISIONES]);
void   imprimirMatriz(const char* titulo, int n, double M[MAX_ESTADOS][MAX_ESTADOS]);
void   imprimirSistemaEstacionario(int m, double P[MAX_ESTADOS][MAX_ESTADOS], double pi[MAX_ESTADOS]);

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
        if (g_capturado)
            printf("3. Mejoramiento de Politicas  [datos cargados: %d estados, %d decisiones]\n",
                   g_m + 1, g_k);
        else
            printf("3. Mejoramiento de Politicas  [ejecute opcion 1 primero para cargar datos]\n");
        if (g_capturado)
            printf("4. Mejoramiento de Politicas con Descuento  [datos cargados: %d estados, %d decisiones]\n",
                   g_m + 1, g_k);
        else
            printf("4. Mejoramiento de Politicas con Descuento  [ejecute opcion 1 primero para cargar datos]\n");
        printf("5. Salir\n");
        printf("=====================================\n");
        printf("Selecciona una opcion: ");
        scanf("%d", &x);
        getchar();

        switch (x) {
            case 1: algoritmo1(); break;
            case 2: algoritmo2(); break;
            case 3: algoritmo3(); break;
            case 4: algoritmo4(); break;
            case 5:
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

void construirMatrizPolitica(int m, int politica[MAX_ESTADOS],
                             double Pdec[MAX_DECISIONES][MAX_ESTADOS][MAX_ESTADOS],
                             double Ppol[MAX_ESTADOS][MAX_ESTADOS]) {
    int filas = m + 1;
    for (int i = 0; i < filas; i++) {
        int d = politica[i];
        for (int j = 0; j < filas; j++) Ppol[i][j] = Pdec[d][i][j];
    }
}

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
    printf("  |  1\n\n  Solucion => vector estacionario pi:\n");
    for (int i = 0; i < n; i++) printf("    pi_%d = %10.6f\n", i, pi[i]);
    printf("\n");
}

double calcularEsperanza(int m, int politica[MAX_ESTADOS],
                         double pi[MAX_ESTADOS], double costos[MAX_ESTADOS][MAX_DECISIONES]) {
    double ec = 0.0;
    for (int i = 0; i <= m; i++) ec += pi[i] * costos[i][ politica[i] ];
    return ec;
}

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

    g_capturado = 1;

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

    int numRestr = numEstados;

    static double A_mat[MAX_RESTR_PPL][MAX_VARS_PPL];
    static double b_rhs[MAX_RESTR_PPL];
    static double c_obj[MAX_VARS_PPL];
    memset(A_mat, 0, sizeof(A_mat));
    memset(b_rhs, 0, sizeof(b_rhs));
    memset(c_obj, 0, sizeof(c_obj));

    for (int v = 0; v < numVars; v++)
        c_obj[v] = g_costos[ var_i[v] ][ var_kv[v] ];

    for (int j = 0; j < numEstados - 1; j++) {
        for (int d = 1; d <= g_k; d++) { int v = var_idx[j][d]; if (v >= 0) A_mat[j][v] += 1.0; }
        for (int v = 0; v < numVars; v++) A_mat[j][v] -= g_Pdec[ var_kv[v] ][ var_i[v] ][j];
        b_rhs[j] = 0.0;
    }
    for (int v = 0; v < numVars; v++) A_mat[numEstados - 1][v] = 1.0;
    b_rhs[numEstados - 1] = 1.0;

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

static long long mcd(long long a, long long b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { long long t = b; b = a % b; a = t; }
    return a;
}

static void imprimirFraccion(double val, int forzar_decimal) {
    if (forzar_decimal) { printf("%.6f", val); return; }
    long long best_num = 0, best_den = 1;
    double best_err = fabs(val);
    for (long long den = 1; den <= 10000; den++) {
        long long num = (long long)round(val * den);
        double err = fabs(val - (double)num / den);
        if (err < best_err) { best_err = err; best_num = num; best_den = den; }
        if (best_err < 1e-9) break;
    }
    if (best_err < 1e-9 && best_den > 1) {
        long long g = mcd(best_num < 0 ? -best_num : best_num, best_den);
        printf("%lld/%lld (%.6f)", best_num / g, best_den / g, val);
    } else {
        printf("%.6f", val);
    }
}

static int resolverPaso1(int m, int politica[MAX_ESTADOS],
                         double Pdec[MAX_DECISIONES][MAX_ESTADOS][MAX_ESTADOS],
                         double costos[MAX_ESTADOS][MAX_DECISIONES],
                         double *gRn, double V[MAX_ESTADOS]) {
    int n = m + 1;
    static double A[MAX_ESTADOS + 1][MAX_ESTADOS + 2];
    memset(A, 0, sizeof(A));

    for (int i = 0; i < n; i++) {
        int k = politica[i];
        A[i][0] = 1.0;
        for (int j = 0; j < m; j++)
            A[i][1 + j] = (i == j ? 1.0 : 0.0) - Pdec[k][i][j];
        A[i][n] = costos[i][k];
    }

    for (int col = 0; col < n; col++) {
        int pivot = -1; double best = 0.0;
        for (int row = col; row < n; row++)
            if (fabs(A[row][col]) > best) { best = fabs(A[row][col]); pivot = row; }
        if (pivot < 0 || best < 1e-12) return 0;
        if (pivot != col)
            for (int j = 0; j <= n; j++) {
                double t = A[col][j]; A[col][j] = A[pivot][j]; A[pivot][j] = t;
            }
        double pv = A[col][col];
        for (int j = 0; j <= n; j++) A[col][j] /= pv;
        for (int row = 0; row < n; row++) {
            if (row == col) continue;
            double f = A[row][col];
            for (int j = 0; j <= n; j++) A[row][j] -= f * A[col][j];
        }
    }

    *gRn = A[0][n];
    for (int i = 0; i < m; i++) V[i] = A[1 + i][n];
    V[m] = 0.0;
    return 1;
}

void algoritmo3() {
    system("cls");
    printf("===== MEJORAMIENTO DE POLITICAS =====\n");

    if (!g_capturado) {
        printf("\n[!] No hay datos cargados.\n");
        printf("    Ejecute primero la Enumeracion Exhaustiva (opcion 1)\n");
        printf("    para ingresar las matrices de transicion y costos.\n");
        return;
    }

    int numEstados = g_m + 1;

    printf("\n========================================\n");
    printf("  PASO 0: Politica de arranque R1\n");
    printf("========================================\n");
    printf("  Ingrese la politica inicial R1 (vector de %d decisiones):\n", numEstados);

    int politica[MAX_ESTADOS];
    for (int i = 0; i < numEstados; i++) {
        int d;
        do {
            printf("    Estado %d -> decision (1..%d): ", i, g_k);
            scanf("%d", &d); getchar();
            if (d < 1 || d > g_k || !g_aplicable[d][i])
                printf("    [!] Decision %d no valida para estado %d. Intente de nuevo.\n", d, i);
        } while (d < 1 || d > g_k || !g_aplicable[d][i]);
        politica[i] = d;
    }

    int iteracion = 1;

    while (1) {
        printf("\n========================================\n");
        printf("  ITERACION %d\n", iteracion);
        printf("========================================\n");
        printf("\n  --- PASO 1: Determinacion del Valor ---\n");
        printf("  Politica R%d: (", iteracion);
        for (int i = 0; i < numEstados; i++)
            printf("%d%s", politica[i], i < numEstados - 1 ? "," : "");
        printf(")\n\n");

        printf("  Sistema (con V_%d = 0):\n\n", g_m);
        for (int i = 0; i < numEstados; i++) {
            int k = politica[i];
            printf("    i=%d:  g(R%d)", i, iteracion);
            for (int j = 0; j < g_m; j++) {
                double coef = (i == j ? 1.0 : 0.0) - g_Pdec[k][i][j];
                if (fabs(coef) < 1e-12) continue;
                if (coef > 0) printf(" + ");
                else          printf(" - ");
                imprimirFraccion(fabs(coef), 0);
                printf("*V_%d", j);
            }
            {
                double coef_m = (i == g_m ? 1.0 : 0.0) - g_Pdec[k][i][g_m];
                if (fabs(coef_m) > 1e-12)
                    printf(" [+ (");
                    imprimirFraccion(coef_m, 0);
                    printf(")*0]");
            }
            printf(" = ");
            imprimirFraccion(g_costos[i][k], 0);
            printf("\n");
        }

        double gRn;
        double V[MAX_ESTADOS];
        if (!resolverPaso1(g_m, politica, g_Pdec, g_costos, &gRn, V)) {
            printf("\n  [!] Sistema singular. No se puede continuar.\n");
            break;
        }

        printf("\n  Resultados:\n");
        printf("    g(R%d) = ", iteracion);
        imprimirFraccion(gRn, 0);
        printf("\n");
        for (int i = 0; i < numEstados; i++) {
            printf("    V_%d   = ", i);
            imprimirFraccion(V[i], 0);
            printf("\n");
        }

        printf("\n  --- PASO 2: Mejoramiento de la Politica ---\n");
        printf("  Para cada estado i, minimizar:  C_{ik} + sum_j P_{ij}(k)*V_j - V_i\n\n");

        int nuevaPolitica[MAX_ESTADOS];
        int cambio = 0;

        for (int i = 0; i < numEstados; i++) {
            printf("  Estado i=%d:\n", i);
            double mejorVal = 1e18;
            int    mejorK   = -1;

            for (int k = 1; k <= g_k; k++) {
                if (!g_aplicable[k][i]) continue;

                double val = g_costos[i][k];
                for (int j = 0; j < numEstados; j++)
                    val += g_Pdec[k][i][j] * V[j];
                val -= V[i];

                printf("    k=%d:  C[%d][%d] + sum_j P_{%d,j}(k)*V_j - V_%d = ",
                       k, i, k, i, i);
                imprimirFraccion(val, 0);

                if (val < mejorVal - 1e-9) {
                    mejorVal = val;
                    mejorK   = k;
                    printf("  <-- minimo hasta ahora");
                }
                printf("\n");
            }
            nuevaPolitica[i] = mejorK;
            printf("    => d_%d(R%d) = k=%d\n", i, iteracion + 1, mejorK);
        }

        printf("\n  Nueva politica R%d: (", iteracion + 1);
        for (int i = 0; i < numEstados; i++)
            printf("%d%s", nuevaPolitica[i], i < numEstados - 1 ? "," : "");
        printf(")\n");

        for (int i = 0; i < numEstados; i++)
            cambio += (nuevaPolitica[i] != politica[i]);

        if (cambio == 0) {
            printf("\n  ========================================\n");
            printf("  PRUEBA DE OPTIMALIDAD: R%d = R%d\n", iteracion + 1, iteracion);
            printf("  => Se detiene el proceso.\n");
            printf("  => R%d es la POLITICA OPTIMA.\n", iteracion);
            printf("  ========================================\n");
            printf("\n  RESUMEN FINAL:\n");
            printf("    Politica optima: R*(");
            for (int i = 0; i < numEstados; i++)
                printf("%d%s", politica[i], i < numEstados - 1 ? "," : "");
            printf(")\n");
            printf("    g(R*) = costo promedio a largo plazo = ");
            imprimirFraccion(gRn, 0);
            printf("\n");
            for (int i = 0; i < numEstados; i++) {
                printf("    V_%d = ", i);
                imprimirFraccion(V[i], 0);
                printf("\n");
            }
            break;
        }

        for (int i = 0; i < numEstados; i++)
            politica[i] = nuevaPolitica[i];
        iteracion++;

        printf("\n  [La politica cambio, continuando iteracion %d...]\n", iteracion);
        printf("  Presione Enter para continuar...\n");
        getchar();
    }
}

/* ============================================================
   ALGORITMO 4: Mejoramiento de Politicas con Descuento
   (Criterio del valor descontado - seccion 4.3.2)
   ============================================================ */

static int resolverSistemaDescuento(int m, int politica[MAX_ESTADOS],
                                    double Pdec[MAX_DECISIONES][MAX_ESTADOS][MAX_ESTADOS],
                                    double costos[MAX_ESTADOS][MAX_DECISIONES],
                                    double alpha,
                                    double V[MAX_ESTADOS]) {
    int n = m + 1;
    static double A[MAX_ESTADOS][MAX_ESTADOS + 1];
    memset(A, 0, sizeof(A));

    for (int i = 0; i < n; i++) {
        int k = politica[i];
        for (int j = 0; j < n; j++)
            A[i][j] = (i == j ? 1.0 : 0.0) - alpha * Pdec[k][i][j];
        A[i][n] = costos[i][k];
    }

    for (int col = 0; col < n; col++) {
        int pivot = -1; double best = 0.0;
        for (int row = col; row < n; row++)
            if (fabs(A[row][col]) > best) { best = fabs(A[row][col]); pivot = row; }
        if (pivot < 0 || best < 1e-12) return 0;
        if (pivot != col)
            for (int j = 0; j <= n; j++) {
                double t = A[col][j]; A[col][j] = A[pivot][j]; A[pivot][j] = t;
            }
        double pv = A[col][col];
        for (int j = 0; j <= n; j++) A[col][j] /= pv;
        for (int row = 0; row < n; row++) {
            if (row == col) continue;
            double f = A[row][col];
            for (int j = 0; j <= n; j++) A[row][j] -= f * A[col][j];
        }
    }

    for (int i = 0; i < n; i++) V[i] = A[i][n];
    return 1;
}

void algoritmo4() {
    system("cls");
    printf("===== MEJORAMIENTO DE POLITICAS CON DESCUENTO =====\n");
    printf("       (Criterio del Valor Descontado - Sec. 4.3.2)\n");

    if (!g_capturado) {
        printf("\n[!] No hay datos cargados.\n");
        printf("    Ejecute primero la Enumeracion Exhaustiva (opcion 1)\n");
        printf("    para ingresar las matrices de transicion y costos.\n");
        return;
    }

    int numEstados = g_m + 1;

    printf("\n============================================================\n");
    printf("  FACTOR DE DESCUENTO\n");
    printf("  alpha = 1 / (1 + epsilon)   donde epsilon = tasa de interes\n");
    printf("  Ejemplo: si epsilon = 0.10 (10%%)  =>  alpha = 1/1.10 = 0.9091\n");
    printf("============================================================\n");

    double alpha;
    int modo;
    printf("\nDesea ingresar (1) epsilon (tasa de interes) o (2) alpha directamente?: ");
    scanf("%d", &modo); getchar();

    if (modo == 1) {
        double epsilon;
        printf("Ingrese epsilon (tasa de interes, ej. 0.10 para 10%%): ");
        epsilon = leerValor();
        if (epsilon <= -1.0) {
            printf("[!] Tasa de interes invalida (debe ser > -1). Se usara epsilon=0.\n");
            epsilon = 0.0;
        }
        alpha = 1.0 / (1.0 + epsilon);
        printf("=> alpha = 1 / (1 + %.6f) = %.8f\n", epsilon, alpha);
    } else {
        printf("Ingrese alpha directamente (0 < alpha < 1, ej. 9/10): ");
        alpha = leerValor();
        if (alpha <= 0.0 || alpha >= 1.0) {
            printf("[!] alpha debe estar en (0,1). Se usara alpha=0.9.\n");
            alpha = 0.9;
        }
        printf("=> alpha = %.8f\n", alpha);
    }

    /* ---- Paso 0: Politica inicial ---- */
    printf("\n========================================\n");
    printf("  PASO 0: Politica de arranque R1\n");
    printf("========================================\n");
    printf("  Ingrese la politica inicial R1 (vector de %d decisiones):\n", numEstados);

    int politica[MAX_ESTADOS];
    for (int i = 0; i < numEstados; i++) {
        int d;
        do {
            printf("    Estado %d -> decision (1..%d): ", i, g_k);
            scanf("%d", &d); getchar();
            if (d < 1 || d > g_k || !g_aplicable[d][i])
                printf("    [!] Decision %d no valida para estado %d. Intente de nuevo.\n", d, i);
        } while (d < 1 || d > g_k || !g_aplicable[d][i]);
        politica[i] = d;
    }

    int iteracion = 1;

    while (1) {
        /* ---- Paso 1: Determinacion del valor con descuento ---- */
        printf("\n========================================\n");
        printf("  ITERACION %d\n", iteracion);
        printf("========================================\n");
        printf("\n  --- PASO 1: Determinacion del Valor (con descuento alpha=");
        imprimirFraccion(alpha, 0);
        printf(") ---\n");
        printf("  Politica R%d: (", iteracion);
        for (int i = 0; i < numEstados; i++)
            printf("%d%s", politica[i], i < numEstados - 1 ? "," : "");
        printf(")\n\n");

        /* Mostrar ecuacion de Bellman expandida por terminos */
        printf("  Ecuacion de Bellman con descuento:\n");
        printf("  V_i(R%d) = C_ik + alpha * { sum_j P_ij(k) * V_j(R%d) }\n\n", iteracion, iteracion);
        printf("  Sistema lineal  (I - alpha*P(R%d)) * V = C(R%d):\n\n", iteracion, iteracion);

        for (int i = 0; i < numEstados; i++) {
            int k = politica[i];

            /* --- Mostrar la ecuacion de Bellman expandida (forma original) --- */
            printf("    i=%d [k=%d]:  V_%d = ", i, k, i);
            imprimirFraccion(g_costos[i][k], 0);
            printf(" + ");
            imprimirFraccion(alpha, 0);
            printf(" * { ");
            for (int j = 0; j < numEstados; j++) {
                if (j > 0) printf(" + ");
                printf("P_%d%d(%d)*V_%d", i, j, k, j);
            }
            printf(" }\n");

            /* Sustituir los valores de P conocidos */
            printf("          V_%d = ", i);
            imprimirFraccion(g_costos[i][k], 0);
            printf(" + ");
            imprimirFraccion(alpha, 0);
            printf(" * { ");
            int primero_t = 1;
            for (int j = 0; j < numEstados; j++) {
                double pij = g_Pdec[k][i][j];
                if (fabs(pij) < 1e-12) continue;
                if (!primero_t) printf(" + ");
                imprimirFraccion(pij, 0);
                printf("*V_%d", j);
                primero_t = 0;
            }
            printf(" }\n");

            /* Mostrar la ecuacion reordenada (forma del sistema lineal) */
            printf("          ");
            int primero_eq = 1;
            for (int j = 0; j < numEstados; j++) {
                double coef = (i == j ? 1.0 : 0.0) - alpha * g_Pdec[k][i][j];
                if (fabs(coef) < 1e-12) continue;
                if (!primero_eq) {
                    if (coef > 0) printf(" + ");
                    else          printf(" - ");
                    imprimirFraccion(fabs(coef), 0);
                } else {
                    if (coef < 0) { printf("-"); imprimirFraccion(fabs(coef), 0); }
                    else            imprimirFraccion(coef, 0);
                    primero_eq = 0;
                }
                printf("*V_%d", j);
            }
            printf("  =  ");
            imprimirFraccion(g_costos[i][k], 0);
            printf("   [Ec. %d]\n\n", i + 1);
        }

        /* Resolver el sistema */
        double V[MAX_ESTADOS];
        if (!resolverSistemaDescuento(g_m, politica, g_Pdec, g_costos, alpha, V)) {
            printf("\n  [!] Sistema singular. No se puede continuar.\n");
            break;
        }

        printf("  Valores V_i(R%d) resueltos:\n", iteracion);
        for (int i = 0; i < numEstados; i++) {
            printf("    V_%d = ", i);
            imprimirFraccion(V[i], 0);
            printf("\n");
        }

        /* ---- Paso 2: Mejoramiento de la politica ---- */
        printf("\n  --- PASO 2: Mejoramiento de la Politica ---\n");
        printf("  Para cada estado i, encontrar k que minimice:\n");
        printf("      C_ik + alpha * sum_j P_ij(k) * V_j(R%d)\n\n", iteracion);

        int nuevaPolitica[MAX_ESTADOS];
        int cambio = 0;

        for (int i = 0; i < numEstados; i++) {
            /* Contar decisiones aplicables */
            int num_aplic = 0;
            for (int k = 1; k <= g_k; k++) if (g_aplicable[k][i]) num_aplic++;

            if (num_aplic == 1) {
                /* Decision unica: mostrar como en el PDF */
                for (int k = 1; k <= g_k; k++) {
                    if (!g_aplicable[k][i]) continue;
                    double val = g_costos[i][k];
                    for (int j = 0; j < numEstados; j++)
                        val += alpha * g_Pdec[k][i][j] * V[j];
                    printf("  i=%d  Decision Unica\n", i);
                    printf("    k=%d:  V_%d = ", k, i);
                    imprimirFraccion(val, 0);
                    printf("  * MINIMO  =>  d_%d(R%d) = %d\n", i, iteracion + 1, k);
                    nuevaPolitica[i] = k;
                }
            } else {
                printf("  i=%d:\n", i);
                double mejorVal = 1e18;
                int    mejorK   = -1;

                for (int k = 1; k <= g_k; k++) {
                    if (!g_aplicable[k][i]) continue;

                    /* Mostrar la forma expandida: C_ik + alpha*{suma} */
                    printf("    k=%d:  V_%d = C_%d%d + ", k, i, i, k);
                    imprimirFraccion(alpha, 0);
                    printf(" * { ");
                    int pt = 1;
                    for (int j = 0; j < numEstados; j++) {
                        double pij = g_Pdec[k][i][j];
                        if (fabs(pij) < 1e-12) continue;
                        if (!pt) printf(" + ");
                        printf("P_%d%d(%d)(", i, j, k);
                        imprimirFraccion(pij, 0);
                        printf(")*V_%d(", j);
                        imprimirFraccion(V[j], 0);
                        printf(")");
                        pt = 0;
                    }
                    printf(" }\n");

                    /* Calcular el valor numerico */
                    double val = g_costos[i][k];
                    for (int j = 0; j < numEstados; j++)
                        val += alpha * g_Pdec[k][i][j] * V[j];

                    printf("          V_%d = ", i);
                    imprimirFraccion(g_costos[i][k], 0);
                    printf(" + ");
                    imprimirFraccion(val - g_costos[i][k], 0);
                    printf(" = ");
                    imprimirFraccion(val, 0);

                    if (val < mejorVal - 1e-9) {
                        mejorVal = val;
                        mejorK   = k;
                        printf("  * MINIMO");
                    }
                    printf("\n");
                }
                nuevaPolitica[i] = mejorK;
                printf("    => d_%d(R%d) = k=%d\n", i, iteracion + 1, mejorK);
            }
            printf("\n");
        }

        printf("  Nueva politica R%d: (", iteracion + 1);
        for (int i = 0; i < numEstados; i++)
            printf("%d%s", nuevaPolitica[i], i < numEstados - 1 ? "," : "");
        printf(")\n");

        /* ---- Prueba de optimalidad ---- */
        for (int i = 0; i < numEstados; i++)
            cambio += (nuevaPolitica[i] != politica[i]);

        if (cambio == 0) {
            printf("\n  ========================================\n");
            printf("  PRUEBA DE OPTIMALIDAD: R%d = R%d\n", iteracion + 1, iteracion);
            printf("  => R%d == R%d  =>  Se encontro la politica optima.\n",
                   iteracion + 1, iteracion);
            printf("  ========================================\n");
            printf("\n  RESUMEN FINAL:\n");
            printf("    alpha  = ");
            imprimirFraccion(alpha, 0);
            printf("\n");
            printf("    Politica optima: R*(");
            for (int i = 0; i < numEstados; i++)
                printf("%d%s", politica[i], i < numEstados - 1 ? "," : "");
            printf(")\n\n");
            printf("    Valores descontados optimos V_i(R*):\n");
            for (int i = 0; i < numEstados; i++) {
                printf("      V_%d = ", i);
                imprimirFraccion(V[i], 0);
                printf("\n");
            }
            printf("\n  Interpretacion: V_i representa el valor presente total\n");
            printf("  del costo esperado partiendo del estado i bajo la politica optima.\n");
            break;
        }

        for (int i = 0; i < numEstados; i++)
            politica[i] = nuevaPolitica[i];
        iteracion++;

        printf("\n  [La politica cambio, continuando con iteracion %d...]\n", iteracion);
        printf("  Presione Enter para continuar...\n");
        getchar();
    }
}

int main() {
    Portada();
    menuPrincipal();
    return 0;
}