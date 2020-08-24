#ifndef MATRIX8_H
#define MATRIX8_H

struct Point3D{
    double x;
    double y;
    double z;
} Point3D;

//2020

long MBfilas=0;
long MBcolumnas=0;
double *MBmatriz;
long MBanchobanda=0;
double xCg, yCg, zCg;

double porcentaje;

double M8GetProgress();

double M8GetVersion();
double debanda(long ib, long jb);
void abanda(long ib,long jb,double valor);
void addbanda(long ib,long jb,double valor);
void M8borrar(double  *matriz, long filas,long columnas);

void M8matxmatSQ(double  *matriz1, double  *matriz2,double  *resultado, long orden);

long M8gauss(double *matriz, double *x, double * ti, long n , int *salir, void (*cbGambas)(int));

void M8numeroXmatriz(double numero, double *matriz,long filas,long columnas);
void M8simetrizarSQ(double *matriz,long n);
void M8transpuesta( double *m1,double *m2, long f1 , long c1);
long M8VerificarSimetria(double *m , long n);
void M8vectorXmatrizSQ(double *vector, double *matriz, double *resultado, long n);
void M8trasponerSQ(double *matriz, long n );

double M8matrizvalue(double *matriz,long fila,long columna,long orden );
void M8borrar(double  *matriz, long filas, long columnas);

long M8solucionarXgaussBanda(double *matriz, double *x, double * ti, long orden , long anchobanda, int *salir, void (*cbGambas)(int));


double M8matriztest(double *matriz,long fila,long columna, long orden, long anchobanda);
double M8ChequearCondicionamiento(double *matriz, long anchobanda, long orden, void (*cbGambas)(int));

long M8JacobiBanda(double *matriz, long anchobanda, double *x, double *ti, double *x1, long orden, double precision);

// 2020 colgado
void SetMatrizRotacion(double *mr );
void setCG(double x, double y , double z );
void rotar3Dcg( double *punto , double *pRotado );
void SetMatrizRotacionLocal(double *mr);

//Public Extern setCG(x As Float, y As Float, z As Float)
void Local3D( double *punto , double *prRotado);


#endif /* MATRIX8_H */
