#include "matrix8.h"
#include <math.h>


double porcentaje=0;
double MRotacion[9];
double MRotacionLocal[9];



double M8GetProgress()
{
	//return porcentaje;
	return 0;
}


//BOOL DoEvents();




double M8GetVersion()
{double v = 0.04;
	return v;
}

void M8matxmatSQ(double  *matriz1, double  *matriz2, double  *resultado, long orden)
{


    // rutina valida  para matrices cuadradas

    long i,j,k,c1,c3;

    i=0;
    j=0;
    for(i=0;i<orden;i++)
    {
        c1=i*orden;
        for(j=0;j<orden;j++)
            resultado[j+c1]=0;//matriz1[j+c1];
    }

    //return 0;

    i=0;
    j=0;
    k=0;
    //s=0;
    for(i=0;i<orden;i++)
    {
        c1=i*orden;
        for(j=0;j<orden;j++)
        {

            for(k=0;k<orden;k++)
            {
                c3=k*orden;
                resultado[j+c1] += matriz1[c1+k] * matriz2[c3+j];

            }
        }
    }

    return ;
}


long M8gauss(double *matriz, double *x, double * ti, long n , int *salir,   void (*cbGambas) (int) )
{

    if(n == 0)return 1;
    int i, i1, j, h, h1;
    //long ok = -1;
    //long contador = 0;
    double c1;
    // triangulacion

    for (i = 0;i<n;i++)
    {
        i1=i*n;
        //*porcentaje = i/n;
        if(*salir == -1)return i;

        // llamo a la rutina de gambas que informa al user de mi progreso y mantiene vivo el FMain
        cbGambas(i);

        if(matriz[i+i1] != 0)
        {
            for (j = (i + 1);j<n;j++)
            {
                if(matriz[j+i1] != 0)
                {
                    c1 = -matriz[j+i1] / matriz[i+i1];
                    for(h = i; h<n ; h++)
                    {
                        h1=h*n;
                        matriz[j+h1] += c1 * matriz[i+h1];
                        //DoEvents //necesario para el ambiente
                    }
                    ti[j] = ti[j] + c1 * ti[i];
                }
            }
        }
         else
            // la posisicon de la diagonal ppal donde es = 0
             return i+i1;
    }

    // sustitucion atrás
    double  p1;

    x[n-1] = ti[n-1] / matriz[n* n-1];

    for (i = n - 2;i>=0;i--)
    {
        p1 = 0;
        for (j = i + 1;j<n;j++)
            p1 += - x[j] * matriz[i+ j*n];

        p1 += ti[i];
        x[i] = p1 / matriz[i+i*n];
    }
    return -1;

}


void M8trasponerSQ(double *matriz, long n )
{
	// Gambas
    double ww=0;
    long i,j;

    for (i = 0;i<n;i++)
        for (j = i ;j<n;j++)
        {
            ww = matriz[i*n+ j];
            matriz[i*n+ j] = matriz[j*n+ i];
            matriz[j*n+ i] = ww;
        }

}

void  M8vectorXmatrizSQ(double *vector, double *matriz, double *resultado, long n)
{
    long i,j;
        for (i = 0;i<n;i++)
            resultado[i] = 0;

        for (i = 0; i<n;i++)
            for (j = 0;j<n;j++)
                resultado[i] = resultado[i] + vector[j] *matriz[j+ i*n];



    }

long  M8VerificarSimetria(double *m , long n)
{	// Gambas
    long fila =0, columna =0;

    for (fila = 0;fila<n-1;fila++)
        for (columna = fila + 1;columna<n;columna++)
            if ( m[fila*n+ columna] != m[columna*n+ fila])
                return 0;

    return-1;

}


void M8transpuesta( double * m1,double *  m2, long f1 , long c1)
// GAMBAS
// CREO LA MATRIZ TRANSPUESTA
// VA A TENER CxF
{
    long f,c;
    for (f = 0;f<f1;f++)
        for( c = 0;c<c1;c++)
            m2[c*c1+ f] = m1[f*f1+ c];

        // verificar los resultados
}

void M8simetrizarSQ(double * matriz,long n)
// GAMBAS
// realiza una simetria de la matriz tomando como que
// los valores sobre y en la diag. son los que están.
{
    long a,b;
    for (a = 0;a<n;a++)
        for (b = a+1;b<n;b++)
            matriz[b*n+ a] = matriz[a*n+ b];
}


void M8numeroXmatriz(double numero, double *matriz,long filas,long columnas)
{	// Gambas
    long a,b;
    for (a = 0; a<filas;a++)
        for (b = 0;b<columnas;b++)
            matriz[a*filas+ b] = matriz[a*filas+ b] * numero;
}


void abanda(long ib,long jb,double valor)
{
    long col;

    col = jb - ib ;

    if ((col >= 0) && (col <= MBanchobanda))
        MBmatriz [ib*MBanchobanda+col] = valor;

    return;
}

void addbanda(long ib,long jb,double valor)
{

    long col;

    col = jb - ib ;

if ((col >= 0) && (col <= MBanchobanda))
        MBmatriz [ib*MBanchobanda+col] += valor;
    return;
}


double debanda(long ib, long jb)
// ojo, llama de 0 to orden-1,0 to orden-1
{
    long col ,fila;

    if (ib <= jb)
    {
        fila =ib;
        col = jb - ib  ;
    }
    else
    {
        fila = jb;
        col = ib - jb ;
    }

    if ((col >= 0) && (col < MBanchobanda))
        return MBmatriz [fila*MBanchobanda+col];

    else
        return 0;
}



void M8borrar(double  *matriz, long filas,long columnas)
{
    // gambas
    long i,j;
    for(i=0;i<filas;i++)
        for(j=0;j<columnas;j++)
            matriz[i*columnas+j]=0;

    return ;
}

long M8solucionarXgaussBanda(double *matriz, double * x, double * ti, long orden , long anchobanda, int *salir, void (*cbGambas) (int))
    //
{
    double aaa,bbb,d,e,p1;
    long j,i,l,t1,paso_qc = 0 ,s,totalqc,n;
    int progreso;
    n=orden;
    s = anchobanda;
    totalqc = orden*2 ;

    // seteo vars globales

    MBanchobanda=anchobanda;
    MBfilas=orden;
    s = anchobanda;
    MBmatriz=matriz;
    progreso=0;

    // triangulacion
    for (i = 0; i< n - 1;i++)
    {//1
        paso_qc += 1;

        progreso=paso_qc*100/totalqc;

        // llamo a la rutina de gambas que informa al user de mi progreso y mantiene vivo el FMain
        cbGambas(progreso);

        if(*salir == -1)return i;

        aaa = debanda(i, i);
        if (aaa != 0)
        {//2
            if (i <= n-1 - s + 1)
            {//3
                t1 = i + s - 1;
            }//3
            else
            {//3
                t1 = n-1;
            }//3

            for (j = i + 1;j<=t1;j++)
            {//3
                bbb = debanda(i, j);
                if (bbb != 0 )
                {//4
                    for (l = j;l<=t1;l++)
                    {//5
                        d = debanda(i, l);
                        e = debanda(j, l);
                        if (d != 0)
                            abanda (j, l, e - d * bbb / aaa);
                    }//5
                    ti[j] += - bbb / aaa * ti[i];
                }//4
            }//3
        }//2
        else // hay un 0 en la diagonal ppal
            return -3;

    }//1

    // sustitucion atrás
    x[n-1] = ti[n-1] / debanda(n-1, n-1);
    for (i = n - 2;i>=0 ;i--)
        //DoEvents //necesario para el ambiente
    {
        paso_qc += 1;
        progreso =paso_qc*100/totalqc;
        // llamo a la rutina de gambas que informa al user de mi progreso y mantiene vivo el FMain
        cbGambas(progreso);

        p1 = 0;
        for (j = n-1;j>=i + 1;j--)
            p1 += x[j] * debanda(i, j);

        x[i] = (ti[i] - p1) / debanda(i, i);
    }

    return -1;
}

double M8matriztest(double *matriz,long fila,long columna,long orden ,long anchobanda)
// test ok , devuelve el valor de la matriz banda
{

    MBmatriz=matriz;
    MBfilas=orden;
    MBanchobanda=anchobanda;
    return debanda(fila,columna);
}
double M8matrizvalue(double *matriz,long fila,long columna,long orden )
// test ok , devuelve el valor de la matriz banda
{

    return matriz[fila*orden+columna];
}

double M8ChequearCondicionamiento(double *matriz,long anchobanda, long orden, void (*cbGambas) (int))
//
// devuelve la relacion entre el menor y el mayor valor de
// la matriz ppal. esto se hace antes de aplicar condiciones de borde
{
    double mayor,menor,actual;

    long a,b;
    MBmatriz=matriz;
    MBfilas=orden;
    MBanchobanda=anchobanda;

    menor = 1E+50;

    mayor = 0;

    for (a = 0;a<orden;a++)
        {
        cbGambas(a*100/orden);
        for (b = 0;b<anchobanda;b++)
        //for (b = 0;b<orden;b++)
            {
                //actual = debanda(a, b);
                actual= matriz[a*anchobanda+b];
                if (actual<0) actual = -actual;

                if (actual!=0)
                {
                    if (actual > mayor)
                        mayor = actual;
                    if (actual < menor)
                        menor = actual;
                }
            }
        }
    if (menor==0)
        return 0;

    return mayor/menor;
}


long M8JacobiBanda(double *matriz,long anchobanda, double *x, double *ti,double *x1, long orden, double precision)
{
    /* intenta solucionar el sistema por el metodo de jacobi

               n+1      -1                    n
              X      = D   . [ B - (L + U) x X  ]

      Donde     D = matriz diagonal
                B = terminos independientes
                L = triangulo inferior de la matriz
                U = triangulo superior de la matriz
                X = incognitas
                n = numero de iteracion

      Voy a empezar por iterar con la inversa de los valores de la matriz ppal

    */
    long i,j,iteraciones;
    double vi, relacion;

    // variables de la matriz banda
    MBmatriz=matriz;
    MBfilas=orden;
    MBanchobanda=anchobanda;

    for (i=0;i<orden;i++){

        x[i]=1/debanda(i,i);        // valores iniciales de iteracion
        // d[i]=x[i];               // matriz diagonal
        //abanda (i,i,0);           // armo L+U  que se guarda en abanda
    }

    //return anchobanda;

    for (iteraciones=0;iteraciones<1000;iteraciones++){
    // o sea como maximo haremos 1000 iteraciones

        vi = x[MBfilas/2];
        // esto lo saque de vector x matriz
        for (i = 0; i<orden;i++){
            x1[i] = 0;
            for (j = 0;j<orden;j++)
                if (i != j) x1[i] += x[j] * debanda(i,j);  // X1 = (L+U).Xn
        }
        // ok, se podria optimizar no operando sobre la banda

        //return -1;

        for (j = 0;j<orden;j++){
                x[j] = ti[j] - x1[j];           // X = B-(L+U).Xn
        }

        for (j = 0;j<orden;j++){
                x1[j] = x[j]/debanda(j,j);      // X1 = D-1 .(B-(L+U).Xn)
        }

        for (j = 0;j<orden;j++){
                x[j] = x1[j];                   // me preparo para la proxima
        }

        // fin de la primer iteracion, veo si me aproximo al resultado final

        relacion = x1[MBfilas/2]-vi;
        if ((relacion<precision) && (iteraciones>100))
            break;

    }// fin blucle iteraciones

    return iteraciones;

}// fin funcion


//2020 Funciones que se me perdieron
void SetMatrizRotacion(double *mr )
{
    int a;
    for (a=0; a<9; a++)
    {
        MRotacion[a] = mr[a];
    }
}

void setCG(double x, double y , double z )
{
    xCg = x;
    yCg = y;
    zCg = z;

}

void rotar3Dcg( double *punto , double *pRotado )
{

}

void SetMatrizRotacionLocal(double *mr)
{
    int a;
    for (a=0; a<9; a++)
    {
        MRotacionLocal[a] = mr[a];
    }
}


void Local3D( double *punto , double *prRotado)
 {
    prRotado[0] = MRotacionLocal[0] * punto[0] + MRotacionLocal[ 1] * punto[1]+ MRotacionLocal[ 2] * punto[2];

    prRotado[1] = MRotacionLocal[3] * punto[0] + MRotacionLocal[4] * punto[1] + MRotacionLocal[5] * punto[2];

    prRotado[2]= MRotacionLocal[6] * punto[0] + MRotacionLocal[7] * punto[1] + MRotacionLocal[8] * punto[2];
}
