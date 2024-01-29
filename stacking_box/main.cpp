#include <iostream>
#include <string>
#include <cstring>
#include <cmath>

/*
  fccanalyze - code to analyze defects in an FCC lattice.
  Code originally taken from graphitelayeranalyze and modified
  from it.

  Also works for BCC if started as bccanalyze or with -bcc option
  or DIA or HCP...

  Original author: Kai Nordlund

  Modified by: Xin Jin
                to only calculate the minimum x and y lattice positions
*/

typedef int logical;

/* Debug defines */

logical debug;

//The '#' in '#var' is called stringizing operator
#define DEBUG(var) if (debug) fprintf(stderr, #var " = %g\n",(double) (var)); fflush(stdout); fflush(stderr);
#define DEBUGS(str) if (debug) fprintf(stderr,str); fflush(stderr); fflush(stdout)
#define DEBUGSR(str,var) if (debug) fprintf(stderr, "%s " #var " = %g\n",str,(double) (var));fflush(stdout); fflush(stderr);
#define DEBUGSS(str,str2) if (debug) fprintf(stderr, "%s " #str2 " = %s\n",str,str2);fflush(stdout); fflush(stderr);
#define DEBUGSRR(str,var,var2) if (debug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g\n",str,(double) (var),(double) (var2));fflush(stdout); fflush(stderr);
#define DEBUGSRRR(str,var,var2,var3) if (debug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g " #var3 " = %g\n",str,(double) (var),(double) (var2),(double) (var3));fflush(stdout); fflush(stderr);
#define DEBUGSRRRR(str,var,var2,var3,var4) if (debug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g " #var3 " = %g " #var4 " = %g\n",str,(double) (var),(double) (var2),(double) (var3),(double) (var4));fflush(stdout); fflush(stderr);
#define DEBUGSRRRRR(str,var,var2,var3,var4,var5) if (debug) fprintf(stderr, "%s " #var " = %g " #var2 " = %g " #var3 " = %g " #var4 " = %g " #var5 " = %g\n",str,(double) (var),(double) (var2),(double) (var3),(double) (var4),(double) (var5));fflush(stdout); fflush(stderr);

/* Other defines */

int MAXAT;
#define MAXDEF 300000
#define MAXVACMOV 20000

#define WSCELLMAX 16

#define True 1
#define False 0
#define pi 3.1415926535897

#define MIN(a,b)  (((a)<(b)? (a) : (b)))
#define MAX(a,b)  (((a)>(b)? (a) : (b)))

#define MIN3(a,b,c)  ( ( (a)<(b)? (a) : (b) ) < (c)? ( (a)<(b)? (a) : (b) ) : (c) )
#define MAX3(a,b,c)  ( ( (a)>(b)? (a) : (b) ) > (c)? ( (a)>(b)? (a) : (b) ) : (c) )

/*
  Function protypes
*/

double dist(double, double, double, double, double, double );
double distmax(double, double, double, double, double, double, double, double );

void analyzedata(int nat,
                 double *unitcella,
                 int *unitcelln,
                 logical printstat,
                 logical weakPeakCondition,
                 double time,
                 double forceminx,
                 double forceminy,
                 double forceminz,
                 logical lasttime);



struct point {
    double x;
    double y;
    double z;
};

logical details=False;  /* Shall we give details of the analysis as well ? */

double xsize,ysize,zsize;
double *x,*y,*z;
int *in,*itype;
logical hasrealindices;

logical pbcz=True;

logical alloy=False;
logical nial=False,cu3au=False,gaas=False,gaasn=False;
logical interf=False;
struct point interfv;
char name[20][5];

logical movinganalysis=False;
logical vicinityanalysis=True;

int lattice;
#define FCC 1
#define BCC 2
#define DIA 3
#define HCP 4

int main(int argc, char **argv) {
    int iarg;
    int i,ii,nat;
    int r1,r2,r3,r4;

    FILE *fp;
    char file[50];
    int nline;

    int stringcol,xdatacol;
    char string[80];

    logical usetime,newtime,firstdatafound;
    logical iscoordline;
    logical printstat=False,discardout=False,disc;
    logical weakPeakCondition=False; //AbJ
    int ntime,natexpectedp,natexpected;
    int ndiscard,ndiscardx,ndiscardy,ndiscardz;
    int ndiscardfar,ndiscardxfar,ndiscardyfar,ndiscardzfar;

    char buf[160],buforig[160],c1;
    char *argbuf,arg[40][40];
    int narg;

    double casctime,prevtime;
    double forceminx=-1e30,forceminy=-1e30,forceminz=-1e30;

    double newxsize,newysize,newzsize;

    FILE *fpdisc;
    char discfile[80];

    double unitcella[4];
    int unitcelln[4];

    printf("Code run as %s\n", argv[0]);

    /*
      Equilibrium unit cell size, many other variables are defined from
      this
   */

    if ((argc>1 && strcmp(argv[1],"-h")==0.0) || argc < 2) {
        printf("Usage: [fcc/bcc/dia/hcp]analyze [options] file [stringcol [string [xdatacol [xsize [ysize [zsize]]]]]]\n");
        printf("\n");
        printf("Options: -debug           debug\n");
        printf("         -d               give detailed information\n");
        printf("         -p               print out coordinate distribution at different times\n");
        printf("         -mov             Do vacancy and interstitial moving analysis\n");
        printf("         -novicinity      Do not carry out defect vicinity analyses (may speed up things a lot)\n");
        printf("         -discardout      Discard atoms outside cell\n");
        printf("         -zsurface        Treat z as free surface\n");
        printf(" \n");
        printf("         -fcc -bcc -dia -hcp   Select lattice to analyze\n");
        printf(" \n");
        printf("         -alloy           Treat material as an alloy\n");
        printf("         -nial            Treat material as a NiAl-like alloy\n");
        printf("         -cu3au           Treat material as a Cu3Au-like alloy\n");
        printf("         -gaas            Treat material as a GaAs-like alloy\n");
        printf("         -gaasn           Treat material as a GaAsN-like alloy\n");
        printf("         -interf h j k    Treat material as an interface alloy\n");
        printf("         -a a             Unit cell a\n");
        printf("         -a3 a b c        Unit cell a b c in three dimensions\n");
        printf("         -n ncells        Number of unit cells, gives a\n");
        printf("         -n3 n1 n2 n3     Number of unit cells in 3D, gives a\n");
        printf("         -forcemin x y z  Force FCC minimum lattice position to be x y z\n");
        printf("         -line theta fii  Draw defect distributions with respect to line going\n");
        printf("                          through origo at angle theta fii\n");
        printf("         -weakP           Use weaker criteria for finding lattice peaks\n");
        printf("\n");
        printf("The target should be centered at (0, 0, 0).\n");
        printf("\n");
        printf("Data is taken as x y z beginning from column xdatacol, from lines in which\n");
        printf("stringcol contains string (element name).\n");
        printf("If stringcol==0, all lines are assumed to have xyz data\n");
        printf("If there are data at the column of xdatacol+4, the data will\n");
        printf("be read as the index of atoms (default: 0)\n");
        printf("\n");
        printf("If run as fccanalyzed, debug information given.\n");
        printf("\n");
        printf("Analyses data for defects in a perfect FCC/BCC/DIA structure.\n");
        printf("If string is not empty, takes the time step from any line not containing string but\n");
        printf("containing the string \'fs\' and analyses separately for each time. If string is empty\n");
        printf("any line with less than three arguments is interpreted to separate times.\n");
        printf("\n");
        printf("HCP analysis done for orthogonal unit cell of size (a,sqrt(3)a,c)\n");
        printf("\n");
        printf("xsize, ysize and zsize give the box size\n");
        printf("However, if there is a line containing \"boxsize\" in the input file,\n");
        printf("the above sizes may be over-written.\n");
        printf("\n");
        printf("Default: 1 Au 2 114.24 114.24 114.24\n");
        printf("\n");
        return(0);
    }

    lattice=FCC;
    if (strncmp(argv[0],"bccanalyze",10)==0) lattice=BCC;
    if (strncmp(argv[0],"diaanalyze",10)==0) lattice=DIA;
    if (strncmp(argv[0],"hcpanalyze",10)==0) lattice=HCP;
    if (strncmp(argv[0],"./bccanalyze",12)==0) lattice=BCC;
    if (strncmp(argv[0],"./diaanalyze",12)==0) lattice=DIA;
    if (strncmp(argv[0],"./hcpanalyze",12)==0) lattice=HCP;

    if (strcmp(argv[0],"fccanalyzed")==0) debug=True;
    if (strcmp(argv[0],"bccanalyzed")==0) debug=True;
    if (strcmp(argv[0],"diaanalyzed")==0) debug=True;
    if (strcmp(argv[0],"hcpanalyzed")==0) debug=True;

    iarg = 1;

    /* Set some basic variables which may be modified by command line  arguments */

    unitcella[1] = unitcella[2] = unitcella[3] = 4.08;
    unitcelln[1] = unitcelln[2] = unitcelln[3] = 0;

    stringcol=1;
    sprintf(string,"XX");
    xdatacol=2;

    xsize=114.24;
    ysize=114.24;
    zsize=114.24;

    newxsize=1e30;

    /***
     * Read the command line flags.
     *
     * strncmp compares at most the first n bytes of str1 and str2, return 0 if str1 == str2;
     * sscanf reads formatted input from a string;
     */
    while (strncmp(argv[iarg],"-",1)==0) {
        if (strcmp(argv[iarg],"-d")==0) details=True;
        else if (strcmp(argv[iarg],"-debug")==0) debug=True;

        else if (strcmp(argv[iarg],"-weakP")==0) weakPeakCondition=True; //AbJ

        else if (strcmp(argv[iarg],"-p")==0) printstat=True;
        else if (strcmp(argv[iarg],"-mov")==0) movinganalysis=True;
        else if (strcmp(argv[iarg],"-novicinity")==0) vicinityanalysis=False;
        else if (strcmp(argv[iarg],"-discardout")==0) discardout=True;
        else if (strcmp(argv[iarg],"-zsurface")==0) pbcz=False;

        else if (strcmp(argv[iarg],"-a")==0) {
            if (iarg + 1 >= argc) {printf("Error: missing parameter(s) in \"-a\".\n"); exit(1);} //Abj
            sscanf(argv[++iarg],"%lg",&(unitcella[1]));
            unitcella[3]=unitcella[2]=unitcella[1];
        }
        else if (strcmp(argv[iarg],"-a3")==0) {
            if (iarg + 3 >= argc) {printf("Error: missing parameter(s) in \"-a3\".\n"); exit(1);} //Abj
            sscanf(argv[++iarg],"%lg",&(unitcella[1]));
            sscanf(argv[++iarg],"%lg",&(unitcella[2]));
            sscanf(argv[++iarg],"%lg",&(unitcella[3]));
        }

        else if (strcmp(argv[iarg],"-n")==0) {
            if (iarg + 1 >= argc) {printf("Error: missing parameter(s) in \"-n\".\n"); exit(1);} //Abj
            sscanf(argv[++iarg],"%d",&unitcelln[1]);
            unitcelln[3]=unitcelln[2]=unitcelln[1];
        }
        else if (strcmp(argv[iarg],"-n3")==0) {
            if (iarg + 3 >= argc) {printf("Error: missing parameter(s) in \"-n3\".\n"); exit(1);} //Abj
            sscanf(argv[++iarg],"%d",&unitcelln[1]);
            sscanf(argv[++iarg],"%d",&unitcelln[2]);
            sscanf(argv[++iarg],"%d",&unitcelln[3]);
        }

        else if (strcmp(argv[iarg],"-fcc")==0) lattice=FCC;
        else if (strcmp(argv[iarg],"-bcc")==0) lattice=BCC;
        else if (strcmp(argv[iarg],"-dia")==0) lattice=DIA;
        else if (strcmp(argv[iarg],"-hcp")==0) lattice=HCP;

        else if (strcmp(argv[iarg],"-alloy")==0) alloy=True;
        else if (strcmp(argv[iarg],"-nial")==0) {
            alloy=True; nial=True; cu3au=False;
        }
        else if (strcmp(argv[iarg],"-cu3au")==0) {
            alloy=True; nial=False; cu3au=True;
        }
        else if (strcmp(argv[iarg],"-gaas")==0) {
            alloy=True; gaas=True;
        }
        else if (strcmp(argv[iarg],"-gaasn")==0) {
            alloy=True; gaasn=True;
        }
        else if (strcmp(argv[iarg],"-interf")==0) {
            if (iarg + 3 >= argc) {printf("Error: missing parameter(s) in \"-interf\".\n"); exit(1);} //Abj
            alloy=True; interf=True;
            sscanf(argv[++iarg],"%lg",&(interfv.x));
            sscanf(argv[++iarg],"%lg",&(interfv.y));
            sscanf(argv[++iarg],"%lg",&(interfv.z));
        }

        else if (strcmp(argv[iarg],"-forcemin")==0 || strcmp(argv[iarg],"-force")==0) {
            if (iarg + 3 >= argc) {printf("Error: missing parameter(s) in \"-forcemin (-force)\".\n"); exit(1);} //Abj
            sscanf(argv[++iarg],"%lg",&forceminx);
            sscanf(argv[++iarg],"%lg",&forceminy);
            sscanf(argv[++iarg],"%lg",&forceminz);
        }
        else {
            fprintf(stderr,"Unknown option %s\n",argv[iarg]);
            exit(0);
        }
        iarg++;

        if (iarg >= argc) break; //Abj
    }

    if (lattice==FCC) {
        printf("Analyzing FCC lattice\n");
    }
    else if (lattice==BCC) {
        printf("Analyzing BCC lattice\n");
    }
    else if (lattice==DIA) {
        printf("Analyzing DIAmond lattice\n");
    }
    else if (lattice==HCP) {
        printf("Analyzing HCP lattice\n");
    }

    if (iarg < argc) { //Abj
        sscanf(argv[iarg++],"%s",file);
        if (details) printf("file %s\n",file);
        if (strlen(file)==0) { printf("No filename given. Exiting\n"); return(0); }
    }


    DEBUGSR("Args",iarg);
    if (iarg < argc) sscanf(argv[iarg++],"%d",&stringcol);
    if (iarg < argc) sscanf(argv[iarg++],"%s",string);
    if (iarg < argc) sscanf(argv[iarg++],"%d",&xdatacol);
    if (iarg < argc) sscanf(argv[iarg++],"%lg",&xsize);
    ysize=xsize; zsize=xsize;
    if (iarg < argc) sscanf(argv[iarg++],"%lg",&ysize);
    if (iarg < argc) sscanf(argv[iarg++],"%lg",&zsize);

    if (unitcelln[1]>0) {
        unitcella[1]=xsize/unitcelln[1];
        unitcella[2]=ysize/unitcelln[2];
        unitcella[3]=zsize/unitcelln[3];
        printf("Number of unitcells %d %d %d => using a=%g %g %g\n",
               unitcelln[1],unitcelln[2],unitcelln[3],
               unitcella[1],unitcella[2],unitcella[3]);
    }

    printf("Program arguments %d %s %d %lg %lg %lg\n",
           stringcol,string,xdatacol,xsize,ysize,zsize);
    fflush(NULL);

    usetime=True;
    if (stringcol==0) {
        usetime=False;
        printf("Stringcol 0, assuming all atoms are of the same time\n");
    }
    if (alloy) {
        printf("Treating material as an alloy\n");
        if (nial) {
            printf("... NiAl-like alloy, to be exact\n");
            printf("NOT YET SUPPORTED, EXITING !!");
            exit(0);
        }
        if (cu3au) {
            printf("... Cu3Au-like alloy, to be exact\n");
        }
        if (gaas) {
            printf("... GaAs-like alloy, to be exact\n");
        }
        if (gaasn) {
            printf("... GaAsN-like alloy, to be exact\n");
        }
        if (interf) {
            printf("... Interface alloy with interface plane %g %g %g\n",
                   interfv.x, interfv.y,interfv.z);
            printf("where plane is centered on the origin\n");
        }
        for(i=0;i<20;i++) strncpy(name[i],"\0",1);
    }

    /***
     * Allocate heap memories for atomic coordinates and index
     */

    ntime=0;
    casctime=0.0;

    MAXAT=10000;

    i=sizeof(double)*MAXAT*3+sizeof(int)*MAXAT;
    if (alloy) i+=sizeof(int)*MAXAT;
    printf("Allocating %d bytes for %d atoms ... ",i,MAXAT);
    x=(double *) malloc(sizeof(double)*MAXAT);
    y=(double *) malloc(sizeof(double)*MAXAT);
    z=(double *) malloc(sizeof(double)*MAXAT);
    in=(int *)  malloc(sizeof(int)*MAXAT);
    if (alloy) itype=(int *) malloc(sizeof(int)*MAXAT);
    printf("... done\n");

    if (strcmp(file, "_")==0) fp=stdin;
    else fp=fopen(file,"r");

    if (fp==NULL) {
        printf("File %s open failed, exiting !\n",file);
        exit(1);
    }

    DEBUGS("File opened\n");

    if (discardout) {
        //sprintf: sends formatted output to a string pointed to str.
        sprintf(discfile,"discarded.out");
        fpdisc=fopen(discfile,"w");
        printf("Discarded atoms written to %s\n",discfile);
    }
    if (! pbcz) {
        printf("Treating z surface as open\n");
    }

    nat=0;
    nline=0;
    firstdatafound=False;
    hasrealindices=False;
    newxsize=xsize; newysize=ysize; newzsize=zsize;
    ndiscard=ndiscardx=ndiscardy=ndiscardz=0;
    ndiscardxfar=ndiscardyfar=ndiscardzfar=0;

    //fgets reads a line from the specified stream and stores it into the string pointed to by str,
    //  it stops when either (n-1) characters, the newline character or the end-of-file are read.
    while(fgets(buf,159,fp) != NULL) {
        nline++;

        /* Fix atom type from first xyz atom line if not specified above */
        /* string is the element type, which is read from the 3rd line */
        /* stringcol indicates the column where the name of element is written */
        /* if element name in a line is found at stringcol, then this line contains coordinates */
        /* if non-specified element name appears at stringcol, it will be ignored as an impurity */
        /* A line with single number gives the total number of atoms in a frame */
        /* The number in front of "fs" represents the time of a frame (same line) */
        /* The three numbers after "boxsize" represents cell sizes (same line) */
        /* in[] is the array of atom indices */
        /* The origin of coordinates should be at the target center */

        /***
         * If no element name is given at the command line,
         * the program takes the first element name at line 3 as the default target element.
         */
        if (nline==3 && strcmp(string,"XX")==0) {
            strncpy(string,buf,2);
            printf("Picked atom name string %s from line %d\n",string,nline);
        }

        /***
         * Split line into arguments by using the strtok.
         * The programs cal walk through the token,
         * and decide if a line is a coordinate line or contains impurity atom or none of the above.
        */

        //copies the string pointed to by buf to buforig
        strcpy(buforig,buf);

        if (strlen(buf)>150) {
            printf("Warning: long line > 150 : %s\n",buf); //Mbj
            continue;
        }

        narg=0;
        iscoordline=False;

        //strtok breaks string str into a series of tokens using the delimiter, and argbuf is the first token
        if ((argbuf=(char *) strtok(buf," \n\t\0")) != NULL) {
            narg++;
            strcpy(arg[narg],argbuf);

            if (narg==stringcol) {
                /* Check whether line is coordinate line */
                //Note the 2D array of char returns a int value, in order to print it, %c needs to be used instead of %s
                c1=arg[stringcol][0];
                if (strlen(arg[stringcol])<3&&c1>='A'&&c1<='Z') {
                    if (!alloy && !(strcmp(arg[stringcol],string)==0))  {
                        printf("Impurity ignored: %s (impurity element) %s (target element) %s",arg[stringcol],string,buforig);
                    }
                    else {
                        iscoordline=True;
                    }
                }
            }

            //Walks through the remaining token, which is from the 2nd one
            while((argbuf=(char *) strtok(NULL," \n\t\0")) != NULL) {
                narg++; strcpy(arg[narg],argbuf);
                if (narg==stringcol) {
                    if (strcmp(arg[stringcol],string)==0) {
                        iscoordline=True;
                    }
                    else {
                        c1=arg[stringcol][0];
                        if (c1>='A'&&c1<='Z'&&strlen(arg[stringcol])<3) {
                            printf("Impurity ignored: %s",buforig);
                        }
                    }
                }
            }
            if (stringcol==0) iscoordline=True;
        }

        /***
         * Determine if there are atom indices in the file, which should be at the column of xdatacol+4.
         * Note: it's OK to not have the atom indice, in which case there should be no more than 3 data following the x coordinate
         */
        if (nline==3 && narg>=xdatacol+4) {
            hasrealindices=True;
            printf("Taking atom indices from column %d\n",xdatacol+4);
        }

        /***
         * If the stringcol is not set to 0, usetime will be true, and if the lines are not coordinate lines:
         *      we will enter into this loop
         *
         *      Find frame time, the total number of atoms in a frame and cell sizes.
         *
         *
         *      If a line only contains one token, it is probably the line of total atom number;
         *      If a line contains "fs", then probably the time in front of "fs" is the frame time casctime;
         *      If a line contains "boxsize", then the following three numbers are the box sizes.
        */
        newtime=False;
        if (usetime) {
            if (! iscoordline) {
                DEBUGSS("Not coord line",buforig);
                /* Check whether this is a new time */
                if (strlen(string)==0 && narg < xdatacol+2) {
                    prevtime=casctime;
                    ntime++; casctime+=1.0;
                    DEBUGSRR("New time found",ntime,casctime);
                }
                if (narg==1) {
                    i=sscanf(arg[1],"%d",&ii);
                    if (i==1) {
                        natexpected=ii; //Mbj
                        natexpectedp=natexpected;
                        printf("\nExpecting to read in %d atoms for next time (line %d)\n",
                               natexpected,nline);
                    }
                }

                //strstr finds the first occurrence of the substring "fs" in buforig
                if (strlen(string) != 0 && strstr(buforig,"fs") != NULL) {
                    newtime=True;
                    prevtime=casctime;
                    ntime++; casctime+=1.0;
                    for (i=1;i<=narg;i++) {
                        if (strcmp(arg[i],"fs") == 0) sscanf(arg[i-1],"%lg",&casctime);
                    }
                }

                for (i=1;i<=narg;i++) {
                    /* Check whether line has box size as well */
                    if (strcmp(arg[i],"boxsize")==0) {
                        if (newxsize!=1e30) {
                            xsize=newxsize;
                            ysize=newysize;
                            zsize=newzsize;
                            if (unitcelln[1]>0) {
                                unitcella[1]=xsize/unitcelln[1];
                                unitcella[2]=ysize/unitcelln[2];
                                unitcella[3]=zsize/unitcelln[3];
                            }
                        }
                        sscanf(arg[i+1],"%lg",&newxsize);
                        sscanf(arg[i+2],"%lg",&newysize);
                        sscanf(arg[i+3],"%lg",&newzsize);
                        printf("Picked next cell size %lg %lg %lg from time line\n",
                               newxsize,newysize,newzsize);
                    }
                }

                DEBUGSRRR("New time found",ntime,casctime,nline);
            }
        }

        /***
         * Read in atomic coordinates and the atomic index (optional).
         *
         * If a line is a coordinate line and there are at least 2 data following the x coordinate:
         *      Enter into the loop
         *
         * Read x, y and z;
         * If there is a column of xdatacol + 4, the atom index will be read from that column, otherwise use the default 0;
         * Increase the heap memory for coordinates and index (optional).
         */

        if (!newtime) {
            /* Read in atoms */
            if (iscoordline && narg>=xdatacol+2) {
                firstdatafound=True;
                r1=r2=r3=r4=1;
                //nat represents the index of the array of x, y, z and in, r1 is the returned value of sscanf which is just 1
                r1=sscanf(arg[xdatacol+0],"%lg",&(x[nat]));
                r2=sscanf(arg[xdatacol+1],"%lg",&(y[nat]));
                r3=sscanf(arg[xdatacol+2],"%lg",&(z[nat]));
                in[nat]=0;
                if (narg>=xdatacol+4) r4=sscanf(arg[xdatacol+4],"%d",&(in[nat]));
                if (alloy) {
                    sscanf(arg[xdatacol+3],"%d",&(itype[nat]));
                    if (itype[nat]<0) itype[nat]=-itype[nat];
                    if (itype[nat]>=0&&itype[nat]<20)
                        strncpy(name[itype[nat]],arg[stringcol],5);
                }

                if (hasrealindices && narg<xdatacol+4) {
                    printf("ERROR: Atom index suddenly missing on line %d :\n   %s\n",
                           nline,buforig);
                    printf("Ignoring this line completely\n\n");
                    fflush(NULL);
                    continue;
                }

                if (r1<1||r1==EOF||r2<1||r2==EOF||r3<1||r3==EOF||r4<1||r4==EOF) {
                    printf("ERROR: invalid line %d :\n   %s\n",nline,buforig);
                    printf("Ignoring this completely\n\n");
                    fflush(NULL);
                    continue;
                }

                if (discardout) {
                    disc=False;
                    if (x[nat]<-newxsize/2 || x[nat]>newxsize/2) { disc=True; ndiscardx++; }
                    if (y[nat]<-newysize/2 || y[nat]>newysize/2) { disc=True; ndiscardy++; }
                    if (z[nat]<-newzsize/2 || z[nat]>newzsize/2) { disc=True; ndiscardz++; }
                    if (x[nat]<-newxsize/2-newxsize/10 ||
                        x[nat]>newxsize/2+newxsize/10)
                    { disc=True; ndiscardxfar++; }
                    if (y[nat]<-newysize/2-newysize/10 ||
                        y[nat]>newysize/2+newysize/10)
                    { disc=True; ndiscardyfar++; }
                    if (z[nat]<-newzsize/2-newzsize/10 ||
                        z[nat]>newzsize/2+newzsize/10)
                    { disc=True; ndiscardzfar++; }

                    if (disc) {
                        fprintf(fpdisc,"%lg %lg %lg %lg %d %d",x[nat],y[nat],z[nat],casctime,nat,in[nat]);
                        if (alloy) fprintf(fpdisc,"%d\n",itype[nat]);
                        else fprintf(fpdisc,"\n");

                        nat--;
                        ndiscard++;
                    }
                }

                nat++;
                if (nat>=MAXAT) {
                    MAXAT*=2;
                    i=sizeof(double)*MAXAT*3+sizeof(int)*MAXAT;
                    if (alloy) i+=sizeof(int)*MAXAT;
                    printf("Reallocating %d bytes for %d atoms ... ",i,MAXAT);
                    x=(double *) realloc(x,sizeof(double)*MAXAT);
                    y=(double *) realloc(y,sizeof(double)*MAXAT);
                    z=(double *) realloc(z,sizeof(double)*MAXAT);
                    in=(int *) realloc(in,sizeof(int)*MAXAT);
                    if (alloy) itype=(int *) realloc(itype,sizeof(int)*MAXAT);
                    printf("... done\n");
                }
            }
        }

        //Only works when a line contains coordinates and "fs"
        //In my opinion, the firstdatafound and newtime can only be both ture when the firstdatafound is true at line i
        // and the newtime is true at line i+1
        if (firstdatafound && newtime) {
            printf(" Read in %d atoms for time %d %lg fs\n",nat,ntime-1,prevtime);
            if (discardout) {
                printf("Discarded %d atoms outside box; x %d, y %d, z %d\n",
                       ndiscard,ndiscardx,ndiscardy,ndiscardz);
                printf("Discarded atoms farther than size/10 out: x %d, y %d, z %d\n",
                       ndiscardxfar,ndiscardyfar,ndiscardzfar);
                printf("Adatoms: %d Sputtered: %d Cutoff %g �\n",ndiscard-ndiscardfar,ndiscardfar,newzsize/10);
            }
            ndiscard=ndiscardx=ndiscardy=ndiscardz=0;
            ndiscardxfar=ndiscardyfar=ndiscardzfar=0;

            if (nat<3) {
                printf("No or too few atoms to handle, looking for next step %d\n",nat);
                continue;
            }
            if (natexpected!=0) {
                if (nat!=natexpectedp && !discardout) {
                    printf("Warning: Read in nat %d does not match predicted %d\n",
                           nat,natexpectedp);
                    /* continue; */
                }
            }

            DEBUGSRRR("Last atom",x[nat-1],y[nat-1],z[nat-1]);
            printf("\n");
            printf("--------------------------------------------------------------------\n");
            printf("--------- Going into defect analysis at time %10g ------------\n",prevtime);
            printf("--------------------------------------------------------------------\n");
            printf("\n");

            analyzedata(nat,unitcella,unitcelln,printstat,weakPeakCondition,prevtime,forceminx,forceminy,forceminz,False);

            printf("\n");
            printf("--------------------------------------------------------------------\n");
            printf("---------- Done with defect analysis at time %10g ------------\n",prevtime);
            printf("--------------------------------------------------------------------\n");
            printf("\n");

            newtime=False;
            firstdatafound=False;
            printf(" Starting reading in atoms for next time %g\n",casctime);

            /* Reset necessary variables */
            nat=0;

        }

    } /* End of readin loop */

    printf("Read in all atoms %d (number of lines) %s\n",nline,buf);

    if (discardout) {
        printf("Discarded %d atoms outside box; x %d, y %d, z %d\n",
               ndiscard,ndiscardx,ndiscardy,ndiscardz);
        printf("Discarded atoms farther than size/10 out: x %d, y %d, z %d\n",
               ndiscardxfar,ndiscardyfar,ndiscardzfar);
        ndiscardfar=ndiscardxfar+ndiscardyfar+ndiscardzfar;
        printf("Adatoms: %d Sputtered: %d Cutoff %g �\n",ndiscard-ndiscardfar,ndiscardfar,newzsize/10);
    }
    ndiscard=ndiscardx=ndiscardy=ndiscardz=0;

    if (nat<3) {
        printf("No or too few atoms to handle, exiting program %d\n",nat);
        return 0;
    }

    if (natexpected!=0) {
        if (nat!=natexpectedp && !discardout) {
            printf("Warning: Read in nat %d does not match predicted %d\n",
                   nat,natexpectedp);
        }
    }

    if (!usetime) casctime=1.0;

    if (newxsize!=1e30) {
        xsize=newxsize;
        ysize=newysize;
        zsize=newzsize;
        if (unitcelln[1]>0) {
            unitcella[1]=xsize/unitcelln[1];
            unitcella[2]=ysize/unitcelln[2];
            unitcella[3]=zsize/unitcelln[3];
        }
    }

    printf("\n");
    printf("--------------------------------------------------------------------\n");
    printf("--------- Going into defect analysis at time %10g ------------\n",casctime);
    printf("--------------------------------------------------------------------\n");
    printf("\n");

    printf("nat = %d\n", nat);
    printf("unitcella = %.5f %.5f %.5f\n", unitcella[1], unitcella[2], unitcella[3]);
    printf("unitcelln = %d %d %d\n", unitcelln[1], unitcelln[2], unitcelln[3]);
    printf("forcemin = %.5e %.5e %.5e\n", forceminx, forceminy, forceminz);

    /***
     * If no flag "-a" / "-a3", "-n" / "-n3" and "-forcemin" is used,
     * unitcella (4.080, 4.080, 4.080), unitcelln (0, 0, 0) and forceminx-y-z (-1e30, -1e30, -1e30) will just use the default values.
     */

    analyzedata(nat,unitcella,unitcelln,printstat,weakPeakCondition,casctime,forceminx,forceminy,forceminz,True);

    printf("\n");
    printf("--------------------------------------------------------------------\n");
    printf("---------- Done with defect analysis at time %10g ------------\n",casctime);
    printf("--------------------------------------------------------------------\n");
    printf("\n");

    return 0;
}

/*

  Actual analysis of data. Because of the enormous amount of data memory
  use has to be kept down, so a linkcell approach is used.

*/
void analyzedata(int nat,
                 double *unitcella,
                 int *unitcelln,
                 logical printstat,
                 logical weakPeakCondition,
                 double time,
                 double forceminx,
                 double forceminy,
                 double forceminz,
                 logical lasttime)
{

    /***
     *  nat: number of atoms;
     *  unitcella: lattice parameters of unit cells, starting from unitcella[1], default: (4.080, 4.080, 4.080);
     *  unitcelln: number of of unit cells, starting from unitcelln[1], default: (0, 0, 0);
     *  printstat: a logical variable to print out coordinate distribution at different times;
     *  time: current simulation time;
     *  forceminx-y-z: Force FCC minimum lattice positions to be x y z, default: (-1e30, -1e30, -1e30);
     *  lasttime: True.
     *
     *  - Divide atoms into voxels (3D bins), the bin size is determined by the unit cell size / NUNIT (order of magnitude of o.1 A);
     *  - Calculation of maximum (xmaxstat), minimum (xminstat) and average (xstatave) numbers of atoms in a bin;
     *  - Determination of a peak:
     *      - If the total number of atoms in a bin is 3 times larger than the average value;
     *      - If the total numbers of atoms in a bin and in the first right bin are 2 times larger than the average value;
     *      - If the total numbers of atoms in a bin and in the first and second right bins are 1.4 times larger than the average value;
     *        (In this process, periodic boundary conditions are used)
     *        (Counting from the center of a peak, its left/right boundary stops at the bin where the intensity is 0.2 times of the peak center intensity)
     *  - Creation of peak intensities (ixpeak) as a function of peak number (nxpeak);
     *  - Calculate the minimum peak (lattice) position, x0lat.     *
     *
     *  Note: There is a limit on the cell size, which is determined by UNITCELLMAX*unitcella (order of magnitude of 500 nm)
     *
     *  The function of WignerSeitzdefects() is at the end.
     *
     *  Modification by Xin:
     *  weakPeakCondition - If the original method cannot find more than one peak and weakPeakCondition is Ture,
     *                      use weaker criteria for the identification of peak.
     *                      (use different peak height and peak valley identification methods)
     */

    int i,j,k;

    double zz,t,r;

    double halfx,halfy,halfz;

    /* Parameters for doing summary statistics over all events */

    int nofVLn[20];
    int nofstrange=0,nofrecognized=0,nofnotrecognized=0;
    int ntotat=0;

    int d;   /* Index over dimensions: 1 x 2 y 3 z */

    //UNITCELLMAX: maximum number of unit cell along each dimension.
    //NUNIT: For each unit cell along a certain dimension, there will be NUNIT bins.
#define UNITCELLMAX 1010
#define NUNIT 40

    /* Arrays for doing peak distribution statistics */

    int xstat[UNITCELLMAX*NUNIT][4];
    logical isxpeak[UNITCELLMAX*NUNIT][4];
    int ixpeak[UNITCELLMAX*4][4];                   /* Peak intensities */
    double xpeak[UNITCELLMAX*4][4];                  /* Values of peak positions */
    double xoffs[4],x0lat[4],dxyzmean[4];
    int xmaxstat[4],xminstat[4];
    double xstatave[4];
    int nxpeak[4],ndxyzmean[4];

    //Adj
    int ystat[UNITCELLMAX*NUNIT];

    /* Temporary help variables which can be independent of dimension */

    double peaksum,peaksummax,wpeaksum,size,xx;
    int n,ix,iy,iz,ixmax,iymax,izmax,ndpeak,ip1,ip2;
    double maxz_pos;

    int ix_min, iz_min, iz_max;

    char statfile[80],dumbfile[80];
    FILE *fp;

    static logical firsttime=True;

    DEBUGS("Start of analyzedata\n");

    printf("Analyze data using size %g %g %g a %g %g %g\n",
           xsize,ysize,zsize,unitcella[1],unitcella[2],unitcella[3]);

    /* ixmax is the total number of bin along one dimension */
    /* xmaxstat is the largest bin containing an atom */
    /* xminstat is the minimum bin containing an atom */
    /* xstatave is the average bin number */
    /* xstat[ix][d] is the counter of atom numbers in each bin along each dimension */
    /* If atom number in a bin is larger than 1/5 of the average number, the bin is considered to be in the range of a peak */
    /* nxpeak is the total number of peaks */
    /* ixpeak[n][d] is the total number of atoms in a peak or peak intensity, where n is the number of the peak and d is the dimension */
    /* xpeak[][] is the peak position */
    /* The center of a peak should be at least 1.4 times higher than the average */
    /* dxyzmean[] is the mean distance between peaks */
    /* forceminx is to force the minimum lattice position */

    if (forceminx==-1e30) {

        if (details) printf("\n********** DIMENSION 1 ***********\n");
        /*
        First do statistics of x y z coordinate distribution to get x0lat
        */
        for (i=0;i<UNITCELLMAX*NUNIT; i++) {
            xstat[i][1]=0;
        }

        /*
        Do statistics of x y z distribution
        */
        xoffs[1] = xsize/2.0;

        ixmax = (int)(xsize/unitcella[1]*NUNIT)+1;

        if (ixmax > UNITCELLMAX*NUNIT) {
            printf("fccanalyze analyzedata() ERROR: ixmax %d > %d ",
                   ixmax,UNITCELLMAX*NUNIT);
            printf("for dimension 1; size %g a %g \n", xsize,unitcella[1]);
            exit(0);
        }

        /* Divide atoms into voxels according to their coordinates */
        for (i=0;i<nat;i++) {
            xx = x[i];
            ix = (int)((xx+xoffs[1])/unitcella[1]*NUNIT+0.5);
            xstat[ix][1]++;
        }

        if (printstat) {
            sprintf(statfile,"xstat.%lg",time);
            //if (d==2) sprintf(statfile,"ystat.%lg",time);
            printf("Writing x y z distribution statistics xyzstat.%lg\n",time);
            fp=fopen(statfile,"w");
            //write the coordinates of voxels and the number of atoms in the voxel
            for(i=0;i<ixmax;i++) {
                fprintf(fp,"%lg %d\n",i*unitcella[1]/NUNIT-xoffs[1],xstat[i][1]);
            }
            fclose(fp);
        }

        /*
        Get maximum, minimum and average value
        */
        xmaxstat[1]=0;
        xminstat[1]=1e8;
        xstatave[1]=0.0;
        for(i=0;i<ixmax;i++) {
            /* printf("%lg %d\n",i*unitcella[d]/NUNIT-xoffs,xstat[i]); */
            if (xstat[i][1]>xmaxstat[1]) xmaxstat[1]=xstat[i][1];
            if (xstat[i][1]<xminstat[1]) xminstat[1]=xstat[i][1];
            xstatave[1]=xstatave[1]+xstat[i][1];
        }
        if (ixmax !=0) xstatave[1]/=ixmax;

        if (details) {
            printf("Distr. along x, max %d ave. %g min %d\n",
                   xmaxstat[1],xstatave[1],xminstat[1]);
        }

        /*
        Recognize and analyze peaks
        */
        DEBUGSR("Recognize peak",ixmax);
        nxpeak[1]=0;
        peaksummax=0;

        //Original peak condition for x
        for(i=0;i<ixmax;i++) {
            xx=xstatave[1];
            /* We must be able to recognize peak over border ! */

            //ip1 and ip2 are the first and second right neighbour bins
            ip1=i+1;
            ip2=i+2;

            //Represents periodic boundary conditions
            if (ip1>=ixmax) ip1=ip1-ixmax;
            if (ip2>=ixmax) ip2=ip2-ixmax;

            if (xstat[i][1]>3*xx ||
                (xstat[i][1]>2*xx && xstat[ip1][1]>2*xx) ||
                (xstat[i][1]>1.4*xx&&xstat[ip1][1]>1.4*xx&&xstat[ip2][1]>1.4*xx)) {

                /* Peak, get its range and middle value */
                peaksum=xstat[i][1];
                wpeaksum=i*xstat[i][1];
                isxpeak[i][1]=True;
                for(j=i-1;j>0;j--) {
                    if (xstat[j][1]<xstatave[1]/5) break;
                    peaksum+=xstat[j][1];
                    wpeaksum+=i*xstat[j][1];
                    isxpeak[j][1]=True;
                }
                for(j=i+1;j<ixmax;j++) {
                    if (xstat[j][1]<xstatave[1]/5) break;
                    peaksum+=xstat[j][1];
                    wpeaksum+=i*xstat[j][1];
                    isxpeak[j][1]=True;
                }

                //Create the peak intensity, where nxpeak[d] is the number of peak
                ixpeak[nxpeak[1]][1]=peaksum;

                if (peaksum>peaksummax) peaksummax=peaksum;

                //Calculate the peak center position with the help of the wpeaksum (peak index weighted by the peaksum)
                if (peaksum!=0.0)
                    xpeak[nxpeak[1]][1]=unitcella[1]*(1.0*wpeaksum/peaksum)/NUNIT-xoffs[1];

                //printf("nxpeak[1] = %d (%d), xpeak = %.5f\n", nxpeak[1], i, xpeak[nxpeak[1]][1]);
                /* DEBUGSRRR("Peak",nxpeak[d],peaksum,xpeak[nxpeak[d]][d]); */

                /* Move i over the peak ! */
                //This is to move to the peak boundary
                i=j;

                nxpeak[1]++;
            }
        }

        //Go through the y dimension

        if (details) printf("\n********** DIMENSION 2 ***********\n");

        /*
            First do statistics of x y z coordinate distribution to get x0lat
            */
        for (i=0;i<UNITCELLMAX*NUNIT; i++) {
            xstat[i][2]=0;
        }

        /*
	        Do statistics of x y z distribution
	        */
        xoffs[2] = ysize/2.0;

        iymax = (int) (ysize/unitcella[2]*NUNIT) + 1;

        if (iymax > UNITCELLMAX*NUNIT) {
            printf("fccanalyze analyzedata() ERROR: iymax %d > %d ",
                   iymax,UNITCELLMAX*NUNIT);
            printf("for dimension 2: size %g a %g \n", ysize, unitcella[2]);
            exit(0);
        }

        //Find the voxel number of the minimum x plane, iy_min
        //(The minimum y lattice position will be searched in this plane.)
        ix_min = (int) ( (xpeak[0][1]+xoffs[1])/ unitcella[1]*NUNIT+0.5);

        /* Divide atoms into voxels according to their coordinates */
        for (i=0;i<nat;i++) {
            if (ix_min == (int)((x[i]+xoffs[1])/unitcella[1]*NUNIT+0.5)) {
                iy = (int)((y[i]+xoffs[2])/unitcella[2]*NUNIT+0.5);
                xstat[iy][2]++;
            }
        }

        if (printstat) {
            sprintf(statfile,"ystat.%lg",time);
            printf("Writing x y z distribution statistics xyzstat.%lg\n",time);
            fp=fopen(statfile,"w");
            //write the coordinates of voxels and the number of atoms in the voxel
            for(i=0;i<iymax;i++) {
                fprintf(fp,"%lg %d\n",i*unitcella[2]/NUNIT-xoffs[2],xstat[i][2]);
            }
            fclose(fp);
        }

        /*
	        Get maximum, minimum and average value
	        */

        xmaxstat[2]=0;
        xminstat[2]=1e8;
        xstatave[2]=0.0;
        for(i=0;i<iymax;i++) {
            if (xstat[i][2]>xmaxstat[2]) xmaxstat[2]=xstat[i][2];
            if (xstat[i][2]<xminstat[2]) xminstat[2]=xstat[i][2];
            xstatave[2]=xstatave[2]+xstat[i][2];
        }
        if (iymax !=0) xstatave[2]/=iymax;

        if (details) {
            printf("Distr. along y (in mini-x plane), max %d ave. %g min %d\n",
                   xmaxstat[2],xstatave[2],xminstat[2]);
        }

        /*
	        Recognize and analyze peaks
	        */
        DEBUGSR("Recognize peak",iymax);
        nxpeak[2]=0;
        peaksummax=0;

        //Original peak condition for y
        for(i=0;i<iymax;i++) {
            xx=xstatave[2];
            /* We must be able to recognize peak over border ! */

            //ip1 and ip2 are the first and second right neighbour bins
            ip1=i+1;
            ip2=i+2;

            //Represents periodic boundary conditions
            if (ip1>=iymax) ip1=ip1-iymax;
            if (ip2>=iymax) ip2=ip2-iymax;

            //The peak condition is a littile bit different with that for x
            if (xstat[i][2]>10*xx ||
                (xstat[i][2]>4*xx && xstat[ip1][2]>4*xx) ||
                (xstat[i][2]>2*xx && xstat[ip1][2]>2*xx && xstat[ip2][2]>2*xx)) {

                /* Peak, get its range and middle value */
                peaksum=xstat[i][2];
                wpeaksum=i*xstat[i][2];
                isxpeak[i][2]=True;
                for(j=i-1;j>0;j--) {
                    if (xstat[j][2]<xstatave[2]/5) break;
                    peaksum+=xstat[j][2];
                    wpeaksum+=i*xstat[j][2];
                    isxpeak[j][2]=True;
                }
                for(j=i+1;j<ixmax;j++) {
                    if (xstat[j][2]<xstatave[2]/5) break;
                    peaksum+=xstat[j][2];
                    wpeaksum+=i*xstat[j][2];
                    isxpeak[j][2]=True;
                }

                //Create the peak intensity, where nxpeak[d] is the number of peak
                ixpeak[nxpeak[2]][2]=peaksum;

                if (peaksum>peaksummax) peaksummax=peaksum;

                //Calculate the peak center position with the help of the wpeaksum (peak index weighted by the peaksum)
                if (peaksum!=0.0)
                    xpeak[nxpeak[2]][2]=unitcella[2]*(1.0*wpeaksum/peaksum)/NUNIT-xoffs[2];

                //printf("nxpeak[1] = %d (%d), xpeak = %.5f\n", nxpeak[1], i, xpeak[nxpeak[1]][1]);
                /* DEBUGSRRR("Peak",nxpeak[d],peaksum,xpeak[nxpeak[d]][d]); */

                /* Move i over the peak ! */
                //This is to move to the peak boundary
                i=j;

                nxpeak[2]++;
            }
        }

        //Go through the z dimension

        if (details) printf("\n********** DIMENSION 3 ***********\n");

        /*
            First do statistics of x y z coordinate distribution to get x0lat
            */
        for (i=0;i<UNITCELLMAX*NUNIT; i++) {
            xstat[i][3]=0;
        }

        /*
	        Do statistics of x y z distribution
	        */
        xoffs[3] = zsize/2.0;

        izmax = (int) (zsize/unitcella[3]*NUNIT) + 1;

        if (izmax > UNITCELLMAX*NUNIT) {
            printf("fccanalyze analyzedata() ERROR: iymax %d > %d ",
                   izmax,UNITCELLMAX*NUNIT);
            printf("for dimension 2: size %g a %g \n", zsize, unitcella[3]);
            exit(0);
        }

        /* Divide atoms into voxels according to their coordinates */
        for (i=0;i<nat;i++) {
            iz = (int)((z[i]+xoffs[3])/unitcella[3]*NUNIT+0.5);
            xstat[iz][3]++;
        }

        /*
	        Get maximum, minimum and average value
	        */
        xmaxstat[3]=0;
        xminstat[3]=1e8;
        xstatave[3]=0.0;
        for(i=0;i<izmax;i++) {
            if (xstat[i][3]>xmaxstat[3]) xmaxstat[3]=xstat[i][3];
            if (xstat[i][3]<xminstat[3]) xminstat[3]=xstat[i][3];
            xstatave[3]=xstatave[3]+xstat[i][3];
        }
        if (izmax !=0) xstatave[3]/=izmax;

        if (details) {
            printf("Distr. along z, max %d ave. %g min %d\n",
                   xmaxstat[3],xstatave[3],xminstat[3]);
        }

        /*
	        Recognize and analyze peaks
	        */
        DEBUGSR("Recognize peak",izmax);
        nxpeak[3]=0;
        peaksummax=0;

        //Original peak condition for z
        for(i=0;i<izmax;i++) {
            xx=xstatave[3];
            /* We must be able to recognize peak over border ! */

            //ip1 and ip2 are the first and second right neighbour bins
            ip1=i+1;
            ip2=i+2;

            //Represents periodic boundary conditions
            if (ip1>=izmax) ip1=ip1-izmax;
            if (ip2>=izmax) ip2=ip2-izmax;

            //The peak condition is same with that for x
            if (xstat[i][3]>3*xx ||
                (xstat[i][3]>2*xx && xstat[ip1][3]>2*xx) ||
                (xstat[i][3]>1.4*xx&&xstat[ip1][3]>1.4*xx&&xstat[ip2][3]>1.4*xx)) {

                /* Peak, get its range and middle value */
                peaksum=xstat[i][3];
                wpeaksum=i*xstat[i][3];
                isxpeak[i][3]=True;
                for(j=i-1;j>0;j--) {
                    if (xstat[j][3]<xstatave[3]/5) break;
                    peaksum+=xstat[j][3];
                    wpeaksum+=i*xstat[j][3];
                    isxpeak[j][3]=True;
                }
                for(j=i+1;j<ixmax;j++) {
                    if (xstat[j][3]<xstatave[3]/5) break;
                    peaksum+=xstat[j][3];
                    wpeaksum+=i*xstat[j][3];
                    isxpeak[j][3]=True;
                }

                //Create the peak intensity, where nxpeak[d] is the number of peak
                ixpeak[nxpeak[3]][3]=peaksum;

                if (peaksum>peaksummax) peaksummax=peaksum;

                //Calculate the peak center position with the help of the wpeaksum (peak index weighted by the peaksum)
                if (peaksum!=0.0)
                    xpeak[nxpeak[3]][3]=unitcella[3]*(1.0*wpeaksum/peaksum)/NUNIT-xoffs[3];

                //printf("nxpeak[1] = %d (%d), xpeak = %.5f\n", nxpeak[1], i, xpeak[nxpeak[1]][1]);
                /* DEBUGSRRR("Peak",nxpeak[d],peaksum,xpeak[nxpeak[d]][d]); */

                /* Move i over the peak ! */
                //This is to move to the peak boundary
                i=j;

                nxpeak[3]++;
            }
        }

        /*
        Get x0lat, i.e. starting position of lattice, as the minimum peak position
        */
        /* First estimate */
        x0lat[1]=xpeak[0][1]; //Minimum position of x plane
        x0lat[2]=xpeak[0][2]; //Minimum position of y rows in the x plane
        x0lat[3]=xpeak[0][3]; //Minimum position of z plane
        maxz_pos=xpeak[nxpeak[3]-1][3]; //Maximum position of z plane
        if (details) printf("Minimum FCC position first estimate %lg %lg %lg %lg\n",x0lat[1],x0lat[2],x0lat[3],maxz_pos);

        /* Get mean distance between peaks */
        /*
        dxyzmean[d]=0.0; ndxyzmean[d]=0;
        for(i=1;i<nxpeak[d];i++) {
            dxyzmean[d]+=xpeak[i][d]-xpeak[i-1][d];
            ndxyzmean[d]++;
        }
        if (ndxyzmean[d]>0) dxyzmean[d]/=(ndxyzmean[d]);
        else {
            printf("WARNING: Can't do Wigner-Seitz analysis, found too few peaks !\n");
            printf("d %d ndxyzmean[d] %d nxpeak[d] %d\n",d,ndxyzmean[d],nxpeak[d]);
            printf("Skipping to next time, let's hope for better luck then\n"); //Mbj
            return;
        }

        if (details) printf("Found %d peaks, dxyzmean %g, peaksummax %g\n",
                            nxpeak[d],dxyzmean[d],peaksummax);
        */

        /* Second estimate using all strong peaks */
        //Currently, I can't really get it
        /***
        xx=0; n=0;
        DEBUGSRR("x0lat",dxyzmean[d],nxpeak[d]);
        for(i=1;i<nxpeak[d];i++) {
            // Get number of peaks between this peak and lowest peak
            ndpeak=(int)((xpeak[i][d]-xpeak[0][d])/dxyzmean[d]+0.5);
            // Estimate lowest peak position from difference
            xx+=xpeak[i][d]-ndpeak*dxyzmean[d];
            n++;
        }
        x0lat[d]=xx/n;
        if (details) printf("Minimum FCC position second estimate %lg\n",x0lat[d]);
        ***/



        if (details) printf("\n********** END OF DIMENSION LOOP ***********\n");

    }/* End of if (forcemin) */

    if (forceminx!=-1e30) {
        x0lat[1]=forceminx;
        x0lat[2]=forceminy;
        x0lat[3]=forceminz;
    }

    printf("Minimum lattice position %lg %lg %lg %lg\n",x0lat[1],x0lat[2],x0lat[3],maxz_pos);

    //for (int i=0; i<nxpeak[3]; i++) {
    //    printf("%d %.5f\n", i, xpeak[i][3]);
    //}

    if (details) printf("\n********** DIMENSION 3 (Minimum z plane) ***********\n");

    for (i=0;i<UNITCELLMAX*NUNIT; i++) {
        xstat[i][1]=0;
    }

    //Find the voxel number of the minimum z plane
    //(The minimum x lattice position will be searched in this plane.)
    iz_min = (int) ( (x0lat[3]+xoffs[3])/ unitcella[3]*NUNIT+0.5);

    /* Divide atoms into voxels according to their coordinates */
    for (i=0;i<nat;i++) {
        if (iz_min == (int)((z[i]+xoffs[3])/unitcella[3]*NUNIT+0.5)) {
            ix = (int)((x[i]+xoffs[1])/unitcella[1]*NUNIT+0.5);
            xstat[ix][1]++;
        }
    }

    xmaxstat[1]=0;
    xminstat[1]=1e8;
    xstatave[1]=0.0;
    for(i=0;i<ixmax;i++) {
        if (xstat[i][1]>xmaxstat[1]) xmaxstat[1]=xstat[i][1];
        if (xstat[i][1]<xminstat[1]) xminstat[1]=xstat[i][1];
        xstatave[1]=xstatave[1]+xstat[i][1];
    }
    if (ixmax !=0) xstatave[1]/=ixmax;

    nxpeak[1]=0;
    peaksummax=0;

    //Peak condition
    for(i=0;i<ixmax;i++) {
        xx=xstatave[1];
        /* We must be able to recognize peak over border ! */

        //ip1 and ip2 are the first and second right neighbour bins
        ip1=i+1;
        ip2=i+2;

        //Represents periodic boundary conditions
        if (ip1>=ixmax) ip1=ip1-ixmax;
        if (ip2>=ixmax) ip2=ip2-ixmax;

        //The peak condition is a little bit different with that of original
        if (xstat[i][1]>10*xx ||
            (xstat[i][1]>4*xx && xstat[ip1][1]>4*xx) ||
            (xstat[i][1]>2*xx && xstat[ip1][1]>2*xx && xstat[ip2][1]>2*xx)) {

            /* Peak, get its range and middle value */
            peaksum=xstat[i][1];
            wpeaksum=i*xstat[i][1];
            isxpeak[i][1]=True;
            for(j=i-1;j>0;j--) {
                if (xstat[j][1]<xstatave[1]/5) break;
                peaksum+=xstat[j][1];
                wpeaksum+=i*xstat[j][1];
                isxpeak[j][1]=True;
            }
            for(j=i+1;j<ixmax;j++) {
                if (xstat[j][1]<xstatave[1]/5) break;
                peaksum+=xstat[j][1];
                wpeaksum+=i*xstat[j][1];
                isxpeak[j][1]=True;
            }

            //Create the peak intensity, where nxpeak[d] is the number of peak
            ixpeak[nxpeak[1]][1]=peaksum;

            if (peaksum>peaksummax) peaksummax=peaksum;

            //Calculate the peak center position with the help of the wpeaksum (peak index weighted by the peaksum)
            if (peaksum!=0.0)
                xpeak[nxpeak[1]][1]=unitcella[1]*(1.0*wpeaksum/peaksum)/NUNIT-xoffs[1];

            /* Move i over the peak ! */
            //This is to move to the peak boundary
            i=j;

            nxpeak[1]++;
        }
    }

    /* First estimate */
    x0lat[1]=xpeak[0][1]; //Minimum position of x row in the minimum z plane
    printf("3 minimum x positions on the minimum z plane: %lg %lg %lg\n",x0lat[1], xpeak[1][1], xpeak[2][1]);

    if (details) printf("\n********** DIMENSION 3 (Maximum z plane) ***********\n");

    for (i=0;i<UNITCELLMAX*NUNIT; i++) {
        xstat[i][1]=0;
    }

    //Find the voxel number of the maximum z plane
    //(The minimum x lattice position will be searched in this plane.)
    iz_max = (int) ( (maxz_pos+xoffs[3])/ unitcella[3]*NUNIT+0.5);

    /* Divide atoms into voxels according to their coordinates */
    for (i=0;i<nat;i++) {
        if (iz_max == (int)((z[i]+xoffs[3])/unitcella[3]*NUNIT+0.5)) {
            ix = (int)((x[i]+xoffs[1])/unitcella[1]*NUNIT+0.5);
            xstat[ix][1]++;
        }
    }

    xmaxstat[1]=0;
    xminstat[1]=1e8;
    xstatave[1]=0.0;
    for(i=0;i<ixmax;i++) {
        if (xstat[i][1]>xmaxstat[1]) xmaxstat[1]=xstat[i][1];
        if (xstat[i][1]<xminstat[1]) xminstat[1]=xstat[i][1];
        xstatave[1]=xstatave[1]+xstat[i][1];
    }
    if (ixmax !=0) xstatave[1]/=ixmax;

    nxpeak[1]=0;
    peaksummax=0;

    //Peak condition
    for(i=0;i<ixmax;i++) {
        xx=xstatave[1];
        /* We must be able to recognize peak over border ! */

        //ip1 and ip2 are the first and second right neighbour bins
        ip1=i+1;
        ip2=i+2;

        //Represents periodic boundary conditions
        if (ip1>=ixmax) ip1=ip1-ixmax;
        if (ip2>=ixmax) ip2=ip2-ixmax;

        //The peak condition is a little bit different with that of original
        if (xstat[i][1]>10*xx ||
            (xstat[i][1]>4*xx && xstat[ip1][1]>4*xx) ||
            (xstat[i][1]>2*xx && xstat[ip1][1]>2*xx && xstat[ip2][1]>2*xx)) {

            /* Peak, get its range and middle value */
            peaksum=xstat[i][1];
            wpeaksum=i*xstat[i][1];
            isxpeak[i][1]=True;
            for(j=i-1;j>0;j--) {
                if (xstat[j][1]<xstatave[1]/5) break;
                peaksum+=xstat[j][1];
                wpeaksum+=i*xstat[j][1];
                isxpeak[j][1]=True;
            }
            for(j=i+1;j<ixmax;j++) {
                if (xstat[j][1]<xstatave[1]/5) break;
                peaksum+=xstat[j][1];
                wpeaksum+=i*xstat[j][1];
                isxpeak[j][1]=True;
            }

            //Create the peak intensity, where nxpeak[d] is the number of peak
            ixpeak[nxpeak[1]][1]=peaksum;

            if (peaksum>peaksummax) peaksummax=peaksum;

            //Calculate the peak center position with the help of the wpeaksum (peak index weighted by the peaksum)
            if (peaksum!=0.0)
                xpeak[nxpeak[1]][1]=unitcella[1]*(1.0*wpeaksum/peaksum)/NUNIT-xoffs[1];

            /* Move i over the peak ! */
            //This is to move to the peak boundary
            i=j;

            nxpeak[1]++;
        }
    }

    /* First estimate */
    x0lat[1]=xpeak[0][1]; //Minimum position of x row in the maximum z plane
    printf("3 minimum x positions on the maximum z plane: %lg %lg %lg\n",x0lat[1], xpeak[1][1], xpeak[2][1]);

    firsttime=False;

}
