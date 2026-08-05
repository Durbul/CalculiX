/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include "CalculiX.h"

void external_solve(ITG *jq, ITG *irow, double *ad, double *au, double *b,
       ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat, double *co,
       ITG *nk, ITG *nactdof, char *jobnamef, ITG *mi, ITG *idofstart,
       ITG *ndofpernode, ITG *idofoffset){

  char jobname[133],matrixfile[160],rhsfile[160],rigfile[160],solfile[160],
    matrixarg[172],rhsarg[172],rigarg[172],solarg[172],
    solverarg[64],precondarg[64],tolarg[64],maxitersarg[64],dofsarg[64],
    line[256];
  char *path,*env,*argv[11];
  pid_t pid;
  int status,len=0,outpipe[2],nargv,haverig;
  double solvetime=-1.,*solfull=NULL;
  FILE *f;
  ITG rows,cols,i,nfull,*iperm=NULL;

  /* null-terminated copy of the blank-padded job name, used to build
     the file paths below; writetrilinos.f applies the same trimming
     rule (first blank) to derive its own .lhs/.rhs/.rig names */

  while((len<132)&&(jobnamef[len]!=' ')) len++;
  memcpy(jobname,jobnamef,len);
  jobname[len]='\0';

  /* writetrilinos.f re-expands the neq active dofs into an nfull-sized
     node-major numbering for MueLu's blocked AMG; iperm(1:neq) maps an
     active dof to its row there, used below to map the solution back
     onto the neq-sized b. It only writes the .rig (rigid-body-mode)
     file for a pure mechanical (3 dofs/node) block, so --nullspace is
     only added to the command line in that case */

  haverig=((*idofstart==1)&&(*ndofpernode==3));

  NNEW(iperm,ITG,*neq);

  FORTRAN(writetrilinos,(jq,irow,ad,au,b,neq,nzs,symmetryflag,inputformat,
			 co,nk,nactdof,jobnamef,mi,idofstart,ndofpernode,
			 idofoffset,iperm,&nfull));

  sprintf(matrixfile,"%s.lhs",jobname);
  sprintf(rhsfile,"%s.rhs",jobname);
  sprintf(rigfile,"%s.rig",jobname);
  sprintf(solfile,"%s.sol",jobname);

  sprintf(matrixarg,"--matrix=%s",matrixfile);
  sprintf(rhsarg,"--rhs=%s",rhsfile);
  sprintf(rigarg,"--nullspace=%s",rigfile);
  sprintf(solarg,"--solution=%s",solfile);

  env=getenv("CCX_EXTERNAL_SOLVER_TYPE");
  snprintf(solverarg,sizeof(solverarg),"--solver=%s",
	   env ? env : ((*symmetryflag==0) ? "cg" : "gmres"));

  env=getenv("CCX_EXTERNAL_SOLVER_PRECOND");
  snprintf(precondarg,sizeof(precondarg),"--precond=%s",env ? env : "amg");

  env=getenv("CCX_EXTERNAL_SOLVER_TOL");
  snprintf(tolarg,sizeof(tolarg),"--tol=%s",env ? env : "1e-8");

  env=getenv("CCX_EXTERNAL_SOLVER_MAXITERS");
  snprintf(maxitersarg,sizeof(maxitersarg),"--maxiters=%s",env ? env : "500");

  sprintf(dofsarg,"--dofs-per-node=%" ITGFORMAT,*ndofpernode);

  path=getenv("CCX_EXTERNAL_SOLVER");
  if(!path) path="external_solver";

  nargv=0;
  argv[nargv++]=path;
  argv[nargv++]=matrixarg;
  argv[nargv++]=rhsarg;
  if(haverig) argv[nargv++]=rigarg;
  argv[nargv++]=solarg;
  argv[nargv++]=solverarg;
  argv[nargv++]=precondarg;
  argv[nargv++]=tolarg;
  argv[nargv++]=maxitersarg;
  argv[nargv++]=dofsarg;
  argv[nargv]=NULL;

  /* pipe stdout (rather than plain fork+exec) so the parent can relay
     it and pick out the SOLVE_TIME_SECONDS line -- the external
     solver's own solve-only timing, comparable to what the in-process
     solvers report */

  if(pipe(outpipe)<0){
    printf(" *ERROR in external_solve: pipe failed\n\n");
    FORTRAN(stop,());
  }

  pid=fork();
  if(pid<0){
    printf(" *ERROR in external_solve: fork failed\n\n");
    FORTRAN(stop,());
  }else if(pid==0){
    close(outpipe[0]);
    dup2(outpipe[1],STDOUT_FILENO);
    close(outpipe[1]);

    /* Kokkos::OpenMP prints a perf warning without this set; don't
       override it if the caller already set their own (overwrite=0) */
    setenv("OMP_PROC_BIND","spread",0);

    execvp(path,argv);
    perror(" *ERROR in external_solve: exec of the external solver failed");
    _exit(127);
  }else{
    close(outpipe[1]);
    f=fdopen(outpipe[0],"r");
    while(fgets(line,256,f)!=NULL){
      fputs(line,stdout);
      sscanf(line,"SOLVE_TIME_SECONDS: %lf",&solvetime);
    }
    fclose(f);
    if((waitpid(pid,&status,0)<0)||(!WIFEXITED(status))||
       (WEXITSTATUS(status)!=0)){
      printf(" *ERROR in external_solve: the external solver failed\n\n");
      FORTRAN(stop,());
    }
    if(solvetime<0.){
      printf(" *ERROR in external_solve: the external solver did not report a solve time\n\n");
      FORTRAN(stop,());
    }
    printf(" EXTERNAL solve time: %f s\n\n",solvetime);
  }

  if((f=fopen(solfile,"r"))==NULL){
    printf(" *ERROR in external_solve: cannot open the solution file %s\n\n",
	   solfile);
    FORTRAN(stop,());
  }

  do{
    if(fgets(line,256,f)==NULL){
      printf(" *ERROR in external_solve: unexpected end of solution file %s\n\n",
	     solfile);
      FORTRAN(stop,());
    }
  }while(line[0]=='%');

  if((sscanf(line,"%" ITGFORMAT " %" ITGFORMAT,&rows,&cols)!=2)||
     (rows!=nfull)||(cols!=1)){
    printf(" *ERROR in external_solve: unexpected dimensions in solution file %s\n\n",
	   solfile);
    FORTRAN(stop,());
  }

  NNEW(solfull,double,nfull);

  for(i=0;i<nfull;i++){
    if(fscanf(f,"%lf",&solfull[i])!=1){
      printf(" *ERROR in external_solve: cannot read entry %" ITGFORMAT " from solution file %s\n\n",
	     i,solfile);
      FORTRAN(stop,());
    }
  }

  fclose(f);

  /* map the expanded-numbering solution back onto the active dofs
     (iperm is 1-based) */

  for(i=0;i<*neq;i++){
    b[i]=solfull[iperm[i]-1];
  }

  SFREE(solfull);
  SFREE(iperm);

  return;
}
